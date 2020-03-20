# !/usr/bin/env python
""" This script reads BAM files and a population VCF file with known alleles and sublineages
(e.g. created by make_population_vcf.py) and detects reads mapped to sublineage-specific alleles.
It is possible to specify various threshold values for attributes of mapped alleles,  below which
the mappings are considered not informative enough, and ignored.
"""

import argparse
import itertools
import logging
import math
import numpy
import sys
from collections import OrderedDict as odict, defaultdict, namedtuple
from os import path
import multiprocessing as concurrent  # multiprocessing.dummy is the interface for threaded concurrency, 1 core only.
from sortedcontainers import SortedListWithKey as slist

import pysam
import utils


def ParseArgs():  # one glorious day we will replace this with docopts-based parsing.
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('bams_in', nargs='+', help='input BAM files')
    arg_parser.add_argument('-v', help='VCF file containing alleles and sublineage incidence info',
                            dest='allele_sl_info', required=True)
    arg_parser.add_argument('-o', help='prefix for output files', dest='out_name', default=None)
    arg_parser.add_argument('-b', help='bed file', dest='bed_file', default=None)
    arg_parser.add_argument('-x', help='list of excluded samples', dest='x_file', default=None)
    arg_parser.add_argument('-s',
                            help='significance level for sublineage contributions (sublineages below this are elided.)',
                            default=0.0075, dest='sig_level', type=float)
    arg_parser.add_argument('-q', help='sequencing error rate appropriate for the samples',
                            default=0.01, dest='seq_error_rate', type=float)
    arg_parser.add_argument('-a', help='approximate estimation. Sometimes faster.',
                            default=False, dest='approx', action='store_true')
    arg_parser.add_argument('--min_ldsnp_coverage',
                            help='minimum proportion of covered sublineage-defining SNPs for sublineage to be viable',
                            default=0.0, dest='min_ldsnp_coverage', type=float)
    arg_parser.add_argument('--min_contribution',
                            help='minimum proportion of reads supporting a LDSNP for it to be counted', default=0.0,
                            type=float)
    arg_parser.add_argument('--min_reads', help='min coverage for region not be ignored', default=5, type=int)
    arg_parser.add_argument('--n_violations',
                            help='number of violations in sublineage-defining SNPs before read is ignored. \
    (default: as implied by allele fidelity - i.e. n_violation = |log(0.001) / log(allele_err_rate)|)',
                            default=-1, type=int)
    arg_parser.add_argument('--valid_allele_support', help='min number of supporting samples for allele to be valid',
                            default=2, type=int)
    arg_parser.add_argument('--debug', help='debugging output', default=False, dest='debug', action='store_true')
    arg_parser.add_argument('--convergence_threshold', help='threshold of convergence to finish iterating',
                            default=1e-6, dest='convergence_threshold', type=float)
    return arg_parser.parse_args()


class BackboneOnlyException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class InconsistentSupportException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def implied_errors_():
    implied_errors = namedtuple('implied_errors', ('log_odds',
                                                   'err_limit_log_odds',
                                                   'penalty_ratio',
                                                   'max_violations'))
    if not hasattr(implied_errors_, 's_values'):
        log_odds = math.log(g_args.seq_error_rate) - math.log1p(-g_args.seq_error_rate) - math.log(3.0)
        ierrs = implied_errors(log_odds,
                               -math.log(1000) if g_args.n_violations < 0 else log_odds * (g_args.n_violations - 1),
                               g_args.seq_error_rate / (3.0 * (1 - g_args.seq_error_rate)),
                               int(-math.log(1000) / log_odds) if g_args.n_violations < 0 else g_args.n_violations)
        implied_errors_.s_values = ierrs
    return implied_errors_.s_values


def alignment_is_ok(align, min_mq=15, map_factor=0.75):
    return align.mapping_quality > min_mq and align.query_alignment_length > align.query_length * map_factor


# @profile
def GuesstimateInitialSublineagePriors(bam_in, sl_alleles, bed_regions, min_reads, min_contrib, lin_cutoff):
    def updateAlleleContributionToSublineages(intervals, base_counts):
        coverage, n_bf = numpy.sum(base_counts), 0
        n_mc = max(int(min_contrib * coverage), min_reads)
        for _, rgn in intervals:
            disc_allele = rgn.data
            if len(disc_allele.allele) == 1:
                bc = base_counts[ord(disc_allele.allele)]
                n_bf += bc
                if bc > n_mc:
                    # bc = bc / float(coverage)   # P(allele)
                    for cur_sl, cur_sl_f in disc_allele.sublineages.viewitems():
                        n_c, s_c = sl_alleles_contrib[cur_sl]
                        sl_alleles_contrib[cur_sl] = (n_c + 1, s_c + bc * cur_sl_f)  # coverage, P(allele)P(sl|allele)
                if n_bf == coverage:
                    break

    def allele_locus(ctg_snp):
        return ctg_snp[0], ctg_snp[1].begin

    sl_alleles_contrib, allele_incidence = defaultdict(lambda: (0, 0)), sl_alleles.allele_incidence
    base_frequency = numpy.zeros(256, dtype=numpy.int32)
    all_snps = slist(((ct, rg) for ct, rgns in allele_incidence.viewitems() for rg in rgns
                      if rg.end == rg.begin + 1 and (not bed_regions or
                                                     any(brg[0] == ct and
                                                         brg[2] >= rg.end and
                                                         brg[1] <= rg.begin for brg in bed_regions))),
                     key=allele_locus)
    for (contig, pos), snps in itertools.groupby(all_snps, key=allele_locus):
        for locus in bam_in.pileup(contig=contig, start=pos, stop=pos + 1, truncate=True,
                                   max_depth=2000000, min_mapping_quality=30, adjust_capq_threshold=50):
            if locus.n < min_reads:
                logger.debug('ignoring position {} due to low coverage.'.format(pos))
                continue

            base_frequency.fill(0)
            for qs in locus.pileups:
                if not qs.is_refskip and alignment_is_ok(qs.alignment):
                    base_frequency[ord('*' if qs.is_del else qs.alignment.query_sequence[qs.query_position])] += 1
            updateAlleleContributionToSublineages(snps, base_frequency)

    if lin_cutoff > 0:
        lc_ok = {sl: ac[0] > lin_cutoff * sl_alleles.sublineage_allele_counts[sl]
                 for sl, ac in sl_alleles_contrib.viewitems()}
        sl_priors = {sl: ac[0] if lc_ok[sl] else 0.0 for sl, ac in sl_alleles_contrib.viewitems()}
    else:
        sl_priors = {sl: ac[0] for sl, ac in sl_alleles_contrib.viewitems()}

    logger.debug('sl_priors before normalisation: {}'.format(sl_priors))

    slp_norm = float(sum(sl_priors.viewvalues()))
    if slp_norm > 0:
        sl_priors = {l: v / slp_norm for l, v in sl_priors.viewitems()}

    logger.info("Initialised sublineage distribution: {}\n".format(utils.present_distribution(sl_priors)))

    bam_in.seek(0)
    return sl_priors


# find likelihood of the the supplied allele set in all sublineages. Uses heuristic based on assumed base-call
# error rate to aggressively penalise smaller subsets
def ComputeHaplotypeLogProbabilities_approx(haplotype, all_alleles, sublineage_samples):
    if not hasattr(ComputeHaplotypeLogProbabilities_approx, 'lin_probs_cache'):
        ComputeHaplotypeLogProbabilities_approx.lin_probs_cache = odict()

    allele_set_str = '|'.join(tuple(str((a.begin, a.end, a.data.allele)) for a in haplotype))
    lin_probs = ComputeHaplotypeLogProbabilities_approx.lin_probs_cache.get(allele_set_str)
    if lin_probs is not None:  # could be empty map though. hence explicit comparison with None
        return lin_probs

    lin_probs, supporting_samples = {}, frozenset(itertools.chain.from_iterable(a.data.samples for a in haplotype))

    if len(supporting_samples) > 0:
        all_allele_samples = frozenset(itertools.chain.from_iterable(ivl.data.samples for ivl in all_alleles))
        error_values = implied_errors_()
        for l, ls in sublineage_samples.viewitems():
            l_samples, all_l_samples = ls & supporting_samples, ls & all_allele_samples
            if len(l_samples) >= min(g_args.valid_allele_support, len(sublineage_samples[l])):  # ignore low evidence
                unsupported_alleles_penalty = 0
                for a in haplotype:
                    ss = a.data.samples
                    if not ss or not l_samples.intersection(ss):
                        unsupported_alleles_penalty += 1
                unsupported_alleles_penalty *= error_values.log_odds
                if unsupported_alleles_penalty >= error_values.err_limit_log_odds:
                    lin_probs[l] = math.log1p(len(l_samples) / (len(all_l_samples) + 1.0e-8) - 1.0) \
                                   + unsupported_alleles_penalty

    # OK to limit cache size because only proximate alleles will be in cache. BAM is traversed in order, so forget past.
    if len(ComputeHaplotypeLogProbabilities_approx.lin_probs_cache) == 0x100000:
        for j in xrange(1024):  # evict older entries
            ComputeHaplotypeLogProbabilities_approx.lin_probs_cache.popitem(False)
    ComputeHaplotypeLogProbabilities_approx.lin_probs_cache[allele_set_str] = lin_probs
    return lin_probs


# find likelihood of the the supplied allele set in all sublineages. Uses heuristic based on assumed base-call
# error rate to aggressively penalise smaller subsets
# @profile
def ComputeHaplotypeLogProbabilities(haplotype, all_alleles, sublineage_samples,
                                     len_penalty=0.0, seen_sl_counts=None, seen_sl_samples=None):
    if not hasattr(ComputeHaplotypeLogProbabilities, 'lin_probs_cache'):
        ComputeHaplotypeLogProbabilities.lin_probs_cache = odict()

    allele_set_str = '|'.join(tuple(str((a.begin, a.end, a.data.allele)) for a in haplotype))
    lin_probs_x = ComputeHaplotypeLogProbabilities.lin_probs_cache.get(allele_set_str)
    if lin_probs_x is not None:  # could be empty map though. hence explicit comparison with None
        sl_tuple = lin_probs_x[1]
        if len(sl_tuple) > 0 and seen_sl_counts is not None:
            for l in sl_tuple:
                seen_sl_counts[l] = seen_sl_counts.get(l, 0) + 1
        return {l: v + len_penalty for l, v in lin_probs_x[0].viewitems()}

    error_values = implied_errors_()

    # check the most common maximal subset for each sublineage, penalising subset score by length penalty
    lin_probs, n_alleles = {}, len(haplotype)
    if n_alleles > 1 and len_penalty >= error_values.err_limit_log_odds:
        seen_subs, seen_sub_samples = dict(), set()
        sub_penalty = len_penalty + error_values.log_odds
        # estimate of supporting samples needed for subset to win
        # req_samples = int(req_samples * ComputeAlleleSetLogProbabilities.ok_odds
        for j in xrange(n_alleles):  # n_alleles iterations, i.e. n!/(n-m)! complexity
            excluded_ivl = haplotype[j]
            del haplotype[j]
            exc_a = all_alleles.bisect_key_left(excluded_ivl.begin), all_alleles.bisect_key_right(excluded_ivl.begin)
            excluded_alleles = tuple(all_alleles.islice(*exc_a))
            del all_alleles[exc_a[0]:exc_a[1]]
            sub_lin_prob = ComputeHaplotypeLogProbabilities(haplotype, all_alleles, sublineage_samples,
                                                            sub_penalty, seen_subs, seen_sub_samples).viewitems()
            for l, v in sub_lin_prob:
                if v > lin_probs.get(l, -1000):
                    lin_probs[l] = v
            haplotype.insert(j, excluded_ivl)
            for i, ea in enumerate(excluded_alleles):
                all_alleles.insert(exc_a[0] + i, ea)
        seen_subs = tuple(l for l, lc in seen_subs if lc == n_alleles)
    else:
        seen_subs, seen_sub_samples = sublineage_samples.viewkeys(), tuple()

    seen_sl = []
    if len(seen_subs) > 0:  # each sublineage consistent with a set is also consistent with all its subsets
        supporting_samples = haplotype[0].data.samples
        for ivl in haplotype[1:]:
            supporting_samples.intersection_update(ivl.data.samples)
            if not supporting_samples:
                break

        if len(supporting_samples) >= len(seen_sub_samples) * error_values.penalty_ratio:
            if seen_sl_counts is None:  # set it to dummy object, thus avoiding checks inside the loop
                seen_sl_counts = {}
            if seen_sl_samples is not None:
                seen_sl_samples.update(supporting_samples)

            for l in seen_subs:
                l_samples = supporting_samples & sublineage_samples[l]
                all_allele_samples = frozenset(itertools.chain.from_iterable(ivl.data.samples for ivl in all_alleles))
                if len(l_samples) >= min(g_args.valid_allele_support, len(sublineage_samples[l])):  # good evidence
                    likelihood = float(len(l_samples)) / len(all_allele_samples & sublineage_samples[l])
                    seen_sl_counts[l] = seen_sl_counts.get(l, 0) + 1
                    seen_sl.append(l)
                    lp = len_penalty + math.log1p(likelihood - 1.0)
                    if lin_probs.get(l, -1000.0) < lp:
                        lin_probs[l] = lp

    # OK to limit cache size because only proximate alleles will be in cache. BAM is traversed in order, so forget past.
    # cache needs to hold all subsets, in theory (small ones come earlier, which is not so great). So, required size,
    # assuming k alleles max, and m tolerance = k * (k-1) * ... (k-m+1). For k=30 and m=5, this is
    # 30*29*28*27*26 = 17,100,720
    if len(ComputeHaplotypeLogProbabilities.lin_probs_cache) == 0x100000:
        for j in xrange(1024):  # evict older entries
            ComputeHaplotypeLogProbabilities.lin_probs_cache.popitem(False)
    lin_probs_u = {l: v - len_penalty for l, v in lin_probs.viewitems()}  # cache unpenalised values
    ComputeHaplotypeLogProbabilities.lin_probs_cache[allele_set_str] = (lin_probs_u, tuple(seen_sl))
    return lin_probs


def ComputeSublineageLikelihoodsForRead(r, sl_alleles, approximate):
    contig = r.reference_name
    all_alleles = slist(sl_alleles.allele_incidence[contig][r.reference_start:r.reference_end], key=lambda i: i.begin)
    if len(all_alleles) == 0:
        raise BackboneOnlyException(r.query_name)  # skip read if it is entirely on backbone-only section

    logger.debug("considering read {}".format(r.query_name))

    have_incidence_info = any(len(sls) > 1 for sls in sl_alleles.sublineage_samples.viewvalues())

    haplotype, alleles_list = slist(key=all_alleles.key), slist(key=all_alleles.key)
    total_alleles, alleles_per_sl = 0, defaultdict(lambda: 0)
    seg_q_start, seg_r_start = 0, r.reference_start
    prev_allele_begin, prev_ivl = -1, None  # SNPs can overlap with start of insertions, and deletions can overlap
    for (cigar_op, cigar_len) in r.cigartuples:
        if cigar_op == 0 or cigar_op == 7 or cigar_op == 8 or cigar_op == 2:  # M, =, X, D
            seg_r_end = seg_r_start + cigar_len
            # the range below might miss some overlapping del/mnp starting before seg_r_start, but unlikely
            for ivl in all_alleles.irange_key(seg_r_start, seg_r_end, inclusive=(True, False)):
                if have_incidence_info:
                    alleles_list.add(ivl)
                q_base = '*' if cigar_op == 2 else r.query_sequence[seg_q_start + ivl.begin - seg_r_start]
                if q_base == ivl.data.allele and (cigar_op != 2 or ivl.end == seg_r_end):
                    if not (ivl.begin == prev_allele_begin and cigar_op == 2):
                        total_alleles += 1
                        for l in ivl.data.sublineages:
                            alleles_per_sl[l] += 1
                        prev_ivl, prev_allele_begin = ivl, ivl.begin
                        if have_incidence_info:
                            haplotype.add(ivl)
                    else:  # merge info for overlapping deletions - PITA!
                        for l in ivl.data.sublineages:
                            if l not in prev_ivl.data.sublineages:
                                alleles_per_sl[l] += 1
                        if have_incidence_info:
                            del_ivl_idx, ao = haplotype.index(prev_ivl), prev_ivl.data.merge_copy(ivl.data.samples)
                            if ao:
                                haplotype[del_ivl_idx] = type(prev_ivl)(prev_ivl.begin, prev_ivl.end, ao)
                            else:
                                del haplotype[del_ivl_idx]
            seg_r_start = seg_r_end
            if cigar_op != 2:
                seg_q_start += cigar_len
        elif cigar_op == 1:  # I
            seg_q_end = seg_q_start + cigar_len
            for ivl in all_alleles.irange_key(seg_r_start, seg_r_start):
                if have_incidence_info:
                    alleles_list.add(ivl)
                a_base = ivl.data.allele,
                if a_base[0] == '&' and r.query_sequence[seg_q_start:seg_q_end] == a_base[1:]:
                    total_alleles += 1
                    for l in ivl.data.sublineages:
                        alleles_per_sl[l] += 1
                    prev_ivl, prev_allele_begin = ivl, ivl.begin
                    if have_incidence_info:
                        haplotype.add(ivl)
            seg_q_start = seg_q_end
        elif cigar_op == 4:  # S
            seg_q_start += cigar_len
        else:
            raise RuntimeError("cigar operation not supported")

    if have_incidence_info:
        if not bool(haplotype):
            raise BackboneOnlyException(r.query_name)

        sublin_samples = sl_alleles.sublineage_samples
        sl_likelihoods = ComputeHaplotypeLogProbabilities_approx(haplotype, alleles_list, sublin_samples) \
            if approximate else ComputeHaplotypeLogProbabilities(haplotype, alleles_list, sublin_samples)
    else:
        errs = implied_errors_()
        if approximate:  # P(alleles|l) = 1 if all alleles belong to l, else 0
            sl_likelihoods = {sl: 0 for sl, a_sl in alleles_per_sl.viewitems() if a_sl >= total_alleles}
        else:  # P(alleles|l) = P(no-errors|l)
            min_snp_count = (total_alleles - errs.max_violations) if errs.max_violations < total_alleles else 0
            sl_likelihoods = {l: errs.log_odds * (total_alleles - sc_l) for l, sc_l in alleles_per_sl.viewitems()
                              if sc_l >= min_snp_count}
    logger.debug('sl_likelihoods: {}'.format(sl_likelihoods))
    if not sl_likelihoods:
        raise InconsistentSupportException(r.query_name)

    max_rw = max(sl_likelihoods.viewvalues())  # offset all values so maximum is computable (i.e. scale all likelihoods)
    offset = -(max_rw + 22) if max_rw < -22 else 0
    return {l: math.exp(sl_likelihoods[l] + offset) for l in sl_likelihoods.viewkeys()}  # detected sublineages only


def ProcessSample(bam_in, regions, sl_alleles, min_reads, approximate):
    reads_range = None
    if not regions:
        regions = [(' ', 0, 3000000000, '*')]
        logger.info("no bed file specified; entire BAM is 1 region")
    elif all(r[0] == regions[0][0] for r in regions[1:]):  # save some work if all regions of interest on one contig
        reads_range = (regions[0][0], min(r[1] for r in regions), max(r[2] for r in regions))

    coverage = {r[3]: 0 for r in regions}

    bam_in.seek(0)
    reads_info, ignored_reads, ill_fitting_reads = [], 0, 0

    for r in bam_in.fetch(*reads_range) if reads_range else bam_in.fetch():
        r_name = r.query_name
        ignored_reads += 1  # we'll reset this if we don't end up ignoring the read
        if not alignment_is_ok(r):
            logger.warning('read {} is ignored because MQ<=15 or too many soft-clips'.format(r_name))
            continue

        # find region that contains the read completely, allowing for the 5bp tolerance.
        rgn = next(itertools.dropwhile(lambda (c, s, e, _): (e + 5) < r.reference_end or s > (r.reference_start + 5) or
                                                            (c != ' ' and r.reference_name != c),
                                       regions),
                   None)
        if not rgn:  # no amplicon region contains the read, even after generously stretching region by 1 bp
            ill_fitting_reads += 1
            logger.debug('read {}[{},{}] does not fit in any amplicon region.'.format(r_name,
                                                                                      r.reference_start,
                                                                                      r.reference_end))
            continue
        logger.debug('processing read {} - start at {}'.format(r_name, r.reference_start))
        try:
            sublineage_likelihoods = ComputeSublineageLikelihoodsForRead(r, sl_alleles, approximate)
        except BackboneOnlyException:
            logger.debug(
                'read {} is on backbone-only section of graph and carries no sublineage information'.format(r_name))
            continue
        except InconsistentSupportException:
            logger.debug('read {} has inconsistent sublineage support and is therefore excluded.'.format(r_name))
            continue

        ignored_reads -= 1
        p_alignment = 1 - math.pow(10, -0.1 * r.mapping_quality)  # MQ = -10 * log_10(p_err)
        reads_info.append((sublineage_likelihoods, p_alignment, rgn[3], r_name))
        coverage[rgn[3]] += p_alignment
        if len(reads_info) % 1000 == 0:
            logger.debug('read in {} reads'.format(len(reads_info)))

    if ignored_reads > 0:
        logger.info('ignored {} reads because they are unsuitable or uninformative, {} are not amplicon constrained.'
                    .format(ignored_reads, ill_fitting_reads))

    ignore_regions = {r: 0 for r, r_c in coverage.viewitems() if r_c < min_reads}
    if ignore_regions:
        logger.info('also ignored {} reads in sparse regions {}.'.format(sum(coverage[r] for r in ignore_regions),
                                                                         tuple(ignore_regions.viewkeys())))
        coverage.update(ignore_regions)

    return reads_info, coverage


def EstimateSublineagePriors(reads_info, regional_coverage, sl_priors):
    n_rounds, weighted_reads = 0, [(r_l, r_p / regional_coverage[r_r]) for r_l, r_p, r_r, r_n in reads_info
                                   if regional_coverage[r_r] > 0 and any(sl_priors.get(sl, 0) > 1e-20 for sl in r_l)]
    while any(lp > 0 for lp in sl_priors.viewvalues()) and n_rounds < 2500:
        n_rounds += 1
        logger.debug('round {} - cur_lin_prevs: {}\n'.format(n_rounds, sl_priors))
        updated_priors = {l: 0.0 for l in sublineages}
        for r_likelihoods, r_p in weighted_reads:
            r_posteriors = {l: sl_priors[l] * s for l, s in r_likelihoods.viewitems()
                            if sl_priors[l] > 0}
            normalising_factor = sum(r_posteriors.viewvalues())
            if normalising_factor > 0:
                normalising_factor = r_p / normalising_factor
                updated_priors.update({l: updated_priors[l] + normalising_factor * s
                                       for l, s in r_posteriors.viewitems()})

        normalising_factor = sum(updated_priors.viewvalues())
        if normalising_factor > 0:
            updated_priors = {l: v / normalising_factor for l, v in updated_priors.viewitems()}
        if any(abs(v - sl_priors[l]) > g_args.convergence_threshold for l, v in updated_priors.viewitems()):
            sl_priors = updated_priors
        else:
            break
    if n_rounds < 1000:
        logger.info("Converged to a sublineage priors after {} iterations.".format(n_rounds))
    else:
        logger.warning("Unable to converge to sublineage priors even after {} iterations!".format(n_rounds))
    return sl_priors


def FilterEstimatedPriors(sublineage_priors, reads_info, regional_coverage, sig_level):
    from random import choice

    sublineage_popularity = {l: 0 for l in sublineage_priors.viewkeys()}
    t_name, t_weight, t_likelihoods = "?", 0, {}
    t_info = sorted(reads_info, key=lambda ri: ri[3])  # t for template; sorted by name, so mates come together
    t_info.append((t_likelihoods, t_weight, '', t_name))  # add a no-op footer to ensure last template gets processed
    for r_likelihoods, r_weight, r_rgn, r_name in t_info:
        if r_name != t_name:
            if r_rgn and regional_coverage[r_rgn] <= 0:
                continue
            max_p = max(t_likelihoods.get(l, 0) * lp for l, lp in sublineage_priors.viewitems())
            if max_p > 0:
                sl_s = tuple(l for l, lp in sublineage_priors.viewitems()
                             if lp >= sig_level and max_p - t_likelihoods.get(l, 0) * lp < 1e-12)
                if sl_s:
                    sl = sl_s[0] if len(sl_s) == 1 else choice(sl_s)  # estimated read origin
                    sublineage_popularity[sl] += t_weight
            t_name, t_weight, t_likelihoods = r_name, r_weight, r_likelihoods.copy()
        else:
            t_weight += r_weight
            for l, l_l in r_likelihoods.viewitems():
                t_likelihoods[l] = t_likelihoods.get(l, 0) + l_l

    min_reads_weight = sum(sublineage_popularity.viewvalues()) * sig_level  # min weight for sublineage to be signal
    kept_sublineages = {sl: sublineage_priors[sl] for sl, w in sublineage_popularity.viewitems() if
                        w >= min_reads_weight}
    if not kept_sublineages:
        logger.warning("No region has total reads weight {} for any sublineage. Blimey!".format(min_reads_weight))
    else:
        sublineage_priors = kept_sublineages

    logger.info("Discarded these sublineages because they generated total read weight < {}: {}".format(
        min_reads_weight,
        utils.present_distribution((l, p) for l, p in sublineage_popularity.viewitems() if l not in kept_sublineages)
    ))

    lp_sum = sum(sublineage_priors.viewvalues())
    sl = lp_sum * sig_level
    for l, p in sublineage_priors.iteritems():
        if p < sl:
            logger.warning("Discarding popular sublineage {} because {} < {}".format(l, p, sl))
            sublineage_priors[l] = 0
            lp_sum -= p
    return utils.present_distribution((l, p / lp_sum) for l, p in sublineage_priors.viewitems())


def EstimateSampleCompositions(sampleBam, bedRegions, sl_alleles, approx, sigLevel, minCoverage, minContrib, minLdsnps):
    logger.info("\nStarting analysis!")
    if not path.isfile(sampleBam + '.bai'):
        pysam.index(sampleBam)

    with pysam.AlignmentFile(sampleBam, 'rb') as bam_in:
        if minContrib < 0:
            sublineage_priors = {l: 1.0 / len(sublineages) for l in sublineages}
        else:
            sublineage_priors = GuesstimateInitialSublineagePriors(bam_in, sl_alleles, bedRegions,
                                                                   minCoverage, minContrib, minLdsnps)
            sublineage_priors.update({l: 0.0 for l in sublineages if l not in sublineage_priors})

        reads_info, coverage = ProcessSample(bam_in, bedRegions, sl_alleles, minCoverage, approx)
        if not reads_info:
            logger.error("No reads to process! Please check the BED regions are compatible with {}".format(sampleBam))
            return (('?', 0),), reads_info

        sublineage_priors = EstimateSublineagePriors(reads_info, coverage, sublineage_priors)
        logger.info("Estimated sublineage priors as {}".format(utils.present_distribution(sublineage_priors)))

        filtered_priors = FilterEstimatedPriors(sublineage_priors, reads_info, coverage, sigLevel)
        logger.info("After filtering out noise and arranging, we have {}".format(filtered_priors))

        return filtered_priors, reads_info


class ConcurrentSampleLogger(logging.LoggerAdapter):

    def __init__(self, logger, extra=None):
        logging.LoggerAdapter.__init__(self, logger, extra)
        from threading import local  # this construction because I was working with threads, innocently. Leaving it in.
        self.name_cache = local()
        self.name_cache.sample = None

    def setSample(self, sample):
        self.name_cache.sample = (sample + "> ") if sample else None

    def process(self, msg, kwargs):
        if self.name_cache.sample:
            msg = self.name_cache.sample + str(msg)
        return logging.LoggerAdapter.process(self, msg, kwargs)


# ==================================================================================================================== #
if __name__ == '__main__':

    g_args = ParseArgs()

    logging.basicConfig(stream=sys.stderr)
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG if g_args.debug else logging.INFO)
    logger = ConcurrentSampleLogger(logger)

    x_samples = frozenset(map(lambda bf: bf[:bf.rindex('.')], g_args.bams_in))
    sl_alleles = utils.readSublineageAlleles(g_args.allele_sl_info, g_args.valid_allele_support, x_samples)
    logger.info("Relevant alleles = {}".format(sum(itertools.imap(len, sl_alleles.allele_incidence.viewvalues()))))

    sublineages = frozenset(sl_alleles.sublineage_samples.viewkeys())
    assert sublineages == frozenset(sl_alleles.sublineage_allele_counts.viewkeys())

    bedRegions = utils.readBedFile(g_args.bed_file) if g_args.bed_file else None


    def analyse_single_bam(bam_in):
        logger.setSample(bam_in)
        filtered_priors, reads_info = EstimateSampleCompositions(bam_in, bedRegions, sl_alleles, g_args.approx,
                                                                 g_args.sig_level, g_args.min_reads,
                                                                 g_args.min_contribution, g_args.min_ldsnp_coverage)

        return bam_in, ','.join('{}={}'.format(l, s) for l, s in filtered_priors)


    results = {}
    results.update(analyse_single_bam(b_in) for b_in in g_args.bams_in)

    with (open(g_args.out_name + '.coinfection_estimates.txt', 'w') if g_args.out_name else sys.stdout) as output:
        for s in g_args.bams_in:  # print results in order
            output.write('\t'.join((s, results[s] + '\n')))
