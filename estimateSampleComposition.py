#!/usr/bin/env python3.7

"""
It reads BAM-formatted files and a population VCF file with known alleles and strains. It detects reads mapped to
strain-specific alleles. You can specify various threshold values for attributes of mapped alleles, below which the
mappings are considered not informative enough, and ignored.

Usage:
    estimateSampleComposition.py [options] <bams_in>
    estimateSampleComposition.py --help

Options:
    -h --help                                       Show this help message and exit.
    -v <allele_sl_info>                             VCF file containing alleles and strain incidence info.
    -o <output_name>                                Prefix for output files [default: ].
    -b <bed_file>                                   Bed file [default: ].
    -x <x_file>                                     List of excluded samples [default: ].
    -s <sig_level>                                  Significance level for strain contributions, strains below this are elided [default: 0.0075].
    -a                                              Approximate estimation. Sometimes faster [default: True].
    -q <seq_error_rate>                             Sequencing error rate appropriate for the samples [default: 0.01].
    --min_ldsnp_coverage MIN_LDSNP_COVERAGE         Minimum proportion of covered strain-defining SNPs for strain to be viable [default: 0.0].
    --convergence_threshold CONVERGENCE_THRESHOLD   Threshold of convergence to finish iterating [default: 1e-6].
    --min_contribution MIN_CONTRIBUTION             Minimum proportion of reads supporting a LDSNP for it to be counted [default: 0.0].
    --min_reads MIN_READS                           Min coverage for region not be ignored [default: 5].
    --n_violations N_VIOLATIONS                     Number of violations in strain-defining SNPs before read is ignored. default: as implied by allele fidelity - i.e. n_violation = |log(0.001) / log(allele_err_rate)| [default: -1].
    --valid_allele_support VALID_ALLELE_SUPPORT     Min number of supporting samples for allele to be valid [default: 2].
    --debug                                         Debugging output.
    --annotated_bam                                 Create annotated BAM-formatted file as output.
"""

## LIBRARIES ##

import csv, itertools, logging, math, numpy, sys, pysam

from intervaltree import IntervalTree
from collections import OrderedDict as odict, defaultdict, namedtuple
from docopt import docopt
from os import path
from pysam import VariantFile
from sortedcontainers import SortedListWithKey as slist


## CLASSES ##

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


## METHODS ##

# read all regions of a reference from a BED file
def readBedFile(f):

    with open(f, 'r') as bed_in:

        bed_rdr = csv.reader(bed_in, delimiter='\t')
        regions = tuple((r[0], int(r[1]), int(r[2]), r[3]) for r in bed_rdr)

    return regions

def present_distribution(distr, sig_fig = 3):

    if isinstance(distr, dict):
        distr = distr.items()

    noise = 10 ** -sig_fig

    return tuple((c, round(v, sig_fig)) for c, v in sorted(distr, key=lambda kv: kv[1], reverse=True) if v >= noise)

# ignore_ref_alleles is useful when comparing two allele sets, since ref alleles are always present in both set but
# only reported where there is an alt allele also. Used in sublineage-proxy-calc.py, CollectStrainAlleles()
def readStrainAlleles(strains_vcf, valid_allele_support, exclude_samples, ignore_ref_alleles=False):

    def common_prefix_length(a, b):
        return sum(1 for _ in itertools.takewhile(lambda ab: ab[0] == ab[1], zip(a, b)))

    vcf = VariantFile(strains_vcf)
    strain_samples = {s_rec['ID']: set(s_rec['GENOMES'].split(';')) for s_rec in filter(lambda rec: rec.key == 'SAMPLE', vcf.header.records)}

    if strain_samples:  # inferred VCF; strain->samples encoded in ##SAMPLE=<ID=slin,GENOMES=s1;s2...> fields

        sample_strains = [next((sl for sl, sl_s in strain_samples.items() if s in sl_s), None) if s not in exclude_samples else None for s in vcf.header.samples]
        strain_samples = dict()
        for key, sample_strain in itertools.groupby(sorted(enumerate(sample_strains), key = lambda x : x[1]), key = lambda x : x[1]):
            strain_samples[key] = set()
            for key, value in list(sample_strain):
                strain_samples[value].add(key)

    else:   # given VCF; samples are strains

        sample_strains = [s.split('.')[0] for s in vcf.header.samples]
        strain_samples = {sl: frozenset((si,)) for si, sl in enumerate(sample_strains)}

    class AlleleOcc:

        @staticmethod
        def distr_entropy(histogram):

            entropy, Z = 0.0, 0.0

            for v in histogram:

                entropy += (math.log(v) * v) if v > 1e-18 else 0
                Z += v

            return (math.log(Z) - entropy / Z) if Z > 0 else 1.0e12

        s_strain_samples = strain_samples
        s_sample_strains = sample_strains
        s_useful_entropy = 0.999 * math.log(len(s_strain_samples))
        s_sig_sample_count = {s: min(len(s_s), valid_allele_support) for s, s_s in strain_samples.items()}

        @staticmethod
        def is_discriminating(histogram):

            return bool(histogram) if len(histogram) < len(AlleleOcc.s_strain_samples) else (AlleleOcc.distr_entropy(histogram.values()) < AlleleOcc.s_useful_entropy)

        @staticmethod
        def counts_to_frequencies(sl_counts):
            return {sl: float(slc) / len(AlleleOcc.s_strain_samples[sl]) for sl, slc in sl_counts.items() if slc >= AlleleOcc.s_sig_sample_count[sl]}

        def __init__(self, allele, strains, samples):

            assert bool(strains)

            self.strains = strains
            self.allele = allele
            self.samples = samples

        def merge_update(self, samples):

            extra = samples - self.samples

            if extra:

                for s in extra:

                    sl = AlleleOcc.s_sample_strains[s]
                    self.strains[sl] = self.strains.get(sl, 0) + 1.0 / len(AlleleOcc.s_strain_samples[sl])

                if AlleleOcc.is_discriminating(self.strains):
                    self.samples.update(extra)

                else:
                    return False

            return True

        def merge_copy(self, samples):

            extra = samples - self.samples
            smp, subl = self.samples.copy(), self.strains.copy()

            if extra:
                smp.update(extra)

                for s in extra:
                    sl = AlleleOcc.s_sample_strains[s]
                    subl[sl] = subl.get(sl, 0) + 1.0 / len(AlleleOcc.s_strain_samples[sl])

                if AlleleOcc.is_discriminating(subl):
                    smp.update(extra)

                else:
                    return None

            return AlleleOcc(self.allele, subl, smp)

    allele_incidence_per_contig = defaultdict(IntervalTree)
    sl_allele_counts = {s: 0 for s in strain_samples.keys()}
    ignored_alleles = set()

    for variant in vcf.fetch():

        allele_occs, allele_strain_counts = defaultdict(set), [{} for _ in range(len(variant.alleles))]
        missing = 0

        for sample_gt in variant.samples.itervalues():

            s_idx, allele_idx = sample_gt.index, sample_gt.allele_indices[0]  # haploid virus, has only 1 allele

            if allele_idx is not None and (not ignore_ref_alleles or allele_idx > 0):

                sl = sample_strains[s_idx]

                if sl is not None:

                    slc_a = allele_strain_counts[allele_idx]
                    slc_a[sl] = slc_a.get(sl, 0) + 1
                    allele_occs[sample_gt.alleles[0]].add(s_idx)

                else:
                    missing += 1


        allele_incidence = allele_incidence_per_contig[variant.contig]
        pos, vref, slc_r = variant.pos - 1, variant.ref, AlleleOcc.counts_to_frequencies(allele_strain_counts[0])

        if not AlleleOcc.is_discriminating(slc_r) or (pos, vref, pos + 1) in ignored_alleles:
            slc_r = None

        if len(vref) == 1 and all(len(va) == 1 for va in variant.alts):   # common case - all SNPs

            if slc_r:

                ivl = next((ivl for ivl in allele_incidence[pos] if ivl.data.allele == vref), None)

                if ivl:

                    if not ivl.data.merge_update(allele_occs[vref]):

                        allele_incidence.remove(ivl)
                        ignored_alleles.add((pos, vref, pos + 1))

                else:

                    allele_incidence.addi(pos, pos + 1, AlleleOcc(vref, slc_r, allele_occs[vref]))

                    for sl in slc_r:
                        sl_allele_counts[sl] += 1  # not exact because in rare instances we discard some alleles

            for va, slc_a in zip(variant.alts, allele_strain_counts[1:]):

                if va not in allele_occs:
                    continue

                slc_a = AlleleOcc.counts_to_frequencies(slc_a)

                if AlleleOcc.is_discriminating(slc_a):

                    allele_incidence.addi(pos, pos + 1, AlleleOcc(va, slc_a, allele_occs[va]))

                    for sl in slc_a:
                        sl_allele_counts[sl] += 1  # not exact because in rare instances we discard some alleles

    class PopulationAllelesInfo:

        def __init__(self, allele_incidence, strain_samples, sample_strains, strain_allele_counts):

            self.allele_incidence = allele_incidence
            self.strain_samples = strain_samples
            self.sample_strains = sample_strains
            self.strain_allele_counts = strain_allele_counts

    return PopulationAllelesInfo(allele_incidence=allele_incidence_per_contig,
                                 strain_samples=strain_samples,
                                 sample_strains=sample_strains,
                                 strain_allele_counts=sl_allele_counts)

def implied_errors_():

    implied_errors = namedtuple('implied_errors', ('log_odds', 'err_limit_log_odds', 'penalty_ratio', 'max_violations'))

    if not hasattr(implied_errors_, "s_values"):

        log_odds = math.log(float(args["-q"])) - math.log1p(-float(args["-q"])) - math.log(3.0)
        ierrs = implied_errors(log_odds, -math.log(1000) if int(args["--n_violations"]) < 0 else log_odds * (int(args["--n_violations"]) - 1), float(args["-q"]) / (3.0 * (1 - float(args["-q"]))), int(-math.log(1000) / log_odds) if int(args["--n_violations"]) < 0 else int(args["--n_violations"]))
        implied_errors_.s_values = ierrs

    return implied_errors_.s_values

def alignment_is_ok(align, min_mq=15, map_factor=0.75):
    return align.mapping_quality > min_mq and align.query_alignment_length > align.query_length * map_factor

# @profile
def GuesstimateInitialStrainPriors(bam_in, sl_alleles, bed_regions, min_reads, min_contrib, lin_cutoff):

    def updateAlleleContributionToStrains(intervals, base_counts):

        coverage, n_bf = numpy.sum(base_counts), 0
        n_mc = max(int(min_contrib * coverage), min_reads)

        for _, rgn in intervals:

            disc_allele = rgn.data

            if len(disc_allele.allele) == 1:

                bc = base_counts[ord(disc_allele.allele)]
                n_bf += bc

                if bc > n_mc:

                    # bc = bc / float(coverage)   # P(allele)
                    for cur_sl, cur_sl_f in disc_allele.strains.items():
                        n_c, s_c = sl_alleles_contrib[cur_sl]
                        sl_alleles_contrib[cur_sl] = (n_c + 1, s_c + bc * cur_sl_f)  # coverage, P(allele)P(sl|allele)

                if n_bf == coverage:
                    break

    def allele_locus(ctg_snp):
        return ctg_snp[0], ctg_snp[1].begin

    sl_alleles_contrib, allele_incidence = defaultdict(lambda: (0, 0)), sl_alleles.allele_incidence
    base_frequency = numpy.zeros(256, dtype=numpy.int32)
    all_snps = slist(((ct, rg) for ct, rgns in allele_incidence.items() for rg in rgns if rg.end == rg.begin + 1 and (not bed_regions or any(brg[0] == ct and brg[2] >= rg.end and brg[1] <= rg.begin for brg in bed_regions))), key=allele_locus)

    for (contig, pos), snps in itertools.groupby(all_snps, key=allele_locus):

        for locus in bam_in.pileup(contig=contig, start=pos, stop=pos + 1, truncate=True, max_depth=2000000, min_mapping_quality=30, adjust_capq_threshold=50):

            if locus.n < min_reads:

                logger.debug('ignoring position {} due to low coverage.'.format(pos))

                continue

            base_frequency.fill(0)

            for qs in locus.pileups:

                if not qs.is_refskip and alignment_is_ok(qs.alignment):
                    base_frequency[ord('*' if qs.is_del else qs.alignment.query_sequence[qs.query_position])] += 1

            updateAlleleContributionToStrains(snps, base_frequency)

    if lin_cutoff > 0:

        lc_ok = {sl: ac[0] > lin_cutoff * sl_alleles.strain_allele_counts[sl] for sl, ac in sl_alleles_contrib.items()}
        sl_priors = {sl: ac[0] if lc_ok[sl] else 0.0 for sl, ac in sl_alleles_contrib.items()}

    else:
        sl_priors = {sl: ac[0] for sl, ac in sl_alleles_contrib.items()}

    logger.debug('sl_priors before normalisation: {}'.format(sl_priors))

    slp_norm = float(sum(sl_priors.values()))

    if slp_norm > 0:
        sl_priors = {l: v / slp_norm for l, v in sl_priors.items()}

    logger.info("Initialised strain distribution: {}\n".format(present_distribution(sl_priors)))

    bam_in.seek(0)

    return sl_priors


# find likelihood of the the supplied allele set in all strains. Uses heuristic based on assumed base-call
# error rate to aggressively penalise smaller subsets
def ComputeHaplotypeLogProbabilities_approx(haplotype, all_alleles, strain_samples):

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

        for l, ls in strain_samples.items():

            l_samples, all_l_samples = ls & supporting_samples, ls & all_allele_samples

            if len(l_samples) >= min(int(args["--valid_allele_support"]), len(strain_samples[l])):  # ignore low evidence

                unsupported_alleles_penalty = 0

                for a in haplotype:

                    ss = a.data.samples

                    if not ss or not l_samples.intersection(ss):
                        unsupported_alleles_penalty += 1

                unsupported_alleles_penalty *= error_values.log_odds

                if unsupported_alleles_penalty >= error_values.err_limit_log_odds:
                    lin_probs[l] = math.log1p(len(l_samples) / (len(all_l_samples) + 1.0e-8) - 1.0) + unsupported_alleles_penalty

    # OK to limit cache size because only proximate alleles will be in cache. BAM is traversed in order, so forget past.
    if len(ComputeHaplotypeLogProbabilities_approx.lin_probs_cache) == 0x100000:

        for j in range(1024):  # evict older entries

            ComputeHaplotypeLogProbabilities_approx.lin_probs_cache.popitem(False)

    ComputeHaplotypeLogProbabilities_approx.lin_probs_cache[allele_set_str] = lin_probs

    return lin_probs


# find likelihood of the the supplied allele set in all strains. Uses heuristic based on assumed base-call
# error rate to aggressively penalise smaller subsets
# @profile
def ComputeHaplotypeLogProbabilities(haplotype, all_alleles, strain_samples, len_penalty=0.0, seen_sl_counts=None, seen_sl_samples=None):

    if not hasattr(ComputeHaplotypeLogProbabilities, 'lin_probs_cache'):

        ComputeHaplotypeLogProbabilities.lin_probs_cache = odict()

    allele_set_str = '|'.join(tuple(str((a.begin, a.end, a.data.allele)) for a in haplotype))
    lin_probs_x = ComputeHaplotypeLogProbabilities.lin_probs_cache.get(allele_set_str)

    if lin_probs_x is not None:  # could be empty map though. hence explicit comparison with None

        sl_tuple = lin_probs_x[1]

        if len(sl_tuple) > 0 and seen_sl_counts is not None:

            for l in sl_tuple:

                seen_sl_counts[l] = seen_sl_counts.get(l, 0) + 1

        return {l: v + len_penalty for l, v in lin_probs_x[0].items()}

    error_values = implied_errors_()

    # check the most common maximal subset for each strain, penalising subset score by length penalty
    lin_probs, n_alleles = {}, len(haplotype)

    if n_alleles > 1 and len_penalty >= error_values.err_limit_log_odds:

        seen_subs, seen_sub_samples = dict(), set()
        sub_penalty = len_penalty + error_values.log_odds

        # estimate of supporting samples needed for subset to win
        # req_samples = int(req_samples * ComputeAlleleSetLogProbabilities.ok_odds

        for j in range(n_alleles):  # n_alleles iterations, i.e. n!/(n-m)! complexity

            excluded_ivl = haplotype[j]

            del haplotype[j]

            exc_a = all_alleles.bisect_key_left(excluded_ivl.begin), all_alleles.bisect_key_right(excluded_ivl.begin)
            excluded_alleles = tuple(all_alleles.islice(*exc_a))

            del all_alleles[exc_a[0]:exc_a[1]]

            sub_lin_prob = ComputeHaplotypeLogProbabilities(haplotype, all_alleles, strain_samples, sub_penalty, seen_subs, seen_sub_samples).items()

            for l, v in sub_lin_prob:

                if v > lin_probs.get(l, -1000):
                    lin_probs[l] = v

            haplotype.add(excluded_ivl)

            for i, ea in enumerate(excluded_alleles):

                all_alleles.add(ea)

        seen_subs = tuple(l for l, lc in seen_subs if lc == n_alleles)

    else:
        seen_subs, seen_sub_samples = strain_samples.keys(), tuple()

    seen_sl = []
    if len(seen_subs) > 0:  # each strain consistent with a set is also consistent with all its subsets

        supporting_samples = haplotype[0].data.samples

        for ivl in haplotype[1:]:

            supporting_samples.intersection_update(ivl.data.samples)

            if not supporting_samples:
                break

        if len(supporting_samples) > 0:

            all_allele_samples = frozenset(itertools.chain.from_iterable(ivl.data.samples for ivl in all_alleles))
            error_values = implied_errors_()

            for l, ls in strain_samples.items():
                l_samples, all_l_samples = ls & supporting_samples, ls & all_allele_samples

                if len(l_samples) >= min(int(args["--valid_allele_support"]), len(strain_samples[l])):  # ignore low evidence

                    unsupported_alleles_penalty = 0

                    for a in haplotype:

                        ss = a.data.samples

                        if not ss or not l_samples.intersection(ss):
                            unsupported_alleles_penalty += 1

                    unsupported_alleles_penalty *= error_values.log_odds

                    if unsupported_alleles_penalty >= error_values.err_limit_log_odds:
                        lin_probs[l] = math.log1p(len(l_samples) / (len(all_l_samples) + 1.0e-8) - 1.0) + unsupported_alleles_penalty

    # OK to limit cache size because only proximate alleles will be in cache. BAM is traversed in order, so forget past.
    # cache needs to hold all subsets, in theory (small ones come earlier, which is not so great). So, required size,
    # assuming k alleles max, and m tolerance = k * (k-1) * ... (k-m+1). For k=30 and m=5, this is
    # 30*29*28*27*26 = 17,100,720
    if len(ComputeHaplotypeLogProbabilities.lin_probs_cache) == 0x100000:

        for j in range(1024):  # evict older entries
            ComputeHaplotypeLogProbabilities.lin_probs_cache.popitem(False)

    lin_probs_u = {l: v - len_penalty for l, v in lin_probs.items()}    # cache unpenalised values
    ComputeHaplotypeLogProbabilities.lin_probs_cache[allele_set_str] = (lin_probs_u, tuple(seen_sl))

    return lin_probs


def ComputeStrainLikelihoodsForRead(r, sl_alleles, approximate):

    contig = r.reference_name
    all_alleles = slist(sl_alleles.allele_incidence[contig][r.reference_start:r.reference_end], key=lambda i: i.begin)

    if len(all_alleles) == 0:

        raise BackboneOnlyException(r.query_name)    # skip read if it is entirely on backbone-only section

    logger.debug("considering read {}".format(r.query_name))

    have_incidence_info = any(len(sls) > 1 for sls in sl_alleles.strain_samples.values())

    haplotype, alleles_list = slist(key=all_alleles.key), slist(key=all_alleles.key)
    total_alleles, alleles_per_sl = 0, defaultdict(lambda: 0)
    seg_q_start, seg_r_start = 0, r.reference_start
    prev_allele_begin, prev_ivl = -1, None  # SNPs can overlap with start of insertions, and deletions can overlap

    for (cigar_op, cigar_len) in r.cigartuples:

        if cigar_op == 0 or cigar_op == 7 or cigar_op == 8 or cigar_op == 2:     # M, =, X, D

            seg_r_end = seg_r_start + cigar_len

            # the range below might miss some overlapping del/mnp starting before seg_r_start, but unlikely
            for ivl in all_alleles.irange_key(seg_r_start, seg_r_end, inclusive=(True, False)):

                if have_incidence_info:

                    alleles_list.add(ivl)

                q_base = '*' if cigar_op == 2 else r.query_sequence[seg_q_start + ivl.begin - seg_r_start]

                if q_base == ivl.data.allele and (cigar_op != 2 or ivl.end == seg_r_end):

                    if not (ivl.begin == prev_allele_begin and cigar_op == 2):

                        total_alleles += 1

                        for l in ivl.data.strains:
                            alleles_per_sl[l] += 1

                        prev_ivl, prev_allele_begin = ivl, ivl.begin

                        if have_incidence_info:
                            haplotype.add(ivl)

                    else:  # merge info for overlapping deletions - PITA!

                        for l in ivl.data.strains:

                            if l not in prev_ivl.data.strains:
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

        elif cigar_op == 1:   # I

            seg_q_end = seg_q_start + cigar_len

            for ivl in all_alleles.irange_key(seg_r_start, seg_r_start):

                if have_incidence_info:
                    alleles_list.add(ivl)

                a_base = ivl.data.allele,

                if a_base[0] == '&' and r.query_sequence[seg_q_start:seg_q_end] == a_base[1:]:

                    total_alleles += 1

                    for l in ivl.data.strains:
                        alleles_per_sl[l] += 1

                    prev_ivl, prev_allele_begin = ivl, ivl.begin

                    if have_incidence_info:
                        haplotype.add(ivl)

            seg_q_start = seg_q_end

        elif cigar_op == 4:     # S
            seg_q_start += cigar_len

        else:
            raise RuntimeError("cigar operation not supported")

    if have_incidence_info:

        if not bool(haplotype):
            raise BackboneOnlyException(r.query_name)

        stra_samples = sl_alleles.strain_samples
        sl_likelihoods = ComputeHaplotypeLogProbabilities_approx(haplotype, alleles_list, stra_samples) if approximate else ComputeHaplotypeLogProbabilities(haplotype, alleles_list, stra_samples)

    else:

        errs = implied_errors_()

        if approximate:  # P(alleles|l) = 1 if all alleles belong to l, else 0
            sl_likelihoods = {sl: 0 for sl, a_sl in alleles_per_sl.items() if a_sl >= total_alleles}

        else:     # P(alleles|l) = P(no-errors|l)
            min_snp_count = (total_alleles - errs.max_violations) if errs.max_violations < total_alleles else 0
            sl_likelihoods = {l: errs.log_odds * (total_alleles - sc_l) for l, sc_l in alleles_per_sl.items() if sc_l >= min_snp_count}

    logger.debug('sl_likelihoods: {}'.format(sl_likelihoods))

    if not sl_likelihoods:

        raise InconsistentSupportException(r.query_name)

    max_rw = max(sl_likelihoods.values())  # offset all values so maximum is computable (i.e. scale all likelihoods)
    offset = -(max_rw + 22) if max_rw < -22 else 0

    return {l: math.exp(sl_likelihoods[l] + offset) for l in sl_likelihoods.keys()}  # detected strains only


def ProcessSample(bam_in, regions, sl_alleles, min_reads, approximate):

    reads_range = None

    if not regions:
        regions = [(' ', 0, 3000000000, '*')]
        logger.info("no bed file specified; entire BAM is 1 region")

    elif all(r[0] == regions[0][0] for r in regions[1:]):   # save some work if all regions of interest on one contig
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

        # find region that contains the read completely, allowing for the 5bp tolerance. cse_
        rgn = next(itertools.dropwhile(lambda region: (region[2] + 5) < r.reference_end or region[1] > (r.reference_start + 5) or (region[0] != ' ' and r.reference_name != region[0]), regions), None)

        if not rgn:  # no amplicon region contains the read, even after generously stretching region by 1 bp

            ill_fitting_reads += 1
            logger.debug('read {}[{},{}] does not fit in any amplicon region.'.format(r_name, r.reference_start, r.reference_end))

            continue

        logger.debug('processing read {} - start at {}'.format(r_name, r.reference_start))

        try:
            strain_likelihoods = ComputeStrainLikelihoodsForRead(r, sl_alleles, approximate)

        except BackboneOnlyException:

            logger.debug('read {} is on backbone-only section of graph and carries no strain information'.format(r_name))

            continue

        except InconsistentSupportException:

            logger.debug('read {} has inconsistent strain support and is therefore excluded.'.format(r_name))

            continue

        ignored_reads -= 1
        p_alignment = 1 - math.pow(10, -0.1 * r.mapping_quality)  # MQ = -10 * log_10(p_err)
        reads_info.append((strain_likelihoods, p_alignment, rgn[3], r_name))
        coverage[rgn[3]] += p_alignment

        if len(reads_info) % 1000 == 0:
            logger.debug('read in {} reads'.format(len(reads_info)))

    if ignored_reads > 0:
        logger.info('ignored {} reads because they are unsuitable or uninformative, {} are not amplicon constrained.'
                      .format(ignored_reads, ill_fitting_reads))

    ignore_regions = {r: 0 for r, r_c in coverage.items() if r_c < min_reads}

    if ignore_regions:
        logger.info('also ignored {} reads in sparse regions {}.'.format(sum(coverage[r] for r in ignore_regions),
                                                                           tuple(ignore_regions.keys())))
        coverage.update(ignore_regions)

    return reads_info, coverage


def EstimateStrainPriors(reads_info, regional_coverage, sl_priors):

    n_rounds, weighted_reads = 0, [(r_l, r_p / regional_coverage[r_r]) for r_l, r_p, r_r, r_n in reads_info if regional_coverage[r_r] > 0 and any(sl_priors.get(sl, 0) > 1e-20 for sl in r_l)]

    while any(lp > 0 for lp in sl_priors.values()) and n_rounds < 2500:

        n_rounds += 1
        logger.debug('round {} - cur_lin_prevs: {}\n'.format(n_rounds, sl_priors))
        updated_priors = {l: 0.0 for l in strains}

        for r_likelihoods, r_p in weighted_reads:

            r_posteriors = {l: sl_priors[l] * s for l, s in r_likelihoods.items() if sl_priors[l] > 0}
            normalising_factor = sum(r_posteriors.values())

            if normalising_factor > 0:

                normalising_factor = r_p / normalising_factor
                updated_priors.update({l: updated_priors[l] + normalising_factor * s for l, s in r_posteriors.items()})

        normalising_factor = sum(updated_priors.values())

        if normalising_factor > 0:
            updated_priors = {l: v / normalising_factor for l, v in updated_priors.items()}

        if any(abs(v - sl_priors[l]) > float(args["--convergence_threshold"]) for l, v in updated_priors.items()):
            sl_priors = updated_priors

        else:
            break

    if n_rounds < 1000:
        logger.info("Converged to a strain priors after {} iterations.".format(n_rounds))

    else:
        logger.warning("Unable to converge to strain priors even after {} iterations!".format(n_rounds))

    return sl_priors


def FilterEstimatedPriors(strain_priors, reads_info, regional_coverage, sig_level):

    from random import choice

    strain_popularity = {l: 0 for l in strain_priors.keys()}
    t_name, t_weight, t_likelihoods = "?", 0, {}
    t_info = sorted(reads_info, key=lambda ri: ri[3])   # t for template; sorted by name, so mates come together
    t_info.append((t_likelihoods, t_weight, '', t_name))    # add a no-op footer to ensure last template gets processed

    for r_likelihoods, r_weight, r_rgn, r_name in t_info:

        if r_name != t_name:

            if r_rgn and regional_coverage[r_rgn] <= 0:
                continue

            max_p = max(t_likelihoods.get(l, 0) * lp for l, lp in strain_priors.items())
            if max_p > 0:
                sl_s = tuple(l for l, lp in strain_priors.items()
                             if lp >= sig_level and max_p - t_likelihoods.get(l, 0) * lp < 1e-12)
                if sl_s:
                    sl = sl_s[0] if len(sl_s) == 1 else choice(sl_s)  # estimated read origin
                    strain_popularity[sl] += t_weight
            t_name, t_weight, t_likelihoods = r_name, r_weight, r_likelihoods.copy()

        else:
            t_weight += r_weight
            for l, l_l in r_likelihoods.items():
                t_likelihoods[l] = t_likelihoods.get(l, 0) + l_l

    min_reads_weight = sum(strain_popularity.values()) * sig_level  # min weight for strain to be signal
    kept_strains = {sl: strain_priors[sl] for sl, w in strain_popularity.items() if w >= min_reads_weight}

    if not kept_strains:
        logger.warning("No region has total reads weight {} for any strain. Blimey!".format(min_reads_weight))

    else:
        strain_priors = kept_strains

    logger.info("Discarded these strains because they generated total read weight < {}: {}".format(min_reads_weight, present_distribution((l, p) for l, p in strain_popularity.items() if l not in kept_strains)))

    lp_sum = sum(strain_priors.values())
    sl = lp_sum * sig_level

    for l, p in strain_priors.items():

        if p < sl:

            logger.warning("Discarding popular strain {} because {} < {}".format(l, p, sl))
            strain_priors[l] = 0
            lp_sum -= p

    return present_distribution((l, p / lp_sum) for l, p in strain_priors.items())

def WriteAnnotatedOutputBam(bam_out_name, bam_in_name, reads_info, filtered_priors, all_strains):

    from copy import copy

    with pysam.AlignmentFile(bam_in_name, 'rb') as bam_in:

        header = copy(bam_in.header)    # shallow copy
        all_priors = {l: 0 for l in all_strains}
        all_priors.update({l: round(p, 3) for l, p in filtered_priors})

        comments = copy(header.get('CO', []))    # shallow copy
        comments.append("Estimated strain priors: {}".format(','.join("{}={:.3f}".format(l, p) for l, p in sorted(all_priors.items()))))
        comments.append("Strain-discriminating reads have posteriors encoded in XL tag.")
        header['CO'] = comments

        programs = copy(header.get('PG', []))    # shallow copy
        last_id = programs[-1]['ID'] if len(programs) > 0 else ''
        programs.append({'ID': last_id + '.coinf-annotation' if last_id else 'coinf-annotation',
                         'PN': sys.argv[0],
                         'Cl': ' '.join(sys.argv[1:]),
                         'PP': last_id,
                         'DS': "Discriminating reads annotated with posterior probabilities of generating strains"})
        header['PG'] = programs

        with pysam.AlignmentFile(bam_out_name, 'wb', header=header) as bam_f:

            r_idx, r_name = 0, reads_info[0][3]

            for r in bam_in.fetch():

                if r.query_name == r_name:

                    r_likelihoods = reads_info[r_idx][0]
                    norm = sum(r_likelihoods.values())

                    if norm > 0:

                        r_likelihoods = ('{}={}'.format(l, round(p / norm, 3))
                                         for l, p in sorted(r_likelihoods.items(), key=lambda _, v: -v))
                        r.set_tag('XL', ','.join(r_likelihoods))

                    r_idx += 1
                    r_name = reads_info[r_idx][3] if r_idx < len(reads_info) else ' '   # space is an invalid read name

                bam_f.write(r)

def EstimateSampleCompositions(sampleBam, bedRegions, sl_alleles, approx, sigLevel, minCoverage, minContrib, minLdsnps):

    logger.info("\nStarting analysis!")

    if not path.isfile(sampleBam + '.bai'):
        pysam.index(sampleBam)

    with pysam.AlignmentFile(sampleBam, 'rb') as bam_in:

        if minContrib < 0:
            strain_priors = {l: 1.0 / len(strains) for l in strains}

        else:
            strain_priors = GuesstimateInitialStrainPriors(bam_in, sl_alleles, bedRegions, minCoverage, minContrib, minLdsnps)
            strain_priors.update({l: 0.0 for l in strains if l not in strain_priors})

        reads_info, coverage = ProcessSample(bam_in, bedRegions, sl_alleles, minCoverage, approx)

        if not reads_info:

            logger.error("No reads to process! Please check the BED regions are compatible with {}".format(sampleBam))

            return (('?', 0),), reads_info

        strain_priors = EstimateStrainPriors(reads_info, coverage, strain_priors)
        logger.info("Estimated strain priors as {}".format(present_distribution(strain_priors)))

        filtered_priors = FilterEstimatedPriors(strain_priors, reads_info, coverage, sigLevel)
        logger.info("After filtering out noise and arranging, we have {}".format(filtered_priors))

        return filtered_priors, reads_info

## MAIN METHOD ##

if __name__ == '__main__':

    args = docopt(__doc__, version="0.6.2")

    logging.basicConfig(filename="estimateSampleComposition.log")
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG if args["--debug"] else logging.INFO)

    bams_in = args["<bams_in>"].split(',')

    if args["-x"] and path.isfile(args["-x"]):
        from re import compile

        x_samples = compile(r"\w+")
        x_samples = frozenset(w for l in open(args["-x"]) for w in x_samples.findall(l)
                              if l.strip() and not l.strip().startswith('#'))
    else:
        x_samples = frozenset(map(lambda bf: bf[:bf.index('.')], bams_in))

    sl_alleles = readStrainAlleles(args["-v"], int(args["--valid_allele_support"]), x_samples)

    logger.info("Relevant alleles = {}".format(sum(map(len, sl_alleles.allele_incidence.values()))))

    strains = frozenset(sl_alleles.strain_samples.keys())
    assert strains == frozenset(sl_alleles.strain_allele_counts.keys())

    bedRegions = readBedFile(args["-b"]) if args["-b"] else None

    def analyse_single_bam(bam_in):

        filtered_priors, reads_info = EstimateSampleCompositions(bam_in, bedRegions, sl_alleles, args["-a"],
                                                                 float(args["-s"]), int(args["--min_reads"]),
                                                                 float(args["--min_contribution"]),
                                                                 float(args["--min_ldsnp_coverage"]))
        if isinstance(args["-o"], str) and bool(reads_info):

            bam_out_name = args["-o"] + '.' + path.split(bam_in)[1].replace(".sorted.bam", ".coinfection_annotated.bam")

            if args["--annotated_bam"]:
                WriteAnnotatedOutputBam(bam_out_name, bam_in, reads_info, filtered_priors, strains)

        return bam_in, ','.join("{}={}".format(l, s) for l, s in filtered_priors)

    results = {}

    results.update(analyse_single_bam(bam_in) for bam_in in bams_in)

    with (open(args["-o"] + ".coinfection_estimates.txt", 'w') if args["-o"] else sys.stdout) as output:

        for s in bams_in:  # print results in order
            output.write('\t'.join((s, results[s] + '\n')))
