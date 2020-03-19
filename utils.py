import math
import csv
import itertools
from collections import defaultdict
from intervaltree import IntervalTree


def static_vars(**kwargs):
    def decorate(func):
        for k in kwargs:
            setattr(func, k, kwargs[k])
        return func
    return decorate


# read all regions of a reference from a BED file
def readBedFile(f):
    with open(f, 'r') as bed_in:
        bed_rdr = csv.reader(bed_in, delimiter='\t')
        regions = tuple((r[0], int(r[1]), int(r[2]), r[3]) for r in bed_rdr)
    return regions


# read FASTA file (returns header and sequence in separate fields)
def readFasta(f):
    with open(f, 'r') as fasta_in:
        contigs = {}
        name, sequence = '', []
        for ln in fasta_in:  # not very elegant, but itertools.takewhile() eats the element after the last matching one
            if len(ln) < 2:
                continue
            if ln[0] == '>':
                if sequence:
                    contigs[name] = ''.join(sequence)
                name, sequence = ln[1:-1].strip(), []
            elif ln[0].isalpha() or ln[0] == '-':
                sequence.append(ln[:-1])
        if sequence:
            contigs[name] = ''.join(sequence)
        return contigs


# ignore_ref_alleles is useful when comparing two allele sets, since ref alleles are always present in both set but
# only reported where there is an alt allele also. Used in sublineage-proxy-calc.py, CollectSublineageAlleles()
def readSublineageAlleles(sublineages_vcf, valid_allele_support, exclude_samples, ignore_ref_alleles=False):
    from pysam import VariantFile

    def common_prefix_length(a, b):
        return sum(1 for _ in itertools.takewhile(lambda ab: ab[0] == ab[1], itertools.izip(a, b)))

    vcf = VariantFile(sublineages_vcf)
    sublineage_samples = {s_rec['ID']: set(s_rec['GENOMES'].split(';'))
                          for s_rec in filter(lambda rec: rec.key == 'SAMPLE', vcf.header.records)}
    if sublineage_samples:  # inferred VCF; sublineage->samples encoded in ##SAMPLE=<ID=slin,GENOMES=s1;s2...> fields
        sample_sublineages = [next((sl for sl, sl_s in sublineage_samples.viewitems() if s in sl_s), None)
                              if s not in exclude_samples else None for s in vcf.header.samples]
        sublineage_samples = {sl: frozenset(i for i, _ in i_sl) for sl, i_sl in
                              itertools.groupby(sorted(enumerate(sample_sublineages), key=lambda (_, sl): sl),
                                                key=lambda (_, sl): sl) if sl}
    else:   # given VCF; samples are sublineages
        sample_sublineages = [s.split('.')[0] for s in vcf.header.samples]
        sublineage_samples = {sl: frozenset((si,)) for si, sl in enumerate(sample_sublineages)}

    class AlleleOcc:
        @staticmethod
        def distr_entropy(histogram):
            entropy, Z = 0.0, 0.0
            for v in histogram:
                entropy += (math.log(v) * v) if v > 1e-18 else 0
                Z += v
            return (math.log(Z) - entropy / Z) if Z > 0 else 1.0e12

        s_sublineage_samples = sublineage_samples
        s_sample_sublineages = sample_sublineages
        s_useful_entropy = 0.999 * math.log(len(s_sublineage_samples))
        s_sig_sample_count = {s: min(len(s_s), valid_allele_support) for s, s_s in sublineage_samples.viewitems()}

        @staticmethod
        def is_discriminating(histogram):
            return bool(histogram) if len(histogram) < len(AlleleOcc.s_sublineage_samples) else \
                (AlleleOcc.distr_entropy(histogram.viewvalues()) < AlleleOcc.s_useful_entropy)

        @staticmethod
        def counts_to_frequencies(sl_counts):
            return {sl: float(slc) / len(AlleleOcc.s_sublineage_samples[sl])
                    for sl, slc in sl_counts.viewitems() if slc >= AlleleOcc.s_sig_sample_count[sl]}

        def __init__(self, allele, sublineages, samples):
            assert bool(sublineages)
            self.sublineages = sublineages
            self.allele = allele
            self.samples = samples

        def merge_update(self, samples):
            extra = samples - self.samples
            if extra:
                for s in extra:
                    sl = AlleleOcc.s_sample_sublineages[s]
                    self.sublineages[sl] = self.sublineages.get(sl, 0) + 1.0 / len(AlleleOcc.s_sublineage_samples[sl])
                if AlleleOcc.is_discriminating(self.sublineages):
                    self.samples.update(extra)
                else:
                    return False
            return True

        def merge_copy(self, samples):
            extra = samples - self.samples
            smp, subl = self.samples.copy(), self.sublineages.copy()
            if extra:
                smp.update(extra)
                for s in extra:
                    sl = AlleleOcc.s_sample_sublineages[s]
                    subl[sl] = subl.get(sl, 0) + 1.0 / len(AlleleOcc.s_sublineage_samples[sl])
                if AlleleOcc.is_discriminating(subl):
                    smp.update(extra)
                else:
                    return None
            return AlleleOcc(self.allele, subl, smp)

    allele_incidence_per_contig = defaultdict(IntervalTree)
    sl_allele_counts = {s: 0 for s in sublineage_samples.viewkeys()}
    ignored_alleles = set()

    for variant in vcf.fetch():
        allele_occs, allele_sublineage_counts = defaultdict(set), [{} for _ in xrange(len(variant.alleles))]
        missing = 0
        for sample_gt in variant.samples.itervalues():
            s_idx, allele_idx = sample_gt.index, sample_gt.allele_indices[0]  # haploid virus, has only 1 allele
            if allele_idx is not None and (not ignore_ref_alleles or allele_idx > 0):
                sl = sample_sublineages[s_idx]
                if sl is not None:
                    slc_a = allele_sublineage_counts[allele_idx]
                    slc_a[sl] = slc_a.get(sl, 0) + 1
                    allele_occs[sample_gt.alleles[0]].add(s_idx)
                else:
                    missing += 1


        allele_incidence = allele_incidence_per_contig[variant.contig]
        pos, vref, slc_r = variant.pos - 1, variant.ref, AlleleOcc.counts_to_frequencies(allele_sublineage_counts[0])
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
            for va, slc_a in itertools.izip(variant.alts, allele_sublineage_counts[1:]):
                if va not in allele_occs:
                    continue
                slc_a = AlleleOcc.counts_to_frequencies(slc_a)
                if AlleleOcc.is_discriminating(slc_a):
                    allele_incidence.addi(pos, pos + 1, AlleleOcc(va, slc_a, allele_occs[va]))
                    for sl in slc_a:
                        sl_allele_counts[sl] += 1  # not exact because in rare instances we discard some alleles
        else:   # indel or MNP. Strip out prefix/suffix, so alt bases are where they would match up read cigars
            for va, slc_a in itertools.izip(variant.alts, allele_sublineage_counts[1:]):
                if va not in allele_occs:
                    continue
                slc_a = AlleleOcc.counts_to_frequencies(slc_a)
                if not AlleleOcc.is_discriminating(slc_a):
                    slc_a = None
                pfx_len = common_prefix_length(vref, va)    # consider vref=TAT, va=TGGGTAT
                sfx_len = common_prefix_length(reversed(vref[pfx_len:]), reversed(va[pfx_len:]))
                p_pos = pos + pfx_len
                if len(vref) == pfx_len + sfx_len:  # insertion
                    if slc_a:
                        va_t = va[pfx_len:len(va)-sfx_len]
                        assert va_t
                        allele_incidence.addi(p_pos, p_pos + 1, AlleleOcc('&'+va_t, slc_a, allele_occs[va]))
                        for sl in slc_a:
                            sl_allele_counts[sl] += 1  # not exact because in rare instances we discard some alleles
                else:
                    sref = vref[pfx_len:len(vref) - sfx_len]
                    assert sref
                    p_end = p_pos + len(sref)
                    if slc_r and (p_pos, vref, p_end) not in ignored_alleles:
                        ivl = next((ivl for ivl in
                                   allele_incidence.search(p_pos, p_end, True) if ivl.data.allele == sref), None)
                        if ivl:
                            if not ivl.data.merge_update(allele_occs[vref]):
                                allele_incidence.remove(ivl)
                                ignored_alleles.add((p_pos, vref, p_end))
                        else:
                            allele_incidence.addi(p_pos, p_end, AlleleOcc(sref, slc_r, allele_occs[vref]))
                            for sl in slc_r:
                                sl_allele_counts[sl] += 1  # not exact because in rare instances we discard some alleles
                    if slc_a:
                        sa = '*' if len(va) == pfx_len + sfx_len else va[pfx_len:len(va)-sfx_len]
                        allele_incidence.addi(p_pos, p_end, AlleleOcc(sa, slc_a, allele_occs[va]))
                        for sl in slc_a:
                            sl_allele_counts[sl] += 1  # not exact because in rare instances we discard some alleles

    class PopulationAllelesInfo:
        def __init__(self, allele_incidence, sublineage_samples, sample_sublineages, sublineage_allele_counts):
            self.allele_incidence = allele_incidence
            self.sublineage_samples = sublineage_samples
            self.sample_sublineages = sample_sublineages
            self.sublineage_allele_counts = sublineage_allele_counts

    return PopulationAllelesInfo(allele_incidence=allele_incidence_per_contig,
                                 sublineage_samples=sublineage_samples,
                                 sample_sublineages=sample_sublineages,
                                 sublineage_allele_counts=sl_allele_counts)


# call this function with tuples of (c, c_p). If you have distribution as a dictionary, call with dict.viewitems()
# this is because dictionaries don't have any order to the keys, while a list of tuple does.
def sampleFromDiscreteDistribution(category_levels, normalised=True):
    from random import random, uniform

    x = random() if normalised else uniform(0, sum(cl[1] for cl in category_levels))
    for c, l in category_levels:
        if l < x:
            x -= l
        else:
            return c
    assert "No way should we have reached here! Maybe some of category level value are -ve??"


def present_distribution(distr, sig_fig=3):
    if isinstance(distr, dict):
        distr = distr.viewitems()
    noise = 10 ** -sig_fig
    return tuple((c, round(v, sig_fig)) for c, v in sorted(distr, key=lambda kv: kv[1], reverse=True) if v >= noise)
