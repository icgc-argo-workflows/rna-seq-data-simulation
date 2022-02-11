#!/usr/bin/env python
import sys
import os
import re
import numpy as np
import numpy.random as npr
import gzip
from intervaltree import Interval, IntervalTree
import glob
import pickle

#TMB = 3*10**-6

def main():

    if len(sys.argv) < 8:
        sys.stderr.write('Usage: %s <transcripts.gtf> <transcripts.fa> <outname_tag> <outdir> <random_seed> <TMB> <phased_germline_var.vcf-pattern> <somatic_va.vcf-pattern>\n' % sys.argv[0])
        sys.exit(1)

    fname_gtf = sys.argv[1]
    fname_fa = sys.argv[2]
    outtag = sys.argv[3]
    outdir = sys.argv[4]
    random_seed = int(sys.argv[5])
    TMB = float(sys.argv[6])
    germline_vcf_pattern = sys.argv[7]
    somatic_vcf_pattern = sys.argv[8]

    transcripts, tx_intervals, total_length = get_transcripts(fname_gtf)

    sys.stderr.write('Loading germline variants\n')
    variants_germline = load_variants_1kg(germline_vcf_pattern)
    sys.stderr.write('Loading somatic variants\n')
    somatic_pickle = os.path.join(os.path.dirname(somatic_vcf_pattern), 'cosmic.pickle')
    if not os.path.exists(somatic_pickle):
        variants_somatic = load_variants_cosmic(somatic_vcf_pattern, tx_intervals, variants_germline)
        pickle.dump(variants_somatic, open(somatic_pickle, 'wb'))
    else:
        variants_somatic = pickle.load(open(somatic_pickle, 'rb'))
    num_variants_somatic = sum([len(variants_somatic[_]) for _ in variants_somatic])
    sys.stderr.write('Parsed %i somatic variants\n' % num_variants_somatic)

    somatic_fraction = int(TMB * total_length) / num_variants_somatic
    sys.stderr.write('Filtering somatic variants to reach TMB\n')
    sys.stderr.write('Probability to keep a somatic mutation: %.6f%%\n' % (somatic_fraction * 100))
    npr.seed(random_seed)
    for chrm in variants_somatic:
        keep = IntervalTree()
        for v in variants_somatic[chrm]:
            if npr.random() < somatic_fraction:
                keep.add(v)
        variants_somatic[chrm] = keep
    num_variants_somatic = sum([len(variants_somatic[_]) for _ in variants_somatic])
    sys.stderr.write('Kept %i somatic variants\n' % num_variants_somatic)

    personalize_transcriptome(fname_fa, outtag, outdir, transcripts, variants_germline, variants_somatic) 

    write_variant_positions(fname_fa, outtag, outdir, variants_germline, variants_somatic)


def mutate_seq(seq, t_pos, var, exons):

    ### check that we are mutating the right bases
    v_pos_b = t_pos + var.begin - exons[0]
    v_pos_e = v_pos_b + len(var.data[0])
    assert seq[v_pos_b:v_pos_e] == var.data[0]

    if v_pos_b == 0:
        return var.data[1] + seq[v_pos_e:], t_pos + len(var.data[1]) - len(var.data[0])
    elif v_pos_b >= len(seq) - len(var.data[0]):
        return seq[:-len(var.data[0])] + var.data[1], t_pos + len(var.data[1]) - len(var.data[0])
    else:
        return seq[:v_pos_b] + var.data[1] + seq[v_pos_e:], t_pos + len(var.data[1]) - len(var.data[0])


def rev_comp(seq):
    D = {'A':'T',
         'T':'A',
         'G':'C',
         'C':'G',
         'N':'N'}
    return ''.join([D[_] for _ in seq][::-1])


def personalize_transcript(header, seq, transcripts, variants_germ, variants_som):

    ### get transcript ID
    t_id = header[1:].split('|')[0]

    ### compute total transcript length and check that it matches up
    exons = np.array(transcripts[t_id][2])
    strand = transcripts[t_id][1]
    s_idx = np.argsort(exons[:, 0])
    exons = exons[s_idx, :]

    t_len = sum(exons[:, 1] - exons[:, 0])
    assert t_len == len(seq)

    if strand == '-':
        seq = rev_comp(seq)

    ### initialize haplotypes
    seq_ht1_g = seq
    seq_ht2_g = seq
    seq_ht1_gs = seq
    seq_ht2_gs = seq

    t_start = exons.min()
    t_pos1_g = 0
    t_pos2_g = 0
    t_pos1_gs = 0
    t_pos2_gs = 0
    var_cnt1_g = 0
    var_cnt2_g = 0
    var_cnt1_gs = 0
    var_cnt2_gs = 0

    germ_het = []
    som_het = []

    for i in range(exons.shape[0]):
        ### check whether we can apply germline variation
        if transcripts[t_id][0] in variants_germ:
            ### get matching germline variants
            for var in variants_germ[transcripts[t_id][0]][exons[i, 0]:exons[i, 1]]:
                gt = var.data[2].split('|')
                if gt[0] == '1':
                    seq_ht1_g, t_pos1_g = mutate_seq(seq_ht1_g, t_pos1_g, var, exons[i, :])
                    seq_ht1_gs, t_pos1_gs = mutate_seq(seq_ht1_gs, t_pos1_gs, var, exons[i, :])
                    var_cnt1_g += 1
                if gt[1] == '1':
                    seq_ht2_g, t_pos2_g = mutate_seq(seq_ht2_g, t_pos2_g, var, exons[i, :])
                    seq_ht2_gs, t_pos2_gs = mutate_seq(seq_ht2_gs, t_pos2_gs, var, exons[i, :])
                    var_cnt2_g += 1
                if not gt in ['1|1', '0|0']:
                    germ_het.append((transcripts[t_id][0], var.begin)) # tuple: chrm, pos
        ### check whether we can apply somatic variation
        if transcripts[t_id][0] in variants_som:
            for var in variants_som[transcripts[t_id][0]][exons[i, 0]:exons[i, 1]]:
                ### skip if we overlap with a germline position
                if transcripts[t_id][0] in variants_germ and len(variants_germ[transcripts[t_id][0]][var.begin]) > 0:
                    continue
                ### only do short deletions
                if len(var.data[0]) > 3:
                    continue
                gt = var.data[2].split('|')
                if gt[0] == '1':
                    seq_ht1_gs, t_pos1_gs = mutate_seq(seq_ht1_gs, t_pos1_gs, var, exons[i, :])
                    var_cnt1_gs += 1
                if gt[1] == '1':
                    seq_ht2_gs, t_pos2_gs = mutate_seq(seq_ht2_gs, t_pos2_gs, var, exons[i, :])
                    var_cnt2_gs += 1
                if not gt in ['1|1', '0|0']:
                    som_het.append((transcripts[t_id][0], var.begin)) # tuple: chrm, pos
        ### move in transcript
        t_pos1_g += exons[i, 1] - exons[i, 0]
        t_pos2_g += exons[i, 1] - exons[i, 0]
        t_pos1_gs += exons[i, 1] - exons[i, 0]
        t_pos2_gs += exons[i, 1] - exons[i, 0]

    if strand == '-':
        seq_ht1_g = rev_comp(seq_ht1_g)
        seq_ht2_g = rev_comp(seq_ht2_g)
        seq_ht1_gs = rev_comp(seq_ht1_gs)
        seq_ht2_gs = rev_comp(seq_ht2_gs)

    return seq_ht1_g, seq_ht2_g, seq_ht1_gs, seq_ht2_gs, var_cnt1_g, var_cnt2_g, var_cnt1_gs, var_cnt2_gs, germ_het, som_het


def write_transcript(header, seq, fh, width=60):
    
    fh.write(header + '\n')
    for i in range(0, len(seq), width):
        fh.write(seq[i:i+width] + '\n')

def write_variant_positions(fname_fa, outtag, outdir, variants_germline, variants_somatic):

    info = '##fileformat=VCFv4.3\n##INFO=<ID=VT,Number=.,Type=String,Description="indicates what type of variant the line represents">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased Genotype">'
    header = '#' + '\t'.join(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'simulation'])

    outbase = os.path.join(outdir, re.sub(r'.fa$', '', os.path.basename(fname_fa))) + '.' + outtag

    ### write germline in vcf
    with open(outbase + '.germline.vcf', 'w') as fh:
        fh.write(info + '\n')
        for chrm in variants_germline:
            fh.write('##contig=<ID=%s>\n' % chrm)
        fh.write(header + '\n')
        for chrm in variants_germline:
            for v in sorted(variants_germline[chrm]):
                fh.write('\t'.join([chrm, str(v.begin + 1), '.', v.data[0], v.data[1], '.', 'PASS', 'VT=SIM', 'GT', v.data[2]]) + '\n') 
    ### write germline in pickle
    pickle.dump(variants_germline, open(outbase + '.germline.pickle', 'wb'))

    ### write somatic in vcf
    with open(outbase + '.somatic.vcf', 'w') as fh:
        fh.write(info + '\n')
        for chrm in variants_somatic:
            fh.write('##contig=<ID=%s>\n' % chrm)
        fh.write(header + '\n')
        for chrm in variants_somatic:
            for v in sorted(variants_somatic[chrm]):
                fh.write('\t'.join([chrm, str(v.begin + 1), '.', v.data[0], v.data[1], '.', 'PASS', 'VT=SIM', 'GT', v.data[2]]) + '\n') 
    ### write somatic in pickle
    pickle.dump(variants_somatic, open(outbase + '.somatic.pickle', 'wb'))


def personalize_transcriptome(fname_fa, outtag, outdir, transcripts, variants_germline, variants_somatic):
    
    header = None
    seq = []

    outbase = os.path.join(outdir, re.sub(r'.fa$', '', os.path.basename(fname_fa))) + '.' + outtag
    fh_ht1 = open(outbase + '.normal_ht1.fa', 'w')
    fh_ht2 = open(outbase + '.normal_ht2.fa', 'w')
    fh_ht1_som = open(outbase + '.tumor_ht1.fa', 'w')
    fh_ht2_som = open(outbase + '.tumor_ht2.fa', 'w')

    var_cnt1_g = 0
    var_cnt2_g = 0
    var_cnt1_gs = 0
    var_cnt2_gs = 0
    
    germline_hets = dict()
    somatic_hets = dict()

    for line in open(fname_fa, 'r'):
        if line.startswith('>'):
            if len(seq) > 0:
                seq_ht1_g, seq_ht2_g, seq_ht1_gs, seq_ht2_gs, cnt1_g, cnt2_g, cnt1_gs, cnt2_gs, germ_het, som_het = personalize_transcript(header, ''.join(seq), transcripts, variants_germline, variants_somatic)
                var_cnt1_g += cnt1_g
                var_cnt2_g += cnt2_g
                var_cnt1_gs += cnt1_gs
                var_cnt2_gs += cnt2_gs
                for p in germ_het:
                    if not p in germline_hets:
                        germline_hets[p] = []
                    germline_hets[p].append(header[1:])
                for p in som_het:
                    if not p in somatic_hets:
                        somatic_hets[p] = []
                    somatic_hets[p].append(header[1:])
                write_transcript(header, seq_ht1_g, fh_ht1)
                write_transcript(header, seq_ht2_g, fh_ht2)
                write_transcript(header, seq_ht1_gs, fh_ht1_som)
                write_transcript(header, seq_ht2_gs, fh_ht2_som)
            seq = []
            header = line.strip()
            continue
        seq.append(line.strip())
    if len(seq) > 0:
        seq_ht1_g, seq_ht2_g, seq_ht1_gs, seq_ht2_gs, cnt1_g, cnt2_g, cnt1_gs, cnt2_gs, germ_het, som_het = personalize_transcript(header, ''.join(seq), transcripts, variants_germline, variants_somatic)
        var_cnt1_g += cnt1_g
        var_cnt2_g += cnt2_g
        var_cnt1_gs += cnt1_gs
        var_cnt2_gs += cnt2_gs
        for p in germ_het:
            if not p in germline_hets:
                germline_hets[p] = []
            germline_hets[p].append(header[1:])
        for p in som_het:
            if not p in somatic_hets:
                somatic_hets[p] = []
            somatic_hets[p].append(header[1:])
        write_transcript(header, seq_ht1_g, fh_ht1)
        write_transcript(header, seq_ht2_g, fh_ht2)
        write_transcript(header, seq_ht1_gs, fh_ht1_som)
        write_transcript(header, seq_ht2_gs, fh_ht2_som)
    fh_ht1.close()
    fh_ht2.close()
    fh_ht1_som.close()
    fh_ht2_som.close()

    pickle.dump(germline_hets, open(outbase + '.tx_hets_germline.pickle', 'wb'))
    pickle.dump(somatic_hets, open(outbase + '.tx_hets_somatic.pickle', 'wb'))

    sys.stderr.write('\ngermline variants introduced in HT1: %i\n' % var_cnt1_g)
    sys.stderr.write('germline variants introduced in HT2: %i\n' % var_cnt2_g)
    sys.stderr.write('somatic variants introduced in HT1: %i\n' % var_cnt1_gs)
    sys.stderr.write('somatic variants introduced in HT2: %i\n' % var_cnt2_gs)


def load_variants_cosmic(fname_vcf_pattern, tx_intervals, variants_germline, loh_frac=0.05):

    total_vars = 0 
    variants = dict()

    fnames_vcf = sorted(glob.glob(fname_vcf_pattern))
    for fname_vcf in fnames_vcf:
        sys.stderr.write('Parsing %s for variant information\n' % fname_vcf)
        if fname_vcf.endswith('gz'):
            fh = gzip.open(fname_vcf, 'rt', encoding='utf-8')
        else:
            fh = open(fname_vcf, 'r')

        for i, line in enumerate(fh):
            if i > 0 and i % 100000 == 0:
                sys.stderr.write('.')
                if i % 1000000 == 0:
                    sys.stderr.write('%i\n' % i)
                sys.stderr.flush()
            ### skip comments    
            if line[0] == '#':
                continue

            sl = line.strip().split('\t')
            ### get data
            chrm = 'chr' + sl[0]
            ref = sl[3]
            alt = sl[4]
            pos = int(sl[1]) - 1

            ### do not parse variants that do not overlap to the transcriptome
            if not chrm in tx_intervals or len(tx_intervals[chrm][pos]) == 0:
                continue

            ### ignore somatic variants that overlap with a germline variant
            if chrm in variants_germline and len(variants_germline[chrm][pos]) > 0:
                continue

            ### randomly select a haplotype for the variant to occur on
            if npr.random_sample() < loh_frac:
                gt = '1|1'
            elif npr.random_sample() < 0.5:
                gt = '0|1'
            else:
                gt = '1|0'

            ### as a simplification, we will only store one (the last) variant in case
            ### multiple variants occur at the same position
            if not chrm in variants:
                variants[chrm] = IntervalTree()
            variants[chrm][pos:pos+len(ref)] = [ref, alt, gt]
            total_vars += 1
        fh.close()

    sys.stderr.write('\nExtracted %i somatic variants (overlapping to the provided transcriptome)\n\n' % total_vars)

    return variants

def load_variants_1kg(fname_vcf_pattern):

    chr_dict = {'ref|NC_000001.11|': 'chr1',
                'ref|NC_000002.12|': 'chr2',
                'ref|NC_000003.12|': 'chr3',
                'ref|NC_000004.12|': 'chr4',
                'ref|NC_000005.10|': 'chr5',
                'ref|NC_000006.12|': 'chr6',
                'ref|NC_000007.14|': 'chr7',
                'ref|NC_000008.11|': 'chr8',
                'ref|NC_000009.12|': 'chr9',
                'ref|NC_000010.11|': 'chr10',
                'ref|NC_000011.10|': 'chr11',
                'ref|NC_000012.12|': 'chr12',
                'ref|NC_000013.11|': 'chr13',
                'ref|NC_000014.9|': 'chr14',
                'ref|NC_000015.10|': 'chr15',
                'ref|NC_000016.10|': 'chr16',
                'ref|NC_000017.11|': 'chr17',
                'ref|NC_000018.10|': 'chr18',
                'ref|NC_000019.10|': 'chr19',
                'ref|NC_000020.11|': 'chr20',
                'ref|NC_000021.9|': 'chr21',
                'ref|NC_000022.11|': 'chr22',
                'ref|NC_000023.11|': 'chrX'}

    total_vars = 0
    variants = dict()

    fnames_vcf = sorted(glob.glob(fname_vcf_pattern))
    for fname_vcf in fnames_vcf:
        sys.stderr.write('Parsing %s for variant information\n' % fname_vcf)
        if fname_vcf.endswith('gz'):
            fh = gzip.open(fname_vcf, 'rt', encoding='utf-8')
        else:
            fh = open(fname_vcf, 'r')

        for i, line in enumerate(fh):
            if i > 0 and i % 10000 == 0:
                sys.stderr.write('.')
                if i % 100000 == 0:
                    sys.stderr.write('%i\n' % i)
                sys.stderr.flush()
            ### skip comments    
            if line[0] == '#':
                continue

            sl = line.strip().split('\t')
            ### get data
            chrm = chr_dict[sl[0]]
            ref = sl[3]
            alt = sl[4]
            pos = int(sl[1]) - 1
            ### NOTE: we have hard-coded the column here to match the 1KG vcfs
            gt = sl[9]

            if not chrm in variants:
                variants[chrm] = IntervalTree()
            variants[chrm][pos:pos+len(ref)] = [ref, alt, gt]
            total_vars += 1
        fh.close()

    sys.stderr.write('\nExtracted %i variants\n\n' % total_vars)

    return variants

def get_transcripts(fname_gtf):

    ### parse all transcripts and store by ID
    transcripts = dict()
    tx_intervals = dict()
    all_pos = dict()
    sys.stderr.write('Parsing %s for transcript coordinates\n' % fname_gtf)
    for i, line in enumerate(open(fname_gtf, 'r')):

        if i > 0 and i % 100000 == 0:
            sys.stderr.write('.')
            if i % 1000000 == 0:
                sys.stderr.write('%i\n' % i)
            sys.stderr.flush()

        ### skip comments
        if line[0] == '#':
            continue

        sl = line.strip().split('\t')

        ### only look at exons
        if sl[2] != 'exon':
            continue

        chrm = sl[0]

        ### get tags
        tags = dict([(_.strip(' ').split(' ')[0], _.strip(' ').split(' ')[1].strip('"')) for _ in sl[8].strip(';').split(';')])
        
        if tags['transcript_id'] in transcripts:
            transcripts[tags['transcript_id']][-1].append([int(sl[3]) - 1, int(sl[4])])
            assert sl[6] == transcripts[tags['transcript_id']][1]
        else:
            transcripts[tags['transcript_id']] = [sl[0], sl[6], [[int(sl[3]) - 1, int(sl[4])]]]

        if not chrm in all_pos:
            all_pos[chrm] = []
        all_pos[chrm].append([int(sl[3]) - 1, int(sl[4])])

        if not chrm in tx_intervals:
            tx_intervals[chrm] = IntervalTree()
        tx_intervals[chrm][int(sl[3]) - 1:int(sl[4])] = True
    
    total_length = 0
    for chrm in all_pos:
        total_length += len(set().union(*[set(range(_[0], _[1])) for _ in all_pos[chrm]]))

    sys.stderr.write('\nExtracted coordinates for %i transcripts\nTotal genome positions in transcriptome: %i\n\n' % (len(transcripts), total_length))

    return transcripts, tx_intervals, total_length


if __name__ == "__main__":
    main()
