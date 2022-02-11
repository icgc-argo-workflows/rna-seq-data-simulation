import glob
import gzip
import numpy as np
import os
import pysam
import re
import sys

REV_MAP = {'A':'T',
           'T':'A',
           'G':'C',
           'C':'G',
           'N':'N',
           'a':'t',
           't':'a',
           'g':'c',
           'c':'g',
           'n':'n'}

### This script takes the headers from the Polyester fasta files and figures out 
### which mate pair each read is and whether it is fwd or rev given its sequence
def main():
    if len(sys.argv) < 5:
        sys.stderr.write('Usage: %s <transcripts.fa> <outbase> <haplotype> <fasta_R1_pattern>\n' % sys.argv[0])
        sys.exit(1)

    tx_fasta = sys.argv[1]
    outbase = sys.argv[2]
    ht = sys.argv[3]
    fa_pattern = sys.argv[4]
    
    ### get transcript sequences
    txs = _load_transcripts(tx_fasta)

    fa1_out = gzip.open(outbase + '_1.fasta.gz', 'wt', encoding='utf-8')
    fa2_out = gzip.open(outbase + '_2.fasta.gz', 'wt', encoding='utf-8')

    ### iterate over fasta pairs
    read_cnt = 1
    fastas = glob.glob(fa_pattern)
    for fasta in sorted(fastas):
        sys.stderr.write('parsing %s\n' % fasta)
        fasta1_in = fasta
        fasta2_in = re.sub(r'_1.fasta.gz$', '', fasta) + '_2.fasta.gz'
        ### we are assuming a certain naming convention here
        ### namely that read pairs are named sample_01_1.fasta.gz and sample_01_2.fasta.gz
        with gzip.open(fasta1_in, 'rt', encoding='utf-8') as fa1, gzip.open(fasta2_in, 'rt', encoding='utf-8') as fa2:
            while True:
                ### read pair
                h1 = fa1.readline().strip()
                h2 = fa2.readline().strip()
                if not h1 or not h2:
                    break
                assert h1 == h2
                s1 = fa1.readline().strip()
                s2 = fa2.readline().strip()

                ### get corresponding transcript
                header, mate1, mate2 = h1[1:].split(';')
                _, txID = header.split('/') 
                tx_seq = txs[txID]

                ### figure out mate pairing
                r1s, r1e = mate1.split(':')[1].split('-')
                r2s, r2e = mate2.split(':')[1].split('-')
                ### we do the weird int(float()) here to fix polyester coordinate outputs such es 1e+05
                r1 = tx_seq[int(float(r1s)) - 1:int(float(r1e))]
                r2 = tx_seq[int(float(r2s)) - 1:int(float(r2e))]
           
                rc1 = _rev_comp(r1)
                rc2 = _rev_comp(r2)

                ### fa1 has mate1 fwd, fa2 has mate2 rev
                if s1 == r1 and s2 == rc2: 
                    fa1_status = 'mate1:fwd'
                    fa2_status = 'mate2:rev'
                ### fa1 has mate2 rev, fa2 has mate1 fwd
                elif s1 == rc2 and s2 == r1:
                    fa1_status = 'mate2:rev'
                    fa2_status = 'mate1:fwd'
                ### fa1 has mate2 fwd, fa2 has mate1 rev
                elif s1 == r2 and s2 == rc1: 
                    fa1_status = 'mate2:fwd'
                    fa2_status = 'mate1:rev'
                ### fa1 has mate1 rev, fa2 has mate2 fwd
                elif s1 == rc1 and s2 == r2:
                    fa1_status = 'mate1:rev'
                    fa2_status = 'mate2:fwd'
                ### this should not happen
                else:
                    assert False, "The read pairs are incompatible"

                fa1_out.write(f'>{ht}-{read_cnt}|' + ';'.join([txID, mate1, mate2, fa1_status]) + '\n')
                fa1_out.write(s1 + '\n')
                fa2_out.write(f'>{ht}-{read_cnt}|' + ';'.join([txID, mate1, mate2, fa2_status]) + '\n')
                fa2_out.write(s2 + '\n')
                read_cnt += 1

    fa1_out.close()
    fa2_out.close()

def _load_transcripts(tx_fa):

    txs = dict()
    seq = []
    for line in open(tx_fa, 'r'):
        if line.startswith('>'):
            if len(seq) > 0:
                txs[header] = ''.join(seq)
            seq = []
            header = line.strip()[1:]
            continue
        seq.append(line.strip())
    if len(seq) > 0:
        txs[header] = ''.join(seq)
    return txs


def _rev_comp(s):
    """computre reverse complement of s"""

    return ''.join([REV_MAP[_] for _ in s][::-1])

if __name__ == "__main__":
    main()
