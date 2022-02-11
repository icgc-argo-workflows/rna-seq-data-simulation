#!/usr/bin/env python

# Pre-filter the annotation for a minimum trranscript length of 
# 400 nt and unique transcript sequences (otherwise polyester 
# has difficulties simulating)

import sys
import re

if len(sys.argv) < 3:
    sys.stderr.write('Usage: %s <transcripts.fa> <annotation.gtf>\n' % sys.argv[0])
    sys.exit(1)
fname_tx_fa = sys.argv[1]
fname_anno_gtf = sys.argv[2]

MINLEN = 400

### pass through transcript sequences and filter for length
### and sequence uniqueness
all_transcripts = set([])
header = ''
seq = []
skipped = 0
kept_tx = set()
kept_ge = set()
out = open(re.sub(r'.fa$', '', fname_tx_fa) + '.uniq.min%int.fa' % MINLEN, 'w')
for line in open(fname_tx_fa, 'r'):
    if line.startswith('>'): 
        if len(seq) > 0:
            tx = ''.join([_.strip() for _ in seq])
            if not tx in all_transcripts and len(tx) >= MINLEN:
                out.write(header + ''.join(seq))
                all_transcripts.add(tx)
                kept_tx.add(header.split('|')[0][1:])
                kept_ge.add(header.split('|')[1])
            else:
                skipped += 1
        ### cut header
        header = '|'.join(line.strip().split('|')[:4]) + '\n'
        seq = []
    else:
        seq.append(line)
if len(seq) > 0:
    tx = ''.join([_.strip() for _ in seq])
    if not tx in all_transcripts and len(tx) >= MINLEN:
        out.write(header + ''.join(seq))
        all_transcripts.add(tx)
        kept_tx.add(header.split('|')[0][1:])
        kept_ge.add(header.split('|')[1])
    else:
        skipped += 1
out.close()
sys.stderr.write('Skipped %i transcripts due to sequence duplicity or shortness (< 400nt) \n' % skipped)

### pass through gtf and only retain kept transcripts
out = open(re.sub(r'.fa$', '', fname_tx_fa) + '.uniq.min%int.gtf' % MINLEN, 'w')
for line in open(fname_anno_gtf, 'r'):
    if line.startswith('#'):
        out.write(line)
        continue
    sl = line.strip().split('\t')

    ### get tags
    tags = dict([(_.strip(' ').split(' ')[0], _.strip(' ').split(' ')[1].strip('"')) for _ in sl[8].strip(';').split(';')])

    if sl[2] == 'gene': 
        if not tags['gene_id'] in kept_ge:
            continue
    else:
        if not tags['transcript_id'] in kept_tx:
            continue
    out.write(line)
out.close()

