#!/usr/bin/env python

# Read an annotation file in GTF format and
# output gene length to stdout. Gene length
# is computed as the number of exon positions
# per gene

import sys
import gzip

if len(sys.argv) < 2:
    sys.stderr.write('Usage: %s <annotation.gtf>\n' % sys.argv[0])
    sys.exit(1)
infname = sys.argv[1]

genes = dict()
sys.stderr.write('collecting gene information\n')
if infname.endswith('gz'):
    fh = gzip.open(infname, 'rt', encoding='utf-8')
else:
    fh = open(infname, 'r')

for i, line in enumerate(fh):
    if i > 0 and i % 100000 == 0:
        sys.stderr.write('.')
        if i % 1000000 == 0:
            sys.stderr.write('%i\n' % i)
        sys.stderr.flush()

    if line[0] == '#':
        continue
    sl = line.strip().split('\t')

    ### only look at exons
    if sl[2] != 'exon':
        continue

    ### get tags
    tags = dict([(_.strip(' ').split(' ')[0], _.strip(' ').split(' ')[1].strip('"')) for _ in sl[8].strip(';').split(';')])

    if not tags['gene_id'] in genes:
        genes[tags['gene_id']] = []
    genes[tags['gene_id']].append([int(sl[3]) -1, int(sl[4])])

sys.stderr.write('\ncompute gene lengths\n')
for g in genes:
    genes[g] = len(set().union(*[set(range(_[0], _[1])) for _ in genes[g]]))
    
sys.stderr.write('\nwrite output\n')
sys.stdout.write('\t'.join(['gene_id', 'length']) + '\n')
for g in genes:
    sys.stdout.write('\t'.join([g, str(genes[g])]) + '\n')
