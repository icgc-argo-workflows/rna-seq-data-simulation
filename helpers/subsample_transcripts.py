#!/usr/bin/env python

import sys
import re
import numpy as np
import numpy.random as npr
npr.seed(23)

if len(sys.argv) < 4:
    sys.stderr.write('Usage: %s <transcripts.fa> <annotation.gtf> <samplesize>\n' % sys.argv[0])
    sys.exit(1)
fname_fa = sys.argv[1]
fname_gtf = sys.argv[2]
samplesize = int(sys.argv[3])

### count genes
genes = []
for line in open(fname_fa, 'r'):
    if line.startswith('>'):
        genes.append(line.split('|')[1])

genes_u = np.unique(genes)
cnt = genes_u.shape[0]
if cnt < samplesize:
    sys.stderr.write('WARNING: sample size (%i) is larger than number of available genes (%i) - using all genes.\n' % (samplesize, cnt)) 
    samplesize = cnt

### generate random subsample
idx = np.sort(npr.choice(np.arange(cnt), samplesize, replace=False))
targets = set(np.array(genes_u)[idx])
assert len(targets) == samplesize

### subsample fasta
outfname_fa = re.sub(r'.fa$', '', fname_fa) + '.%i_genes.fa' % samplesize
with open(outfname_fa, 'w') as out:
    for i, line in enumerate(open(fname_fa, 'r')):
        if line.startswith('>'):
            curr_gene = line.split('|')[1]
        if curr_gene in targets:
            out.write(line)

### subsample gtf
outfname_gtf = re.sub(r'.gtf$', '', fname_gtf) + '.%i_genes.gtf' % samplesize
with open(outfname_gtf, 'w') as out:
    out.write('##subsample from %s of size %i genes\n' % (fname_gtf, samplesize))
    for line in open(fname_gtf, 'r'):
        if line.startswith('#'):
            out.write(line)
            continue
        sl = line.strip().split('\t')
        tags = dict([(_.strip().split(' ')[0], _.split(' ')[1].strip('"')) for _ in sl[8].strip(';').split(';')])
        if tags['gene_id'] in targets:
            out.write(line)
