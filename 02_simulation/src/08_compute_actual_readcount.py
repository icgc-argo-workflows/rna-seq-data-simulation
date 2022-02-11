import sys
import gzip
import re

if len(sys.argv) < 2:
    sys.stderr.write('Usage: %s <fastq.gz>\n' % sys.argv[0])
    sys.exit(1)
fname = sys.argv[1]

counts = dict()
for i, line in enumerate(gzip.open(fname, 'rt', encoding='utf-8')):
    if i > 0 and i % 1000000 == 0:
        sys.stderr.write('.')
        if i % 10000000 == 0:
            sys.stderr.write('%i\n' % i)
        sys.stderr.flush()
    if not line[0] == '@':
        continue
    tx = line.split('|')[1]
    try:
        counts[tx] += 1
    except KeyError:
        counts[tx] = 1

with open(re.sub(r'.fastq.gz$', '', fname) + '.cnt', 'w') as out:
    for tx in sorted(counts):
        out.write('\t'.join([tx, str(counts[tx])]) + '\n')
