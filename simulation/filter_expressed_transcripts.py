import sys
import pandas as pd
import re
import numpy as np

if len(sys.argv) < 3:
    sys.stderr.write('Usage: %s <transcript.fa> <metadata.tsv>\n' % sys.argv[0])
    sys.exit(1)

fname_tx = sys.argv[1]
fname_meta = sys.argv[2]

metadata = pd.read_csv(fname_meta, sep='\t', index_col=0)
kidx = np.where(metadata['reads'] > 0)[0]
kids = set(metadata.index[kidx])

output = False
with open(re.sub(r'.fa$', '', fname_tx) + '.expressed.fa', 'w') as out:
    for line in open(fname_tx, 'r'):
        if line.startswith('>'):
            txid = line[1:].strip().split('|')[0]
            output = txid in kids
        if output:
            out.write(line)
