import sys
import os
import re
import pandas as pd
import numpy as np
import glob

if len(sys.argv) < 2:
    sys.stderr.write('Usage: %s <metadata>\n' % sys.argv[0])
    sys.exit(1)
fname_meta = sys.argv[1]

basedir = os.path.dirname(fname_meta)
fnames_cnt = glob.glob(os.path.join(basedir, 'reads', '*_1.cnt'))

metadata = pd.read_csv(fname_meta, sep='\t', index_col=0)
metadata.index.rename('transcript_id', inplace=True)
for fname in fnames_cnt:
    sample = os.path.basename(fname).split('.')[0][:-2] + '_readcount'
    tmp = pd.read_csv(fname, sep='\t', index_col=0, header=None)
    metadata = metadata.join(tmp.rename(columns={1:sample}))
    metadata[sample] = metadata[sample].fillna(0).astype('int')
metadata.to_csv(re.sub(r'.tsv$', '', fname_meta) + '.augmented.tsv', sep='\t') 
