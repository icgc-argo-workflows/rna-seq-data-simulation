#!/usr/bin/env python

import sys
import os

### arg parse
if len(sys.argv) < 6:
    sys.stderr.write('Usage: %s <transcripts.fa> <fold_changes.csv> <read_counts.csv> <batchsize> <outdir>\n' % sys.argv[0])
    sys.exit(1)
fname_tx = sys.argv[1]
fname_fc = sys.argv[2]
fname_rc = sys.argv[3]
batchsize = int(sys.argv[4])
outdir = sys.argv[5]

### check that outdir exists and create if not
if not os.path.exists(outdir):
    os.makedirs(outdir)

thresh = batchsize
batch = 0
with open(fname_tx, 'r') as fhandle_tx, open(fname_fc, 'r') as fhandle_fc, open(fname_rc, 'r') as fhandle_rc:
    curr_outdir = os.path.join(outdir, f'batch_{batch:05}')
    if not os.path.exists(curr_outdir):
        os.makedirs(curr_outdir)
    tx_out = open(os.path.join(curr_outdir, os.path.basename(fname_tx)), 'w')
    fc_out = open(os.path.join(curr_outdir, os.path.basename(fname_fc)), 'w')
    rc_out = open(os.path.join(curr_outdir, os.path.basename(fname_rc)), 'w')
    tx = fhandle_tx.readline()
    for cnt, fc in enumerate(fhandle_fc):
        rc = fhandle_rc.readline()
        if cnt == thresh:
            thresh += batchsize
            batch += 1
            tx_out.close()
            fc_out.close()
            rc_out.close()
            curr_outdir = os.path.join(outdir, f'batch_{batch:05}')
            if not os.path.exists(curr_outdir):
                os.makedirs(curr_outdir)
            tx_out = open(os.path.join(curr_outdir, os.path.basename(fname_tx)), 'w')
            fc_out = open(os.path.join(curr_outdir, os.path.basename(fname_fc)), 'w')
            rc_out = open(os.path.join(curr_outdir, os.path.basename(fname_rc)), 'w')
        fc_out.write(fc)
        rc_out.write(rc)
        tx_out.write(tx)
        tx = fhandle_tx.readline()
        while not tx.startswith('>') and len(tx) > 0:
            tx_out.write(tx)
            tx = fhandle_tx.readline()
tx_out.close()
fc_out.close()
rc_out.close()

