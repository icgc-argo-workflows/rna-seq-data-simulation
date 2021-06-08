#!/bin/bash

set -e

tx_fa=/cluster/work/grlab/projects/projects2020-ICGC-ARGO/simulation/annotation/gencode.v37.pc_transcripts.fa
tx_gtf=/cluster/work/grlab/projects/projects2020-ICGC-ARGO/simulation/annotation/gencode.v37.annotation.gtf

python pre_filter_transcripts.py $tx_fa $tx_gtf 
