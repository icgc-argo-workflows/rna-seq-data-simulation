#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/projects2020-ICGC-ARGO/simulation/annotation
tx_fa= ${basedir}/gencode.v37.pc_transcripts.fa
tx_gtf=${basedir}/gencode.v37.annotation.gtf
tx_gtf_filt=${tx_gtf%.gtf}.uniq.min400.gtf

### perform transcript filtering
if [ ! -f ${tx_gtf_filt} ]
then
    python pre_filter_transcripts.py $tx_fa $tx_gtf 
fi

### collect gene lengths
for tx in $tx_gtf $tx_gtf_filt 
do
    if [ ! -f ${tx%.gtf}.gene_lens.tsv ]
    then
        python get_gene_length.py ${tx} > ${tx%.gtf}.gene_lens.tsv
    fi
done
