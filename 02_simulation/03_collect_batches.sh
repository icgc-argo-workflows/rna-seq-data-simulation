#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "Usage: $0 <config_file>"
    exit 1
fi
config=$1

### load config - providing all capitalized variables
source $config

### cluster params
threads=1
mem=8000
pmem=$((${mem} / ${threads}))

for i in ${!GENOMES[@]}
do
    echo collecting ${GENOMES[$i]}_${ICGC_SAMPLES[$i]}
    outdir=${BASEDIR}/${SIM_ID}/${GENOMES[$i]}_${ICGC_SAMPLES[$i]}/reads
    mkdir -p ${outdir}
    for ht in ht1 ht2
    do
        for ct in normal tumor
        do
            batchbase=${BASEDIR}/${SIM_ID}/${GENOMES[$i]}_${ICGC_SAMPLES[$i]}/batches/${ct}_${ht}
            for j in $(seq 1 ${NUM_SAMPLES})
            do
                tx_fa=${BASEDIR}/${SIM_ID}/${GENOMES[$i]}_${ICGC_SAMPLES[$i]}/*.${ct}_${ht}.expressed.fa
                fa_pattern=${batchbase}/'batch_*'/sample_0${j}_1.fasta.gz
                outbase=${outdir}/${ct}_${ht}_sample_0${j}
                if [ ! -f ${outbase}_1.fasta.gz -o ! -f ${outbase}_2.fasta.gz ]
                then
                    echo "python $(pwd)/src/06_collect_and_annotate_fasta.py $tx_fa $outbase $ht \"$fa_pattern\"" | bsub -G ms_raets -M ${mem} -n ${threads} -We 20:00 -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -J ${GENOMES[$i]}_${ht}_${ct} -oo ${outbase}.collect.lsf.log
                fi
            done
        done
    done
    wait
done
