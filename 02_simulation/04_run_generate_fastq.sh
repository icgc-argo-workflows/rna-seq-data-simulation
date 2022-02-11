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
mem=6000
threads=1
pmem=$((${mem} / ${threads}))

tx_base=$(basename $ANNO_GTF)
tx_base=${tx_base%.gtf}
for i in ${!GENOMES[@]}
do
    currbase=${BASEDIR}/${SIM_ID}/${GENOMES[$i]}_${ICGC_SAMPLES[$i]}
    for ht in ht1 ht2
    do
        for ct in normal tumor
        do
            for j in $(seq 1 ${NUM_SAMPLES})
            do
                for mate in 1 2
                do
                    tx_fasta=${currbase}/${tx_base}.${GENOMES[$i]}.${ct}_${ht}.expressed.fa 
                    fasta=${currbase}/reads/${ct}_${ht}_sample_0${j}_${mate}.fasta.gz
                    log=${fasta%.fasta.gz}.fa2fq.lsf.log
                    if [ ! -f ${fasta%.fasta.gz}.fastq.gz ]
                    then
                        echo "python $(pwd)/src/07_fasta2fastq.py ${fasta} ${QUALITY_STATS} ${tx_fasta}" | bsub -G ms_raets -M ${mem} -n ${threads} -We 20:00 -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -J ${GENOMES[$i]}_${ht}_${ct} -oo $log
                    else
                        echo ${fasta%.fasta.gz}.fastq.gz exists
                    fi
                done
            done
        done
    done
done
