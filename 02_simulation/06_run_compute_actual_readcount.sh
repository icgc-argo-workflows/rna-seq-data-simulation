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

for i in ${!GENOMES[@]}
do
    currbase=${BASEDIR}/${SIM_ID}/${GENOMES[$i]}_${ICGC_SAMPLES[$i]}
    for ct in normal tumor
    do
        for j in $(seq 1 ${NUM_SAMPLES})
        do
            fastq=${currbase}/reads/${ct}_sample_0${j}_1.fastq.gz
            if [ ! -f ${fastq%.fastq.gz}.cnt ]
            then
                #sbatch -c ${threads} --time=20:00:00 --mem=${mem} --output=/dev/null --job-name=cnt --wrap "python $(pwd)/src/08_compute_actual_readcount.py ${fastq}"
                python $(pwd)/src/08_compute_actual_readcount.py ${fastq}
            else
                echo ${fastq%.fastq.gz}.cnt exists
            fi
        done
    done
done
