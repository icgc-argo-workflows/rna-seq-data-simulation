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

outdir=${BASEDIR}/${SIM_ID}/reads
mkdir -p $outdir

for i in ${!GENOMES[@]}
do
    echo collecting ${GENOMES[$i]}_${ICGC_SAMPLES[$i]}
    currbase=${BASEDIR}/${SIM_ID}/${GENOMES[$i]}_${ICGC_SAMPLES[$i]}
    for ct in normal tumor
    do
        for j in $(seq 1 ${NUM_SAMPLES})
        do
            for mate in 1 2
            do
                if [ ! -f ${currbase}/reads/${ct}_sample_0${j}_${mate}.fastq.gz ]
                then
                    echo collecting ${GENOMES[$i]}_${ICGC_SAMPLES[$i]} $ct sample $j mate $mate
                    (cat ${currbase}/reads/${ct}_ht1_sample_0${j}_${mate}.fastq.gz ${currbase}/reads/${ct}_ht2_sample_0${j}_${mate}.fastq.gz > ${currbase}/reads/${ct}_sample_0${j}_${mate}.fastq.gz) && rm ${currbase}/reads/${ct}_ht1_sample_0${j}_${mate}.fast[aq].gz ${currbase}/reads/${ct}_ht2_sample_0${j}_${mate}.fast[aq].gz
                fi
                if [ ! -f ${outdir}/${GENOMES[$i]}_${ICGC_SAMPLES[$i]}_${ct}_sample_0${j}_${mate}.fastq.gz ]
                then
                    echo linking sample
                    ln -s ${currbase}/reads/${ct}_sample_0${j}_${mate}.fastq.gz ${outdir}/${GENOMES[$i]}_${ICGC_SAMPLES[$i]}_${ct}_sample_0${j}_${mate}.fastq.gz
                fi
            done
        done
    done
done

