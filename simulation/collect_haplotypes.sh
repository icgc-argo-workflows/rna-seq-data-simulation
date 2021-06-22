#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/projects2020-ICGC-ARGO/simulation/results/polyester/full_transcriptome_v1

genomes=("HG00100" "HG00101" "HG00102" "HG00103" "HG00104")
icgc_samples=("8b608b7a-d4d5-11e4-a827-bd90bc314f72" "102fdc84-e1f6-11e4-b21b-b8030ae66825" "e790f616-dd3e-11e4-b051-c140c254ba06" "31615db2-e1f9-11e4-9c66-c7ab578961f5" "1e789636-fb5b-11e4-ab35-d5f0b89b4d05")

outdir=${basedir}/reads
mkdir -p $outdir

for i in ${!genomes[@]}
do
    echo collecting ${genomes[$i]}_${icgc_samples[$i]}
    currbase=${basedir}/${genomes[$i]}_${icgc_samples[$i]}
    for ct in normal tumor
    do
        for j in $(seq 1 3)
        do
            for mate in 1 2
            do
                echo collecting ${genomes[$i]}_${icgc_samples[$i]} $ct sample $j mate $mate
                (cat ${currbase}/reads/${ct}_ht1_sample_0${j}_${mate}.fastq.gz ${currbase}/reads/${ct}_ht2_sample_0${j}_${mate}.fastq.gz > ${currbase}/reads/${ct}_sample_0${j}_${mate}.fastq.gz) && rm ${currbase}/reads/${ct}_ht1_sample_0${j}_${mate}.fast[aq].gz ${currbase}/reads/${ct}_ht2_sample_0${j}_${mate}.fast[aq].gz
                echo linking sample
                ln -s ${currbase}/reads/${ct}_sample_0${j}_${mate}.fastq.gz ${outdir}/${genomes[$i]}_${icgc_samples[$i]}_${ct}_sample_0${j}_${mate}.fastq.gz
            done
        done
    done
done

