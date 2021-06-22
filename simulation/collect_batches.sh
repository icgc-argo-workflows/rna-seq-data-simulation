#!/bin/bash

set -e

threads=1
mem=8000
pmem=$((${mem} / ${threads}))

basedir=/cluster/work/grlab/projects/projects2020-ICGC-ARGO/simulation/results/polyester/full_transcriptome_v1

genomes=("HG00100" "HG00101" "HG00102" "HG00103" "HG00104")
icgc_samples=("8b608b7a-d4d5-11e4-a827-bd90bc314f72" "102fdc84-e1f6-11e4-b21b-b8030ae66825" "e790f616-dd3e-11e4-b051-c140c254ba06" "31615db2-e1f9-11e4-9c66-c7ab578961f5" "1e789636-fb5b-11e4-ab35-d5f0b89b4d05")

for i in ${!genomes[@]}
do
    echo collecting ${genomes[$i]}_${icgc_samples[$i]}
    outdir=${basedir}/${genomes[$i]}_${icgc_samples[$i]}/reads
    mkdir -p ${outdir}
    for ht in ht1 ht2
    do
        for ct in normal tumor
        do
            batchbase=${basedir}/${genomes[$i]}_${icgc_samples[$i]}/batches/${ct}_${ht}
            for j in $(seq 1 3)
            do
                tx_fa=${basedir}/${genomes[$i]}_${icgc_samples[$i]}/*.${ct}_${ht}.expressed.fa
                fa_pattern=${batchbase}/'batch_*'/sample_0${j}_1.fasta.gz
                outbase=${outdir}/${ct}_${ht}_sample_0${j}
                if [ ! -f ${outbase}_1.fasta.gz -o ! -f ${outbase}_2.fasta.gz ]
                then
                    echo "python $(pwd)/collect_and_annotate_fastq.py $tx_fa $outbase \"$fa_pattern\"" | bsub -G ms_raets -M ${mem} -n ${threads} -We 20:00 -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -J ${genomes[$i]}_${ht}_${ct} -oo ${outbase}.collect.lsf.log
                fi
            done
        done
    done
    wait
done
