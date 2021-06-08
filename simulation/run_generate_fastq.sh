#!/bin/bash

set -e

mem=6000
threads=1
pmem=$((${mem} / ${threads}))

basedir=/cluster/work/grlab/projects/projects2020-ICGC-ARGO/simulation/results/polyester/full_transcriptome_v1
stats=/cluster/work/grlab/projects/projects2020-ICGC-ARGO/simulation/annotation/05caea96-dd31-11e4-8808-f884c254ba06.aligned.stats.tsv.pickle

genomes=("HG00100" "HG00101" "HG00102" "HG00103" "HG00104")
icgc_samples=("8b608b7a-d4d5-11e4-a827-bd90bc314f72" "102fdc84-e1f6-11e4-b21b-b8030ae66825" "e790f616-dd3e-11e4-b051-c140c254ba06" "31615db2-e1f9-11e4-9c66-c7ab578961f5" "1e789636-fb5b-11e4-ab35-d5f0b89b4d05")

for i in ${!genomes[@]}
do
    currbase=${basedir}/${genomes[$i]}_${icgc_samples[$i]}
    for ht in ht1 ht2
    do
        for ct in normal tumor
        do
            if [ "${genomes[$i]}_${ct}" != "HG00100_tumor" ]
            then
                continue
            fi
            for j in $(seq 1 3)
            do
                for mate in 1 2
                do
                    tx_fasta=${currbase}/gencode.v37.pc_transcripts.uniq.min400nt.${genomes[$i]}.${ct}_${ht}.expressed.fa 
                    fasta=${currbase}/reads/${ct}_${ht}_sample_0${j}_${mate}.fasta.gz
                    log=${fasta%.fasta.gz}.fa2fq.lsf.log
                    if [ ! -f ${fasta%.fasta.gz}.fastq.gz ]
                    then
                        echo "python $(pwd)/fasta2fastq.py ${fasta} ${stats} ${tx_fasta}" | bsub -G ms_raets -M ${mem} -n ${threads} -We 20:00 -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -J ${genomes[$i]}_${ht}_${ct} -oo $log
                    else
                        echo ${fasta%.fasta.gz}.fastq.gz exists
                    fi
                done
            done
        done
    done
done
