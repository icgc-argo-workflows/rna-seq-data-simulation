#!/bin/bash

set -e

image=/cluster/work/grlab/share/singularity_images/other/Singularity_polyester-v1.26.0.sif
mem=6000
threads=1
pmem=$(($mem / $threads))

### polyester errors out when N=1, we instead have duplicated the fold-change columns and are simulating with N=3 and R=1
N=3 # num samples
R=1 # num replicates
bias=rnaf

genomes=("HG00100" "HG00101" "HG00102" "HG00103" "HG00104")
icgc_samples=("8b608b7a-d4d5-11e4-a827-bd90bc314f72" "102fdc84-e1f6-11e4-b21b-b8030ae66825" "e790f616-dd3e-11e4-b051-c140c254ba06" "31615db2-e1f9-11e4-9c66-c7ab578961f5" "1e789636-fb5b-11e4-ab35-d5f0b89b4d05")

for i in ${!genomes[@]}
do
    for ht in ht1 ht2
    do
        for ct in normal tumor
        do
            batchbase=/cluster/work/grlab/projects/projects2020-ICGC-ARGO/simulation/results/polyester/full_transcriptome_v1/${genomes[${i}]}_${icgc_samples[${i}]}/batches/${ct}_${ht}
            for batch in $(find $batchbase -maxdepth 1 -mindepth 1 -type d -name batch\* | sort)
            do
                log1=${batch}/polyester.log
                log2=${log1%.log}.lsf.log
                donefile=${log1%.log}.done
                if [ ! -f ${donefile} ]
                then
                    echo "/usr/bin/time -v singularity exec -B /cluster/work/grlab -B /cluster/home/akahles $image $(pwd)/../polyester/bin/run_polyester.R --fasta_input ${batch}/*.fa --output_dir ${batch}/ --num_samples $N --num_replicates $R --lib_sizes 1 --fold_changes ${batch}/*factors_*.csv --read_counts ${batch}/*.read_count.csv --bias $bias && touch $donefile 2>&1 > $log1" | bsub -G ms_raets -M ${mem} -n ${threads} -We 20:00 -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -J $(basename $batch) -oo $log2 
                else
                    echo "${batch} already complete"
                fi
            done
        done
    done
done
