#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/projects2020-ICGC-ARGO/simulation/results/polyester/full_transcriptome_v1
mem=30000
threads=1
pmem=$((${mem}/${threads}))

genomes=("HG00100" "HG00101" "HG00102" "HG00103" "HG00104")
icgc_samples=("8b608b7a-d4d5-11e4-a827-bd90bc314f72" "102fdc84-e1f6-11e4-b21b-b8030ae66825" "e790f616-dd3e-11e4-b051-c140c254ba06" "31615db2-e1f9-11e4-9c66-c7ab578961f5" "1e789636-fb5b-11e4-ab35-d5f0b89b4d05")
rand_seeds=(23 24 25 26 27)

for i in ${!genomes[@]}
do
    outdir=${basedir}/${genomes[$i]}_${icgc_samples[$i]}
    mkdir -p ${outdir}
    log1=${outdir}/create_simulation_setup.log
    log2=${outdir}/create_simulation_setup.lsf.log
    echo "cd /cluster/home/akahles/git/projects/2021/icgc_argo_rna_data_simulation/simulation; ./create_simulation_setup.sh ${genomes[$i]} ${icgc_samples[$i]} ${rand_seeds[$i]} 2>&1 > ${log1}" | bsub -G ms_raets -M ${mem} -n ${threads} -We 20:00 -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -J sim_${i} -oo $log2
done
