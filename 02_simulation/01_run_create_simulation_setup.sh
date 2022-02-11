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

### force absolute path
case "$config" in
    /*) ;;
    *) config=$(pwd)/${config} ;;
esac

### cluster params
mem=30000
threads=1
pmem=$((${mem}/${threads}))

outbase=${BASEDIR}/${SIM_ID}
for i in ${!GENOMES[@]}
do
    outdir=${outbase}/${GENOMES[$i]}_${ICGC_SAMPLES[$i]}
    mkdir -p ${outdir}
    log1=${outdir}/create_simulation_setup.log
    log2=${outdir}/create_simulation_setup.lsf.log
    echo "cd $(pwd)/src; ./01_create_simulation_setup.sh ${config} ${GENOMES[$i]} ${ICGC_SAMPLES[$i]} ${RAND_SEEDS[$i]} $TMB 2>&1 > ${log1}" | bsub -G ms_raets -M ${mem} -n ${threads} -We 20:00 -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -J sim_${i} -oo $log2
done
