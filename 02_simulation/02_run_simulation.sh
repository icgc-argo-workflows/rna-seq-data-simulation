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
pmem=$(($mem / $threads))

for i in ${!GENOMES[@]}
do
    for ht in ht1 ht2
    do
        for ct in normal tumor
        do
            batchbase=${BASEDIR}/${SIM_ID}/${GENOMES[${i}]}_${ICGC_SAMPLES[${i}]}/batches/${ct}_${ht}
            for batch in $(find $batchbase -maxdepth 1 -mindepth 1 -type d -name batch\* | sort)
            do
                log1=${batch}/polyester.log
                log2=${log1%.log}.lsf.log
                donefile=${log1%.log}.done
                if [ ! -f ${donefile} ]
                then
                    echo "/usr/bin/time -v singularity exec $SINGULARITY_PARAMS $POLYESTER_IMAGE $(pwd)/../polyester/bin/run_polyester.R --fasta_input ${batch}/*.fa --output_dir ${batch}/ --num_samples $NUM_SAMPLES --num_replicates $NUM_REPLICATES --lib_sizes 1 --fold_changes ${batch}/*factors_*.csv --read_counts ${batch}/*.read_count.csv --bias $BIAS && touch $donefile 2>&1 > $log1" | bsub -G ms_raets -M ${mem} -n ${threads} -We 20:00 -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -J $(basename $batch) -oo $log2 
                else
                    echo "${batch} already complete"
                fi
            done
        done
    done
done
