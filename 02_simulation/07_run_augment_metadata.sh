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

for fname in ${BASEDIR}/${SIM_ID}/*/*.metadata_complete.tsv
do
    if [ ! -f ${fname%tsv}augmented.tsv ]
    then
        echo processing $fname
        python $(pwd)/src/09_augment_metadata.py ${fname}
    fi
done
