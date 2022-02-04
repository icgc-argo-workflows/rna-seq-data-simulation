#!/bin/bash

set -e
genome=/cluster/work/grlab/projects/ICGC/genomes/hg19_hs37d5/genome.fa
#bamfile=/cluster/work/grlab/projects/ICGC/alignments/alignments_ICGC_2015-07-15/035824ea-ffa5-11e4-9160-83b89f3ab522.aligned.bam
bamfile=/cluster/work/grlab/projects/ICGC/alignments/alignments_ICGC_2015-07-15/05caea96-dd31-11e4-8808-f884c254ba06.aligned.bam
bambase=$(basename $bamfile)
outdir=/cluster/work/grlab/projects/projects2020-ICGC-ARGO/simulation/annotation
outfile_stats=${outdir}/${bambase%.bam}.stats.tsv
outfile_vcf=${outdir}/${bambase%.bam}.calls.vcf

threads=8

. "/cluster/home/akahles/anaconda3/etc/profile.d/conda.sh"
conda activate bcftools

if [ ! -f ${outfile_vcf} ]
then
    echo "extracting variant positions from $bamfile"
    bcftools mpileup --ff 256 -f ${genome} --threads $threads -a AD,DP --ignore-RG ${bamfile} | bcftools call -mv --threads $threads > ${outfile_vcf}
fi

if [ ! -f ${outfile_vcf%.vcf}.Q40.vcf ]
then
    echo "filtering variants for acceptable quality"
    bcftools filter -i 'QUAL>40' ${outfile_vcf} > ${outfile_vcf%.vcf}.Q40.vcf
fi
conda deactivate

echo "generating stats for alignment file $bamfile"
python $(pwd)/extract_stats_from_bam.py $bamfile ${outfile_vcf%.vcf}.Q40.vcf $outfile_stats

