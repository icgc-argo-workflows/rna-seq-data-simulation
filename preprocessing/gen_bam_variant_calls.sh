#!/bin/bash

set -e

. "/cluster/home/akahles/anaconda3/etc/profile.d/conda.sh"

conda activate bcftools

 bcftools mpileup --ff 256 -f ../../genomes/hg19_hs37d5/genome.fa -a AD,DP ~/test.bam | bcftools call -mv > ~/test.vcf
