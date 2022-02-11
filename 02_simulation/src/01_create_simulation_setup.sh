#!/bin/bash

set -e

if [ -z "$5" ]
then
    echo "Usage: $0 <config> <1kg_genome_id> <icgc_sample_id> <random_seed> <TMB>" 
    exit 1
fi
config=$1
sample1kg=$2
icgc_sample=$3
random_seed=$4
tmb=$5

### load config - providing all capitalized variables
source ${config}

germline_vcf_pattern=${GERMLINE_VCF_DIR}/${sample1kg}/'ALL.chr*.shapeit2_integrated_v1a.GRCh38.20181129.phased.chr_mapped.'${sample1kg}.vcf.gz

### generate output directory
outdir=${BASEDIR}/${SIM_ID}/${sample1kg}_${icgc_sample}
mkdir -p $outdir

### create personal transcriptome
echo personalizing the transcriptome for $sample1kg and somatic variation
python $(pwd)/02_simulate_personal_transcriptome.py $ANNO_GTF $ANNO_FA $sample1kg $outdir $random_seed $tmb "$germline_vcf_pattern" "$SOMATIC_VCF_PATTERN"

### generate expression distribution, including read counts and fold changes
echo simulating expression distribution, read counts and fold changes
python $(pwd)/03_compute_coverage_distribution.py $GENE_COUNTS_ICGC $icgc_sample $GENE_LENS_ICGC $ANNO_GTF $ANNO_FA $sample1kg $outdir $random_seed

outbase=$(basename $ANNO_FA)
outbase=${outbase%.fa}

### we need to filter the transcript fasta files to remove non-expressed transcripts, as Polyester cannot handle them
for ht in ht1 ht2
do
    for ct in tumor normal
    do
        echo filtering for expressed transcripts $ct $ht
        python $(pwd)/04_filter_expressed_transcripts.py ${outdir}/${outbase}.${sample1kg}.${ct}_${ht}.fa ${outdir}/${outbase}.${sample1kg}.metadata_complete.tsv
    done
done

echo ""

### generate batches
for ht in ht1 ht2
do
    for ct in tumor normal
    do
        batchdir=${outdir}/batches/${ct}_${ht}
        mkdir -p $batchdir
        echo generating batches $ht $ct in $batchdir
        python $(pwd)/05_create_batches.py ${outdir}/${outbase}.${sample1kg}.${ct}_${ht}.expressed.fa ${outdir}/${outbase}.${sample1kg}.factors_${ct}_${ht}.csv ${outdir}/${outbase}.${sample1kg}.read_count.csv $BATCHSIZE $batchdir
    done
done

