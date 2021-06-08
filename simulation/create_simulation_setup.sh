#!/bin/bash

set -e

if [ -z "$3" ]
then
    echo "Usage: $0 <1kg_genome_id> <icgc_sample_id> <random_seed>" 
    exit 1
fi
sample1kg=$1
icgc_sample=$2
random_seed=$3

BATCHSIZE=1000

anno_gtf=/cluster/work/grlab/projects/projects2020-ICGC-ARGO/simulation/annotation/gencode.v37.pc_transcripts.uniq.min400nt.gtf
transcripts=/cluster/work/grlab/projects/projects2020-ICGC-ARGO/simulation/annotation/gencode.v37.pc_transcripts.uniq.min400nt.fa
germline_vcf_pattern=/cluster/work/grlab/projects/projects2019_camda/genome/1000g_variants_hg38_split/${sample1kg}/'ALL.chr*.shapeit2_integrated_v1a.GRCh38.20181129.phased.chr_mapped.'${sample1kg}.vcf.gz
somatic_vcf_pattern='/cluster/work/grlab/projects/projects2020-ICGC-ARGO/simulation/variants/Cosmic*.vcf.gz'
outdir=/cluster/work/grlab/projects/projects2020-ICGC-ARGO/simulation/results/polyester/full_transcriptome_v1/${sample1kg}_${icgc_sample}
mkdir -p $outdir

### these files are on hg19 (from ICGC PCAWG)
gene_counts=/cluster/work/grlab/projects/ICGC/counts/counts_ICGC_combined_gene.STAR.tsv.gz
gene_lens=/cluster/work/grlab/projects/ICGC/annotation/gencode.v19.annotation.hs37d5_chr.gene_lens.tsv

### create personal transcriptome
echo personalizing the transcriptome for $sample1kg and somatic variation
python $(pwd)/simulate_personal_transcriptome.py $anno_gtf $transcripts $sample1kg $outdir $random_seed "$germline_vcf_pattern" "$somatic_vcf_pattern"

### generate expression distribution, including read counts and fold changes
echo simulating expression distribution, read counts and fold changes
python $(pwd)/compute_coverage_distribution.py $gene_counts $icgc_sample $gene_lens $anno_gtf $transcripts $sample1kg $outdir $random_seed

outbase=$(basename $transcripts)
outbase=${outbase%.fa}

### we need to filter the transcript fasta files to remove non-expressed transcripts, as Polyester cannot handle them
for ht in ht1 ht2
do
    for ct in tumor normal
    do
        echo filtering for expressed transcripts $ct $ht
        python $(pwd)/filter_expressed_transcripts.py ${outdir}/${outbase}.${sample1kg}.${ct}_${ht}.fa ${outdir}/${outbase}.${sample1kg}.metadata_complete.tsv
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
        python $(pwd)/create_batches.py ${outdir}/${outbase}.${sample1kg}.${ct}_${ht}.expressed.fa ${outdir}/${outbase}.${sample1kg}.factors_${ct}_${ht}.csv ${outdir}/${outbase}.${sample1kg}.read_count.csv $BATCHSIZE $batchdir
    done
done

