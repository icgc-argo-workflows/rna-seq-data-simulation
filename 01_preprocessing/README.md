Preprocessing scripts
=====================

The goal of these scripts is to prepare the inputs for the simulation setup and to generate additional inputs required for somulation.

- ``run_extract_stats.sh`` - This script pre-processes a given alignment file to extract error statistics from the given alignment. This is a three-step process. First, a simple variant-calling procedure is run on the bam file to extract any relevant variation (to be later ignored in mismatch counting). Second, the called simple variants are filtered to retain only high-quality ones. Last, given an alignment file in BAM format and the pre-called set of variants from the same file, we estimate quality statistics (mismatches, indels and quality strings) from the file and store the parameters to be used for simulation later on. The three steps are accomplished by the following tools:
   - step 1 (simple variant calling): ``bcftools``
   - step 2 (variant filtering): ``bcftools``
   - step 3 (collection of stats): ``extract_stats_from_bam.py``
- ``run_prepare_transcripts.sh`` - This script pre-filters the transcript sequences given in fasta for a minimum length of 400nt and removes identical sequences (via the script ``pre_filter_transcripts.py``). Both is necessary to prevent Polyester from having problems simulating. In addition, this wrapper calls ``get_gene_length.py``, to collect the lengths of all annotated genes, where length is defined as the number of exon positions within the gene.

