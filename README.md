# ICGC ARGO RNA-Seq Data Simulation
This repository contains code for the generation of simulated RNA-Seq reads that can be used for the benchmarking of different parts of the ICGC ARGO RNA-Seq analysis pipeline.

The code has been used to generate read samples with the following properties:

- fully expressed transcriptome with realistic expression distribution (based on empirical estimate taken from PCAWG samples)
- germline variation taken from 1000 genomes project individual(s)
- somatic variation sampled from COSMIC database
- differential transcript expression (both dependent and independent of genotype)

The addition of further simulation properties is planned.

Currently, the simulation is separated into two stages:

1. Data preparation
2. Simulation

In the following, we will provide further description for each of the steps.

### 1. Data preparation
Before the simulation can be started, several steps of data preparation are necessary. 

#### Preparing the transcript input
As the Polyester simulation procedure makes certain assumptions on the input data, it is necessary to appropriately pre-filter the annotation. The script [pre_filter_transcripts.py](01_preprocessing/pre_filter_transcripts.py), which is called by the wrapper [run_prepare_transcripts.sh](01_preprocessing/run_prepare_transcripts.sh), retains only transcripts with a minimum length of 400nt and removes identical sequences. The script accepts the followig parameters:
```
Usage: prepare_transcripts.py <transcripts.fa> <annotation.gtf>
```
In addition to filtering, the preparation wrapper also calls the script [get_gene_length.py](01_preprocessing/get_gene_length.py) to determine the gene length of all annotated genes, where gene length is defined as the number of exonic positions across all annotated transcripts.

#### Extracting quality statistics
As part of this pre-precessing, we use `bcftools` to extract a list of simple variant calls. We will use these calls to ignore the respective positions when estimating the read-error distribution. In a second step, these variants are being filtered (again using `bcftools`) to retain only confident variation. In a third step the necessary quality information is extracted from the alignment of a real RNA-Seq sample (using the script [extract_stats_from_bam.py](01_preprocessing/extract_stats_from_bam.py). All three steps are carried out by the wrapper script [run_extract_stats.sh](01_preprocessing/run_extract_stats.sh) ). After running all thre steps, the following statistics have been collected:

- mismatch distribution (over read position and quality value)
- quality value distribution over read length
- indel/distribution (over read position and quality value)

These entities are used in simulation later on. 

### 2. Simulation
All scripts for the simulation part are collected in the ``simulation`` directory. The simulation uses the program Polyester as an engine for read generation. The step of read generation, however, is embedded into several steps of pre- and post-processing, that gather the relevant data needed for personalized transcriptome simulation. In the following, we will outline each of the subsequent steps.

Before we describe the individual steps, we will briefly outline the configuration file. Each simulation run is determined by a configuration file that contains all relevant settings and paths to inputs. It also allows to provide a unique name for the simulation, so multiple simulation runs can be told apart. The configuration file is structured as a simple bash file that assigned a set of bash variables that are sourced in subsequent scripts. The config file contains different sections, describing central simulation parameters, input and output data, as well as paths to singularity images of Polyester and ICGC containers. An example config file can be found in the `config` directory.

#### 2.1 Creating the simulation setup
As a first step, the necessary setup for the simulation needs to be prepared. The wrapper script [01_run_create_simulation_setup.sh](02_simulation/01_run_create_simulation_setup.sh) parses the config file and creates a simulation setup for each of the given input genomes (via [src/01_create_simulation_setup.sh](02_simulation/src/01_create_simulation_setup.sh). The setup process contains the following steps:

- the script [src/02_simulate_personal_transcriptome.py](02_simulation/src/02_simulate_personal_transcriptome.py) is used to generate a personalized transcriptome, integrating the provided germline and somatic variation.
- the script [src/03_compute_coverage_distribution.py](02_simulation/src/03_compute_coverage_distribution.py]) simulates a gene expression distribution modeled after a given real sample. The stats for that sample have been collected in the preprocessing steps. In addition, the script also models the weights of transcript usage within a gene, following an empirical distribution collected on the GTEx cohort. Based on the provided germline and somatic variants, also allele-specific expression is modeled. A subset of the genes is subject to differential gene expression between tumor and normal. All weights are collected per condition and are then aggregated to form the final read count for each transcript. All factors and counts are reported in a joint metadata table. 
- the script [src/04_filter_expressed_transcripts.py](02_simulation/src/04_filter_expressed_transcripts.py) filters the generated transcript expressions for transcripts with expression of 0. As Polyester cannot handle empty transcripts and they would create no reads in the first place, they are removed from the setup.
- the script [src/05_create_batches.py](02_simulation/src/05_create_batches.py) distributes the transcripts into uniform batches for parallel simulation. This is especially helpful when simulation whole transcriptome settings.

#### 2.2 Read simulation
For each of the batches and conditions collected in the previous step, the script [02_run_simulation.sh](02_simulation/02_run_simulation.sh) calls the polyester run script [run_polyester.R](polyester/bin/run_polyester.R), which has the following signature.
```
./polyester/bin/run_polyester.R --fasta_input polyester/test_inputs/five_transcripts.fa --bias rnaf --output_dir polyester/test_simulation_outputs --num_samples <N> --num_replicates <R>
```
##### Arguments:
- `--fasta_input`: fasta input
- `--output_dir`: output directory (default=./simulation_outputs)
- `--num_samples`: number of samples (default: 50)
- `--num_replicates`: number of replicates (default: 3)
- `--lib_sizes`: factor for multiplying for each sample (default: 1)
- `--fold_changes`: matrix with `ncol==num_samples` and `nrow==number of all transcripts in fasta` (will create a txt input option soon)
- `--bias`: 'none', 'rnaf', or 'cndaf'
- `--frag_GC_bias`: A txt file storing the matrix (will create a txt input option soon)
- `--gcbias`: numeric vector ranging from 0 to 7 of length sum(num_rep) (please don't use this because Polyester got a [bug](https://github.com/GoekeLab/icgc_argo_rna_data_simulation/blob/8003355afd3a63ed298bbfc7b3c2f7a2c2af9991/bin/run_polyester.R#L36) with this)

#### 2.3 Batch collection
After successful simulation, the results of the individual simulation batches are aggregated using the script [03_collect_batches.sh](02_simulation/src/03_collect_batches.sh). As a result, a single fasta file per sample and condition remains. 

#### 2.4 Fasta to fastq conversion
As Polyester only simulates reads without quality values, in this step the error distribution estimated in the pre-processing step is used together with the simulated fasta files to simulate errors into the reads and assign appropriate quality strings. Thereby the quality strings are simulated first and the errors are introduced based on the assigned quality value per position. The conversion step is invoked via the script [04_run_generate_fastq.sh](04_run_generate_fastq.sh), which in turn calls the conversion tool [src/07_fasta2fastq.py](src/07_fasta2fastq.py) for each of the fasta files.

#### 2.5 Haplotype collection
Maternal and paternal haplotypes are simulated independently. In this step. using the script [05_collect_haplotypes.sh](05_collect_haplotypes.sh), the two haplotypes are joint into a single read-file. The read files generated by this step constitute the result of the simulation.

#### 2.6 Read counting
To provide a sanity-check for simulation and as a ground-truth for evaluation, this step extracts from each simulated read its location and generates the read counts per transcript that have been generated. The wrapper script [06_run_compute_actual_readcount.sh](06_run_compute_actual_readcount.sh) calls the counter [src/08_compute_actual_readcount.py](src/08_compute_actual_readcount.py) on each of the fastq files.

#### 2.7 Metadata completion
In this last step, the script [07_run_augment_metadata.sh](07_run_augment_metadata.sh) (in turn calling [src/09_augment_metadata.py](src/09_augment_metadata.py) on each fastq file) adds the ground-truth read counts to the metadata sheet.
