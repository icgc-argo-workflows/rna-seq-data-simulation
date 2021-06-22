# ICGC ARGO RNA-Seq Data Simulation
This repository contains code for the generation of simulated RNA-Seq reads that can be used for the benchmarking of different parts of the ICGC ARGO RNA-Seq analysis pipeline.

The code has been used to generate read samples with the following properties:

- fully expressed transcriptome with realistic expression distribution (based on empirical estimate taken from PCAWG samples)
- germline variation taken from 1000 genomes project individual(s)
- somatic variation sampled from COSMIC database
- differential transcript expression (both dependent and independent of genotype)

The addition of further simulation properties is planned.

### Data preparation
Before the simulation can be started, several steps of data preparation are necessary. 

#### Pre-filtering the transcript input
As the Polyester simulation procedure makes certain assumptions on the input data, it is necessary to appropriately pre-filter the annotation. The script [pre_filter_transcripts.py](https://github.com/icgc-argo-rna-wg/data-simulation/blob/main/preprocessing/pre_filter_transcripts.py), which is called by the wrapped [run_pre_filter_transcripts.sh](https://github.com/icgc-argo-rna-wg/data-simulation/blob/main/preprocessing/run_pre_filter_transcripts.sh), retains only transcripts with a minimum length of 400nt and removes identical sequences. The script accepts the followig parameters:
```
Usage: pre_filter_transcripts.py <transcripts.fa> <annotation.gtf>
```

### Generating simple variant calls
The script [gen_bam_variant_calls.sh](https://github.com/icgc-argo-rna-wg/data-simulation/blob/main/preprocessing/gen_bam_variant_calls.sh) accepts a BAM file as input and uses `bcftools` to extract a list of simple variant calls. We will use these calls to ignore the respective positions when estimating the read-error distribution. 

#### Extracting quality statistics
In a third preprocessing step the necessary quality information is extracted from the alignment of a real RNA-Seq sample. The script [extract_stats_from_bam.py](https://github.com/icgc-argo-rna-wg/data-simulation/blob/main/preprocessing/extract_stats_from_bam.py) (invoked by the wrapped [run_extract_stats.sh](https://github.com/icgc-argo-rna-wg/data-simulation/blob/main/preprocessing/run_extract_stats.sh) ), takes a BAM file as well as the above variant calls as input and generates the following statistics:

- mismatch distribution (over read position and quality value)
- quality value distribution over read length
- indel/distribution (over read position and quality value)

These enitities are used in simulation later on. 

### Running the simulation
This simulation uses the program Polyester as an engine for read generation. The wrapper for using Polyester can be easily tested with the folliwing command. 

```
./polyester/bin/run_polyester.R --fasta_input polyester/test_inputs/five_transcripts.fa --bias rnaf --output_dir polyester/test_simulation_outputs --num_samples 2 --num_replicates 2
```
The above command generates the [test_simulation_output directory](https://github.com/GoekeLab/icgc_argo_rna_data_simulation/tree/main/polyester/
test_simulation_outputs)
#### Arguments include:
- `--fasta_input`: fasta input
- `--output_dir`: output directory (default=./simulation_outputs)
- `--num_samples`: number of samples (default: 50)
- `--num_replicates`: number of replicates (default: 3)
- `--lib_sizes`: factor for multiplying for each sample (default: 1)
- `--fold_changes`: matrix with `ncol==num_samples` and `nrow==number of all transcripts in fasta` (will create a txt input option soon)
- `--bias`: 'none', 'rnaf', or 'cndaf'
- `--frag_GC_bias`: A txt file storing the matrix (will create a txt input option soon)
- `--gcbias`: numeric vector ranging from 0 to 7 of length sum(num_rep) (please don't use this because Polyester got a [bug](https://github.com/GoekeLab/icgc_argo_rna_data_simulation/blob/8003355afd3a63ed298bbfc7b3c2f7a2c2af9991/bin/run_polyester.R#L36) with this)
