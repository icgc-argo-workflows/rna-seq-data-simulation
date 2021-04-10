# icgc_argo_rna_data_simulation
Data simulation code for ICGC ARGO RNA Tech working group
### Running the script for Polyester (./bin/run_polyester.R)
```
./bin/run_polyester.R --fasta_input test_inputs/five_transcripts.fa --bias rnaf --output_dir test_simulation_outputs --num_samples 2 --num_replicates 2
```
this generates the [test_simulation_output directory](https://github.com/GoekeLab/icgc_argo_rna_data_simulation/tree/main/test_simulation_outputs)
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
