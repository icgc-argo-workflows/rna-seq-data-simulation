Preprocessing scripts
=====================

The goal of these scripts is to prepare the inputs for the simulation setup and to generate additional inputs required for somulation.

- ``pre_filter_transcripts.py`` - pre-filters the transcript sequences given in fasta for a minimum length of 400nt and removes identical sequences. Both is necessary to prevent Polyester from having problems simulating. ``run_pre_filter_transcripts.sh`` is the matching run-script
- ``gen_bam_variant_calls.sh`` - given an alignment file in BAM format, we use bcftools for very simple extraction of variants. We will use these variants as ignore-positions during the estimation of quality statistics, as this would bias our estimate
- ``extract_stats_from_bam.py`` - given an alignment file in BAM format and the pre-called set of variants from the same file, we estimate quality statistics (mismatches, indels and quality strings) from the file and store the parameters to be used for simulation later on

