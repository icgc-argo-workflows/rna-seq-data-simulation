#!/usr/bin/env Rscript

## LOAD LIBRARIES  
library(optparse)
library(Biostrings)
library(polyester)

## PARSE COMMAND-LINE PARAMETERS 
option_list <- list(
  make_option(c("-i", "--fasta_input"), type="character", default=NULL, metavar="path", help="fasta input"),
  make_option(c("-o", "--output_dir"), type="character", default="./simulation_outputs", metavar="path", help="output directory"),
  make_option(c("-s", "--num_samples"), type="integer"  , default=50, metavar="integer", help=""),
  make_option(c("-r", "--num_replicates"), type="integer"  , default=3, metavar="integer", help=""),
  make_option(c("-l", "--lib_sizes" ), type="integer", default=1, metavar="integer" , help=""),
  make_option(c("-f", "--fold_changes"), type="matrix", default=NULL, metavar="path", help=""),
  make_option(c("-b", "--bias"), type="character", default='none', metavar="string" , help=""),
  # make_option(c("-g", "--gcbias"), type="numeric", default=FALSE, metavar="numeric", help="numeric vector ranging from 0 to 7 of length sum(num_rep)"),
  make_option(c("-j", "--frag_GC_bias"), type="character", default='none', metavar="path", help="A txt file storing the matrix")
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

## PREPARE PARAMETERS FOR simulate_experiment()
fasta_input <- opt$fasta_input
output_dir <- opt$output_dir
num_transcripts <- length(readDNAStringSet(fasta_input))
num_samples <- opt$num_samples
num_replicates <- opt$num_replicates
num_reps <- c(rep(num_replicates,num_samples))
size <- NULL ##set to default
lib_sizes <- c(rep(opt$lib_sizes,sum(num_reps)))
fold_changes <- matrix(sample(1:4,sum(c(rep(num_transcripts,num_samples))),replace=T),nrow=num_transcripts) ##nrow corresponds to the number of transcripts
bias <- opt$bias
# gcbias<-as.numeric(sample(0:7,sum(num_reps),replace = TRUE)) ##Polyester got a bug that doesn't run even the input is numeric
# for some reason gcbias cannot be recognized as numeric: https://github.com/alyssafrazee/polyester/blob/29263d1a15ee7e33adef5d2d6aeeb93a6cb73d92/R/simulate_experiment.R#L422
frag_GC_bias <- opt$frag_GC_bias

## RUN simulate_experiment()
if (frag_GC_bias == 'none'){
  simulate_experiment(fasta=fasta_input,
                      outdir=output_dir,
                      num_reps=num_reps,
                      meanmodel=TRUE,
                      paired=TRUE,
                      seed=123,
                      reportCoverage = TRUE, 
                      readlen = 100, #empirical error_models can only accept length < 102
                      error_model='illumina5',
                      distr='normal',
                      ##variable parameters
                      # gcbias=gcbias,
                      size=size, #currently NULL -> fixed as default
                      fold_changes=fold_changes,
                      lib_sizes=lib_sizes,
                      bias=bias)
}else{ ##frag_GC_bias='none' throws an error
  simulate_experiment(fasta=fasta_input,
                      outdir=output_dir,
                      num_reps=num_reps,
                      meanmodel=TRUE,
                      paired=TRUE,
                      seed=123,
                      reportCoverage = TRUE, 
                      readlen = 100, #empirical error_models can only accept length < 102
                      error_model='illumina5',
                      distr='normal',
                      ##variable parameters
                      # gcbias=gcbias,
                      size=size, #currently NULL -> fixed as default
                      fold_changes=fold_changes,
                      lib_sizes=lib_sizes,
                      bias=bias,
                      frag_GC_bias=frag_GC_bias)
}



