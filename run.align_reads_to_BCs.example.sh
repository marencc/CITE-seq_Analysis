#!/bin/bash
# Author: SJ Riesenfeld
# Feb 24, 2018
# Cite-seq analysis

## TEMPLATE FOR RUNNING CITE-SEQ ALIGNMENT:

#$ -cwd
#$ -q broad
#$ -P regevlab
#$ -N sjrAlign
#$ -l h_vmem=10g
#$ -l h_rt=03:00:00
#$ -l h='!hw-uger-*'
#$ -e log_align/
#$ -o log_align/
## Submit as a job array (-t min-max) [with max 10 concurrent jobs (-tc 10)]
#$ -t 1-231 ## this * ${do_num} below should be >= the total number of seqs
## #$ -tc 20 ## ignore this - let cluster decide how many concurrent tasks

# Load required software
source /broad/software/scripts/useuse ## I don't understand what this code does
reuse UGER
reuse Python-2.7
reuse .biopython-1.64-python-2.7.1-sqlite3-rtrees
skip="0" # num sequences already processed, in case there was a previous interruption

## PER TASK:
do_num="500000" ## number of sequences a single task takes care of

#
iter=$SGE_TASK_ID
# iter="1" # for testing only
seq_num=$(($(($skip + $(($iter-1))*$do_num)) +1)) # total number of sequences to process

pyscript="/ahg/regevdata/projects/CiteSeq/scripts_ADT_SJR/align_reads_to_BCs_genlen.py"
ref="/ahg/regevdata/projects/CiteSeq/Reference/AB_BC.ADT.fa"
max_dist=2

### SET THESE PARAMETERS FOR EACH SAMPLE; THESE ARE JUST EXAMPLE SETTINGS
topdir="/ahg/regevdata/projects/CiteSeq/example_ADT_dir/"
sample_tag="CITE_seq_ADT_S1" # optional: tag to be used to label input file and output directory
input_fq="${topdir}/fastq_ADT/${sample_tag}.unaligned_mc_tagged_polyA_filtered.fastq" # input fastq file
outdir="${topdir}/BAM/${sample_tag}.sjrAlign/" # output directory name

function error_exit
{
        echo "$1" 1>&2
        exit 1
}
function my_mkdir
{
    if [ ! -d "$1" ]; then
        mkdir -p "$1" || error_exit "Cannot mkdir $1"
    fi
}
function v_exe
{
    echo "$1"
    eval "$1" || error_exit "Cannot execute command: $1"
}

my_mkdir "${outdir}"

# Run the R script
v_exe "python ${pyscript} -max_dist ${max_dist} -seq_num_start ${seq_num} -do_num ${do_num} -o ${outdir} -ref_fasta ${ref} -reads_fastq ${input_fq}"
