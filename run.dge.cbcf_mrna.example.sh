#!/bin/bash
# Author: SJ Riesenfeld
# Feb 24, 2018
# Cite-seq analysis

#$ -cwd
#$ -q broad
#$ -P regevlab
#$ -N DGE_ADT
#$ -l h_vmem=48g
#$ -l h_rt=36:00:00
#$ -l h='!hw-uger-*'
#$ -e log_dge_adt/
#$ -o log_dge_adt/

# Load required software
source /broad/software/scripts/useuse ## I don't understand what this code does
reuse UGER
reuse Picard-Tools
reuse .java-jdk-1.8.0_92-x86-64

DS_SCRIPTS_DIR="/ahg/regevdata/projects/CiteSeq/scripts_DS/Drop-seq_tools-1.12/"
SAMPLE_TAG="CITE_seq_ADT_S1" # used for naming and finding files

## PER TASK:

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

## Parameter settings for Drop-seq Tools DigitalExpression
CBC_MRNA_F="CBC_mRNA/barcodes_mm10_hg19_combined.no_GEM.${SAMPLE_TAG}.tsv"
INPUT_BAM="BAM/${SAMPLE_TAG}.sjrAligned_GE_tagged.bam"
#
READ_MQ=0
EDIT_DIST=2
MIN_BC_READ_THRESHOLD=1
MIN_NUM_GENES_PER_CELL=1
OUTDIR="outs/"
REPORTDIR="reports/"
JVM_HEAP_SZ="8g"

# make output directories
my_mkdir "${OUTDIR}"
my_mkdir "${REPORTDIR}"

## Create output file naming tag to help keep track of parameters used
TAG="${SAMPLE_TAG}.mq_${READ_MQ}_cbc_mrna_ed${EDIT_DIST}"
DGE_F="${OUTDIR}/${TAG}.sjrAligned_GE_tagged.dge.txt.gz"
SUMM_F="${REPORTDIR}/${TAG}.sjrAligned_GE_tagged.dge.summary.txt"

# Run the DGE program

v_exe "${DS_SCRIPTS_DIR}/DigitalExpression -m ${JVM_HEAP_SZ} I=${INPUT_BAM} O=${DGE_F} SUMMARY=${SUMM_F} CELL_BC_FILE=${CBC_MRNA_F} READ_MQ=${READ_MQ} EDIT_DISTANCE=${EDIT_DIST}"

echo "Done"

