#! /bin/bash
# Author: SJ Riesenfeld
# Feb 24, 2018
# Cite-seq analysis

## EXAMPLE SCRIPT

#$ -cwd
#$ -q broad
#$ -P regevlab
#$ -l h_vmem=5g
#$ -l h_rt=12:00:00
#$ -pe smp 8
#$ -binding linear:8
#$ -e log.bclfastq2.0.1/
#$ -o log.bclfastq2.0.1/

source /broad/software/scripts/useuse
reuse .bcl2fastq2-2.17.1.14
reuse Python-2.7
reuse UGER
# export PATH=/seq/regev_genome_portal/kcoprod/galaxy/tools/external_libs/cellranger-2.0.1:$PATH
bcl2fastq --use-bases-mask=Y26,I6,Y57 --create-fastq-for-index-reads --minimum-trimmed-read-length=6 --mask-short-adapter-reads=6 --ignore-missing-positions --ignore-missing-bcls --ignore-missing-filter -p 4 -d 4 -r 4 -w 4 -R /ahg/regev_nextseq/Data03/171205_NB501583_0260_AHYH5YBGX3/ --interop-dir=/ahg/regevdata/projects/CiteSeq/citeseq_ADTlib/InterOp --output-dir=/ahg/regevdata/projects/CiteSeq/citeseq_ADTlib/ --sample-sheet=/ahg/regevdata/projects/CiteSeq/scripts_AB_SJR/samplesheet_blc2fastq.example.csv
