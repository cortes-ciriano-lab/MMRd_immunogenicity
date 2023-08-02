#!/bin/bash

#############################################
## Pipeline for mapping RNA-seq short read data
#############################################

## Path to tools
STAR=/path/to/STAR-2.7.10a/bin/Linux_x86_64_static/STAR # Version 2.7.10a

####https://github.com/nhoffman/argparse-bash
##wget https://raw.githubusercontent.com/nhoffman/argparse-bash/master/argparse.bash
##chmod +x argparse.bash
ARGPARSE=/path/to/argparse_bash/argparse.bash
source $ARGPARSE || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-r1','--read1', required=True, help='Read 1 fastq.gz file')
parser.add_argument('-r2','--read2', required=True, help='Read 2 fastq.gz file')
parser.add_argument('-g','--genome_dir', required=True, help='Genome directory required by STAR')
parser.add_argument('-o', '--out_folder', default='$PWD', help='Out folder (Full path)')
parser.add_argument('-s', '--sample', default="Sample", help='Sample id [default: Sample]')
parser.add_argument('-t', '--threads', default=4, help='Threads [default: 4]')
EOF

echo Parameters:
echo Required read 1: "$READ1"
echo Required read 2: "$READ2"
echo Reference genome dir: "$GENOME_DIR"
echo Out folder: "$OUT_FOLDER"
echo Sample id: "$SAMPLE"
echo Threads: "$THREADS"


######################
# 1. Align fastq files
######################
OUT_BAMS=$OUT_FOLDER
mkdir -p $OUT_BAMS

cd $OUT_BAMS
$STAR --genomeDir $GENOME_DIR --runThreadN $THREADS \
    --readFilesIn $READ1 $READ2 \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNoverReadLmax 0.1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --readFilesCommand zcat \
    --alignMatesGapMax 1000000 \
    --outFilterScoreMinOverLread 0.33 \
    --outFilterMatchNminOverLread 0.33 \
    --limitSjdbInsertNsj 1200000 \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix $OUT_FOLDER/${SAMPLE}
    


