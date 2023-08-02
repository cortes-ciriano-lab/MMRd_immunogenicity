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
parser.add_argument('-r','--ref', required=True, help='Reference genome path - fasta')
parser.add_argument('-o', '--out_folder', default='$PWD', help='Out folder (Full path)')
parser.add_argument('-t', '--threads', default=4, help='Threads [default: 4]')
EOF

######################
# 1. Building the reference and indexes required for mapping
######################

GENOMEDIR=$OUT_FOLDER
mkdir -p $GENOMEDIR

$STAR --runMode genomeGenerate --genomeDir $GENOMEDIR --genomeFastaFiles $REF --runThreadN $THREADS


