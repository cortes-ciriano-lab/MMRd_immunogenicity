#!/bin/bash

#--------------------------
# Script for the alignment of the RNA-seq data against mm10 ref genome
#--------------------------

# Databases
REFGENOME=/path/to//mm10/GRCm38.p6.genome.fa
GENOME_DIR=/path/to//mm10/

# Pipeline for mapping
MAPPING=/path/to/RNAseq_mapping.sh

# Path to fastq files
FASTQ=/path/to/all_fastq_files.txt

# Samples to be analysed
# These ID must match with the IDs in the fastq files
SAMPLES=/path/to//Westcott_exomes_July_2022/metadata/rnaseq_IDs.txt

# Main out folder
out=/path/to/alignment_RNAseq
mkdir -p $out

# Alignment logs folder
logs=$out/logs
mkdir -p $logs

# If genomeDir is not created for our reference genome, this must be run before STAR alignment
GENOME_GENERATE=/path/to/RNAseq_genomeGenerate.sh
run0="bash $GENOME_GENERATE -r $REFGENOME -o $GENOME_DIR -t 4"
bsub -M 50G -n 4 -J Star_genome_generate -o ${logs}/Star_genome_generate.o.log -e ${logs}/Star_genome_generate.e.log "$run0"

# Run alignment
for i in $(cut -f2 -d' ' $SAMPLES | grep -v 'Sample' | sort | uniq);do
    echo $i
    
    # Sample sub-folder
    out2=$out/$i
    mkdir -p $out2

    # Fastq subfolder
    out3=$out2/fastq_merged
    mkdir -p $out3
    
    # Reads (fastq files)
    R1=$(grep "$i" $SAMPLES | cut -f1 -d' ' | grep -f - $FASTQ | grep '_1_sequence.fastq')
    R2=$(grep "$i" $SAMPLES | cut -f1 -d' ' | grep -f - $FASTQ | grep '_2_sequence.fastq')

    # Merging all fastq files
    run1="cat $(echo $R1) | bgzip > $out3/${i}.merged.R1_sequence.fastq.gz"    
    run2="cat $(echo $R2) | bgzip > $out3/${i}.merged.R2_sequence.fastq.gz"
    
    # Alignment with STAR (v2.7.1a) with parameters as GTEx consortium, Yizhak et al and Muyas et al.
    out4=$out2/star_results
    mkdir -p $out4
    
    run3="bash $MAPPING -t 4 -r1 $out3/${i}.merged.R1_sequence.fastq.gz -r2 $out3/${i}.merged.R2_sequence.fastq.gz -g $GENOME_DIR -o $out4 -s $i"

    rm -rf ${logs}/RNAseq.align.$i.o.log ${logs}/RNAseq.align.$i.e.log
    bsub -M 50G -n 4 -J RNAseq.align.$i -o ${logs}/RNAseq.align.$i.o.log -e ${logs}/RNAseq.align.$i.e.log "$run1 && $run2 && $run3"
done





