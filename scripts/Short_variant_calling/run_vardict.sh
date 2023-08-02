#!/bin/bash

#-------------------------------------------
#   Script to run Vardict in WES paired samples
#   https://github.com/AstraZeneca-NGS/VarDict
#-------------------------------------------

# Steps described in:
# https://github.com/AstraZeneca-NGS/VarDict
VARDICT=/path/to/VarDict-1.8.2/bin

# Databases
DBSNP=/path/to/mm10/dbsnp.vcf.gz
REFGENOME=/path/to/mm10/GRCm38.p6.genome.fa

# Path to tumor-normal pairs
# First column is the tumour sample ID
# Second column is the matched normal ID
SAMPLE_PAIRS=/path/to/TumorNormal_Pairs.tsv

# Path to all bam files
# Bam files and sample IDs should match
BAMS=/path/to/Westcott_exomes_July_2022/all_bams.txt

# Bed file for the target enrichment. Lift-over of the original bed file
BED=/path/to/bed_files/SureSelect_mouse_v1.mm10.liftover.UCSC.bed

# Out folder for PoNs
out=/path/to/VarDict_results
mkdir -p $out

# Logs folder
logs=$out/logs
mkdir -p ${logs}

#----------------------
# VarDict
#----------------------
IFS=$'\n'
for pair in $(grep -v '^Tumor' $SAMPLE_PAIRS | sort | uniq); do
    
    #-------------------------------------------
    #   Get parameters and path to files
    #-------------------------------------------
    
    # Tumor id
    tumor_id=$(echo $pair | cut -f1)
    tumor_bam=$(grep "$tumor_id\." $BAMS)
    
    # Normal id
    normal_id=$(echo $pair | cut -f2)
    normal_bam=$(grep "$normal_id\." $BAMS)
    
    # Out folder
    full_id=$tumor_id.vs.$normal_id
    out2=$out/$full_id
    mkdir -p $out2
    
    # Vardict run
    AF_THR="0.005" # minimum allele frequency
    run="$VARDICT/VarDict -G $REFGENOME -f $AF_THR -N $tumor_id -b \"$tumor_bam|$normal_bam\" -c 1 -S 2 -E 3 -g 4 $BED | $VARDICT/testsomatic.R | $VARDICT/var2vcf_paired.pl -N \"$tumor_id|$normal_id\" -f $AF_THR > $out2/$full_id.vardict.vcf"

    echo $full_id
    rm -f $logs/vardict.${tumor_id}.${normal_id}.o.log $logs/vardict.${tumor_id}.${normal_id}.e.log
    echo "$run"
    bsub -M 10G -o $logs/vardict.${tumor_id}.${normal_id}.o.log -e $logs/vardict.${tumor_id}.${normal_id}.e.log  -J vardict.${tumor_id} "cd $out2 && $run"   
    echo

done
