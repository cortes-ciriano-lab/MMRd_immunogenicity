#!/bin/bash

#-------------------------------------------
#   Script to run Strelka2 in WES paired samples
#-------------------------------------------

# Databases
DBSNP=/path/to/mm10/dbsnp.vcf.gz
REFGENOME=/path/to/mm10/GRCm38.p6.genome.fa

# Path to installed version of GATK
strelka2=/path/to/strelka-2.9.10.centos6_x86_64

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
out=/path/to/strelka2_results
mkdir -p $out

# Logs folder
logs=$out/logs
mkdir -p ${logs}

#----------------------
# Strelka2
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
    
    # Strelka run
    rm -f ${out2}/runWorkflow.py
    $strelka2/bin/configureStrelkaSomaticWorkflow.py \
        --normalBam $normal_bam \
        --tumorBam $tumor_bam \
        --exome \
        --callRegions $BED \
        --referenceFasta $REFGENOME \
        --runDir $out2
    
    rm -f $logs/strelka2.${tumor_id}.${normal_id}.o.log $logs/strelka2.${tumor_id}.${normal_id}.e.log
    bsub -M 5G  -n 8  -o $logs/strelka2.${tumor_id}.${normal_id}.o.log -e $logs/strelka2.${tumor_id}.${normal_id}.e.log  -J Strelka2.${tumor_id} "${out2}/runWorkflow.py  -m local -j 8 -g 15"   

done
