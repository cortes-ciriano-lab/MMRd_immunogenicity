#!/bin/bash

#--------------------------
# Script for the alignment of the WES data against mm10 ref genome
#--------------------------

# Tools
PICARD=/nfs/research/icortes/SOFTWARE/picard-2.25.4/picard.jar
SAMTOOLS=/nfs/research/icortes/SOFTWARE/samtools-1.11/bin/samtools

# Samples to be analysed
BAMS=/nfs/research/icortes/DATA/Westcott_exomes_July_2022/alignment_RNAseq/rnaseq.raw_bams.txt

# Main out folder
out=/nfs/research/icortes/DATA/Westcott_exomes_July_2022/alignment_RNAseq/
mkdir -p $out

# Alignment logs folder
logs=$out/logs
mkdir -p $logs

# Run alignment
for bam in $(cat $BAMS);do
   
    # Sample sub-folder
    out2=$(dirname $bam)
    mkdir -p $out2
    
    id=$(basename $bam | sed 's/\.bam//')
    
    # Step 1. AddOrReplaceReadGroups
    run1="java -jar $PICARD AddOrReplaceReadGroups I=$bam \
        O=$out2/${id}.RG.sorted.bam \
        SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample \
        && $SAMTOOLS index $out2/${id}.RG.sorted.bam"
    
    rm -rf ${logs}/RNAseq.process.$id.o.log ${logs}/RNAseq.process.$id.e.log
    bsub -M 30G -J RNAseq.process.$id -o ${logs}/RNAseq.process.$id.o.log -e ${logs}/RNAseq.process.$id.e.log "$run1"
done





