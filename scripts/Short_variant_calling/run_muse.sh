#!/bin/bash

#-------------------------------------------
#   Script to run Muse in WES paired samples
#-------------------------------------------


# Databases
DBSNP=/path/to/mm10/dbsnp.vcf.gz
REFGENOME=/path/to/mm10/GRCm38.p6.genome.fa

# Path to installed version of MUSE
MUSE=/path/to/MuSE.v1/MuSEv1

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
out=/path/to/Muse_results
mkdir -p $out

# Logs folder
logs=$out/logs
mkdir -p ${logs}

#----------------------
# Muse
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
    
    # Muse step 1
    run1="${MUSE} call -O $out2/$full_id \
        -l ${BED} \
        -f ${REFGENOME} \
        ${tumor_bam} \
        ${normal_bam}"

    # Muse step 2
    run2="${MUSE} sump -I $out2/$full_id.MuSE.txt \
    -E -O $out2/$full_id.MuSE.vcf \
    -D ${DBSNP}"
    
    echo $full_id $tumor_bam $normal_bam
    rm -f $logs/MUSE.${tumor_id}.${normal_id}.o.log $logs/MUSE.${tumor_id}.${normal_id}.e.log
    bsub -M 15G  -n 1  -o $logs/MUSE.${tumor_id}.${normal_id}.o.log -e $logs/MUSE.${tumor_id}.${normal_id}.e.log  -J MUSE.${tumor_id} "$run1 && $run2" 

done
