#!/bin/bash

#-------------------------------------------
#   Script to run Mutect2 in WES paired samples
#  
#-------------------------------------------

# Databases
DBSNP=/path/to/mm10/dbsnp.vcf.gz
REFGENOME=/path/to/mm10/GRCm38.p6.genome.fa

# Path to installed version of MUSE
SCRIPT=/path/to/GATK_Mutect2_Calling.sh

# Panel of Normals
PON=/path/to/Mutect2_results/PoN/final_PoN/Westcott_WES_Mouse.pon.vcf.gz

# Path to tumor-normal pairs
# First column is the tumour sample ID
# Second column is the matched normal ID
SAMPLE_PAIRS=/path/to/TumorNormal_Pairs.tsv

# Path to all bam files
# Bam files and sample IDs should match
BAMS=/path/to/Westcott_exomes_July_2022/all_bams.txt

# Bed file for the target enrichment. Lift-over of the original bed file
BED=/path/to/bed_files/SureSelect_mouse_v1.mm10.liftover.UCSC.bed

# Out folder
out=/path/to/Mutect2_results/Calls
mkdir -p $out

#----------------------
# Mutect2 Pipeline
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
    
    # Mutect2 pipeline
    run="bash ${SCRIPT} -t $tumor_bam \
        -n $normal_bam \
        -t_id $tumor_id \
        -n_id $normal_id \
        -bed $BED \
        -r $REFGENOME \
        -pon $PON \
        -o $out2"
    
    echo $run | sh

done
