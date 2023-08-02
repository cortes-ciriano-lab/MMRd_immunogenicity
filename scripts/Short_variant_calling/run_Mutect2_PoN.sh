#!/bin/bash

#-------------------------------------------
#   Script to run Mutec2
#   Building Panel of normals
#-------------------------------------------

SCRIPT=/path/to/GATK_Mutect2_PoN.sh
REFGENOME=/path/to/mm10/GRCm38.p6.genome.fa
out=/path/to/Mutect2_results/PoN
BED=/path/to/bed_files/SureSelect_mouse_v1.mm10.liftover.UCSC.bed

# Path to bam files
BAMS=/path/to/all_WES_tails.normal_bams.txt

bash $SCRIPT \
        -normals $BAMS \
        -b $BED \
        -r $REFGENOME \
        -o $out \
        -p Westcott_WES_Mouse


