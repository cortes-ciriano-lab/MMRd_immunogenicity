#!/bin/bash

#--------
# Step 1 to run GATK-PureCN Copy Number calling pipeline
#--------

SCRIPT=/path/to/Create_Intervals_PoN.sh
NORMALS=/path/to/all_WES_tails.normal_bams.txt
REFGENOME=/path/to/mm10/GRCm38.p6.genome.fa
BED=/path/to/bed_files/SureSelect_mouse_v1.mm10.liftover.UCSC.no_chrM.bed
BLACKLIST=/path/to/mm10/blacklist_regions/ENCFF547MET.bed # https://www.encodeproject.org/files/ENCFF547MET/
MAPPABILITY=/path/to/mm10/mappability/k100.umap.collapsed.bed # wget https://bismap.hoffmanlab.org/raw/mm10.umap.tar.gz
out=/path/to/GATK_PoN

bash $SCRIPT -normals $NORMALS \
-r $REFGENOME \
-bed $BED \
-mappability $MAPPABILITY \
-bedout $BLACKLIST \
-o $out \
-p Westcott_July2022

