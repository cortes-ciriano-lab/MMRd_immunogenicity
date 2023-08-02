#!/bin/bash

#-------------------------------------------
# Script to get the mutation base counts in the RNA-seq data
#-------------------------------------------

# Databases
REFGENOME=/path/to/mm10/GRCm38.p6.genome.fa

# Python script
SCRIPT=/path/to/ExtractBaseCountsFromBam.py

# Path to normal bam files
BAMS=/path/to/rnaseq_bams.txt

# Final vcf files
VCFS=/path/to/final_concensus_vcfs.txt

# Out folder
out=/path/to/RNAseq_base_counts
mkdir -p $out

# Logs folder
logs=$out/logs
mkdir -p ${logs}

#----------------------
# Mutect2
#----------------------
IFS=$'\n'
for bam in $(cat $BAMS); do
    
    #-------------------------------------------
    #   Get base counts for all bases
    #-------------------------------------------
    ID=$(basename $bam | cut -f1 -d'.') # For a bam file like: Sample_id.rnaseq.bam
    VCF=$(grep ${ID} $VCFS) # Getting the VCF to check. The ID must match between RNA-seq and WES-vcf data
    
    if [ ! -z "$VCF" ]
    then
        ID2=$(basename $VCF | cut -f1 -d'.')
        out2=$out/${ID2}
        mkdir -p $out2
        mkdir -p $out2/temp
        
        # Get the genomic coordinates of the PASS mutations in the concensus vcf (somatic variants)
        grep 'PASS' $VCF | grep -v '^#' | awk -F'\t' -v OFS='\t' '{print $1, $2-1, $2}' | sort | uniq | sort -k1,1 -k2,2n > $out2/${ID2}.${ID}.somatic.bed
        
        # Command line to run the python script that extracts the base counts
        run="python $SCRIPT --bam $bam \
        --ref $REFGENOME \
        --nprocs 8 \
        --out_file $out2/${ID2}.${ID}.base_counts.rnaseq.tsv \
        --bin 200000 \
        --bed $out2/${ID2}.${ID}.somatic.bed \
        --tmp_dir $out2/temp"
    
        rm -f ${logs}/${ID}.BaseCounts_RNA.o.log ${logs}/${ID}.BaseCounts_RNA.e.log     
        bsub -M 5G -n 8 -o ${logs}/${ID}.BaseCounts_RNA.o.log -e ${logs}/${ID}.BaseCounts_RNA.e.log  -J BaseCounts_RNA_${ID} "$run && rm -rf $out2/temp"
    fi
    
done
