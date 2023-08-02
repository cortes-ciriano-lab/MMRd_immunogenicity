#!/bin/bash


#--------
# Step 2 to run GATK-PureCN Copy Number calling pipeline
#--------

### PureCN calling
SCRIPT=/path/to/GATK_PureCN_cnv.sh

REFGENOME=/path/to/mm10/GRCm38.p6.genome.fa
out=/path/to/Westcott_exomes_July_2022/somatic_variant_calling/copy_number/GATK_PureCN_results

INTERVALS=/path/to/GATK_PoN/Preprocessed_intervals.interval_list # Originated with Create_Intervals_PoN.sh
PON=/path/to/GATK_PoN/pon/Westcott2023.cnvponC.pon.hdf5 # Originated with Create_Intervals_PoN.sh
REF_DICT=/path/to/mm10/GRCm38.p6.genome.dict
DBSNP=/path/to/mm10_dbsnp/dbsnp.mm10.intersected.snp.vcf.gz
INTERVALS_ANNOTATED=/path/to/GATK_PoN/Westcott2023.annotated_intervals.tsv # Originated with Create_Intervals_PoN.sh

# Path to all bam files
bams=/path/to/WES_tumour_bam_files.txt

# Tsv file with known cancer purities
# If not known, they will be estimated
purities=/path/to/known_cancer_purities.tsv

for BAM in $(grep -v 'Tail' $bams);do
        sample=$(basename $BAM | cut -f1 -d'.')
        
        purity=$(grep "$sample" $purities | awk -F'\t' '{print $NF}' )

        out2=$out/${sample}
        mkdir -p $out2

        if [ -z "$purity" ]
        then
                run="bash $SCRIPT -bam $BAM \
                -r $REFGENOME \
                -r_dict $REF_DICT \
                -r_string mm10 \
                -i_ann $INTERVALS_ANNOTATED \
                -pon $PON \
                -g $DBSNP \
                -i $INTERVALS \
                -o $out2 \
                -p ${sample}"
        else
                run="bash $SCRIPT -bam $BAM \
                -r $REFGENOME \
                -r_dict $REF_DICT \
                -r_string mm10 \
                -i_ann $INTERVALS_ANNOTATED \
                -pon $PON \
                -g $DBSNP \
                -i $INTERVALS \
                --purity ${purity} \
                -o $out2 \
                -p ${sample}"
        fi

        echo $run | sh
done


