#!/bin/bash

#-----
# Script to merge all somatic mutation callers
# Getting concensus calls
#-----

# Databases
REFGENOME=/path/to/mm10/GRCm38.p6.genome.fa

# Path to installed version of GATK
GATK=/path/to/GATK4/gatk-4.1.8.0/gatk
BEDTOOLS=/path/to/bedtools/bedtools

# Somatic-combiner tool
SOMATITC_COMBINER=/path/to/somatic-combiner/example/somaticCombiner.jar # https://github.com/mingyi-wang/somatic-combiner

# VCF2TSV
VCF2TSV=/path/to/VcfToTsv # https://github.com/imgag/ngs-bits

# Annovar info
ANNOVAR=/path/to/annovar # https://annovar.openbioinformatics.org/en/latest/
Mousedb=/path/to/annovar/mm10db

# Bed file for the target enrichment. Lift-over of the original bed file
BED=/path/to/bed_files/SureSelect_mouse_v1.mm10.liftover.UCSC.bed

# Paths to 
MUTECT_VCFS=/path/to/mutect2_somatic_calls.txt
STRELKA_INDEL_VCFS=/path/to/all_strelka_indels_vcfs.paths.txt
STRELKA_SNV_VCFS=/path/to/all_strelka_snvs_vcfs.paths.txt
VARDICT_VCFS=/path/to/vardict_somatic_calls.txt
MUSE_VCFS=/path/to/Muse.vcf_paths.txt

out=/path/to/Merged_calls_short_variants_summary
mkdir -p $out

out1=$out/vcf_pass
mkdir -p $out1

temp=$out/individual_vcf
mkdir -p $temp

out2=$out/annovar_pass
mkdir -p $out2

tsv_folder=$out/vcf2tsv
mkdir -p $tsv_folder

logs=$out/logs
mkdir -p $logs

#-----------------
# Step 1
#-----------------

for vcf in $(cat $MUTECT_VCFS); do

    ID=$(echo $vcf | awk -F'/' '{print $(NF-1)}')
    
    mutect_file=$vcf
    strelka_INDEL_file=$(grep ${ID} $STRELKA_INDEL_VCFS)
    strelka_SNV_file=$(grep ${ID} $STRELKA_SNV_VCFS)
    vardict_file=$(grep ${ID} $VARDICT_VCFS)
    muse_file=$(grep ${ID} $MUSE_VCFS)

    echo $ID
    TUMOR=$(echo $ID | cut -f1 -d'.')
    NORMAL=$(echo $ID | cut -f3 -d'.')

    # MuSe calls
    run0="$BEDTOOLS intersect -a ${muse_file} -b $BED -header | grep 'PASS\|#' > $temp/${ID}.MuSe.vcf"

    # Left align mutect calls
    run1="$BEDTOOLS intersect -a ${mutect_file} -b $BED -header | sed 's/$TUMOR/TUMOR/g' | sed 's/$NORMAL/NORMAL/g' | grep 'PASS\|#' > $temp/${ID}.mutect2.reheader.vcf"
    run1="${run1} && $GATK IndexFeatureFile --input $temp/${ID}.mutect2.reheader.vcf"
    run1="${run1} && $GATK LeftAlignAndTrimVariants -R ${REFGENOME} -V $temp/${ID}.mutect2.reheader.vcf -O $temp/${ID}.mutect2.LeftAligned.vcf"

    # Left align strelka indels
    run2="$BEDTOOLS intersect -a $strelka_SNV_file -b $BED -header | grep 'PASS\|#' > $temp/${ID}.strelka2.SNV.vcf"
    run2="${run2} && $BEDTOOLS intersect -a $strelka_INDEL_file -b $BED -header | grep 'PASS\|#'  > $temp/${ID}.strelka2.INDEL.vcf"
    run2="${run2} && $GATK LeftAlignAndTrimVariants -R ${REFGENOME} -V $temp/${ID}.strelka2.INDEL.vcf -O $temp/${ID}.strelka2.INDEL.LeftAligned.vcf"

    # Left align VarDict calls
    run3="grep -v '<DEL>' $vardict_file | grep -v '<INV>' | grep -v '<DUP>' |  $BEDTOOLS intersect -a stdin -b $BED -header | grep 'PASS\|#' | grep 'Somatic\|#' | sed 's/$TUMOR/TUMOR/g' | sed 's/$NORMAL/NORMAL/g' > $temp/${ID}.vardict.reheader.vcf"
    run3="${run3} && $GATK LeftAlignAndTrimVariants -R ${REFGENOME} -V $temp/${ID}.vardict.reheader.vcf -O $temp/${ID}.vardict.LeftAligned.vcf"

    # Merge vcf files
    # run somatic combiner for SNVs and Indels seperately
    run4="java -jar ${SOMATITC_COMBINER} -M $temp/${ID}.mutect2.LeftAligned.vcf -u $temp/${ID}.MuSe.vcf -s $temp/${ID}.strelka2.SNV.vcf -S $temp/${ID}.strelka2.INDEL.LeftAligned.vcf -D $temp/${ID}.vardict.LeftAligned.vcf -o $out1/${ID}.mutect2_strelka2_vardict_muse.merge.unfiltered.vcf"
    run4="${run4} && grep '^#' $out1/${ID}.mutect2_strelka2_vardict_muse.merge.unfiltered.vcf > $temp/${ID}.header.vcf && \
            grep -v '^#' $out1/${ID}.mutect2_strelka2_vardict_muse.merge.unfiltered.vcf | grep 'PASS' | grep -v 'NumCallers=1' | awk '{if (length(\$4) == length(\$5)) {print \$0}}' | grep 'Mutect2' > $temp/${ID}.pass.snvs.vcf && \
            grep -v '^#' $out1/${ID}.mutect2_strelka2_vardict_muse.merge.unfiltered.vcf | grep 'PASS' | grep -v 'NumCallers=1' | awk '{if (length(\$4) != length(\$5)) {print \$0}}' > $temp/${ID}.pass.indels.vcf && \
            cat $temp/${ID}.pass.snvs.vcf $temp/${ID}.pass.indels.vcf | sort -k1,1 -k2,2n | cat $temp/${ID}.header.vcf - > $out1/${ID}.mutect2_strelka2_vardict_muse.merge.pass.vcf"

    # Remove unnecesary vcfs
    run5="rm -f $temp/${ID}*reheader* $temp/${ID}*INDEL.vcf $temp/${ID}*idx $temp/${ID}.pass.snvs.vcf $temp/${ID}.pass.indels.vcf $temp/${ID}.header.vcf"

    # Vcfs to tsv
    run6="$VCF2TSV -in $out1/${ID}.mutect2_strelka2_vardict_muse.merge.unfiltered.vcf -out $tsv_folder/${ID}.mutect2_strelka2_vardict_muse.merge.unfiltered.tsv"
    run6="${run6} && $VCF2TSV -in $out1/${ID}.mutect2_strelka2_vardict_muse.merge.pass.vcf -out $tsv_folder/${ID}.mutect2_strelka2_vardict_muse.merge.pass.tsv"

    # Annovar
    run7="perl $ANNOVAR/convert2annovar.pl -format vcf4old \
         $out1/${ID}.mutect2_strelka2_vardict_muse.merge.pass.vcf > $out2/${ID}.pass.avinput"
    
    run8="perl $ANNOVAR/table_annovar.pl $out2/${ID}.pass.avinput \
        $Mousedb --buildver mm10 \
        -out $out2/${ID}.pass.annovar \
        -remove -protocol refGene,cytoBand -operation g,r \
        -nastring . -csvout -polish"

    rm -f ${logs}/${ID}.annovar.o.log ${logs}/${ID}.annovar.e.log
    bsub -M 10G -o ${logs}/${ID}.annovar.o.log -e ${logs}/${ID}.annovar.e.log -J ${ID}.annovar "$run0 && $run1 && $run2 && $run3 && $run4 && $run5 && $run6 && $run7 && $run8"

done



