#!/bin/bash

#############################################
# Pipeline to run GATK with WES parameters
# GATK variant calling + PureCN
# conda activate /nfs/research/icortes/fmuyas/anaconda3/envs/PureCN_v2.0.1
#############################################

## Path to tools
# Path to installed version of GATK
GATK=/path/to/GATK4/gatk-4.1.8.0/gatk

####https://github.com/nhoffman/argparse-bash
##wget https://raw.githubusercontent.com/nhoffman/argparse-bash/master/argparse.bash
##chmod +x argparse.bash
ARGPARSE=/path/to/argparse_bash/argparse.bash
source $ARGPARSE || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-bam','--bam', required=True, help='Tumour bam file')
parser.add_argument('-i','--intervals', required=True, help='Interval file obtained by gatk when running PoN')
parser.add_argument('-i_ann','--annotated_intervals', required=True, help='Annotated interval file obtained by gatk when running PoN')
parser.add_argument('-pon','--pon', required=True, help='Panel of normals')
parser.add_argument('-g','--gnomad', required=True, help='Gnomad vcf file, SNPs in gnomad overlapping with the exome kit')
parser.add_argument('-r','--refgenome', required=True, help='Reference genome - fasta')
parser.add_argument('-r_dict','--refgenome_dict', required=True, help='Reference genome - dict')
parser.add_argument('-r_string','--refgenome_string', default='hg38', help='Reference genome - string [hg38,hg19]')
parser.add_argument('-o', '--out_folder', default='$PWD', help='Out folder (Full path)')
parser.add_argument('-p', '--prefix', default='PanelOfNormals', help='Prefix out files')
parser.add_argument('-purity', '--purity', default='Unknown', required=False, help='Expected purity. Do not specify value if unknown')

EOF


echo "#-----------------"
echo "# Script to run GATK CNV pipeline"
echo "#-----------------"
echo

echo Parameters:
echo  List of normal bam files: "$NORMALS"
echo  Interval file: "$INTERVALS"
echo  Annotated intervals: "$ANNOTATED_INTERVALS"
echo  Panel of normals: "$PON"
echo  Gnomad vcf file: "$GNOMAD"
echo  Reference genome: "$REFGENOME"
echo  Out folder: "$OUT_FOLDER"
echo  Prefix of the out files: "$PREFIX"
echo  Expected purity: "$PURITY"


hold_jobs=""

# Create out folder
out=$OUT_FOLDER
logs=$out/logs
mkdir -p $out
mkdir -p $logs

sample=$PREFIX

#-----------------------------------------------------
# 1. GATK read count for each bam
#-----------------------------------------------------

echo "#-----------------"
echo "# Read count"
echo "#-----------------"
echo

out_hdf5=$out/hdf5_samples
mkdir -p $out_hdf5

# Wait bsub command
wait_cmd="'$hold_jobs'"


run1="$GATK CollectReadCounts \
    -L ${INTERVALS} \
    -R $REFGENOME \
    -imr OVERLAPPING_ONLY \
    --minimum-mapping-quality 30 \
    -I $BAM \
    -O $out_hdf5/${sample}.hdf5"

# Get samples files for next step
touch $out_hdf5/${sample}.hdf5

rm -f ${logs}/GATK_RC_${sample}.o.log ${logs}/GATK_RC_${sample}.e.log
echo bsub -M 20G -o ${logs}/GATK_RC_${sample}.o.log -e ${logs}/GATK_RC_${sample}.e.log  -J GATK_RC_${sample} "\"$run1\"" > $out/cmd.logs

hold_jobs="ended(GATK_RC_${sample})"

#-----------------------------------------------------
# 2. GATK read count denoise + ploting
#-----------------------------------------------------

out_denoised=$out/denoised
mkdir -p $out_denoised


# Waiting bsub command
wait_cmd="'$hold_jobs'"

# Denoise
run2="$GATK DenoiseReadCounts \
    -I $out_hdf5/${sample}.hdf5 \
    --count-panel-of-normals $PON \
    --annotated-intervals $ANNOTATED_INTERVALS \
    --standardized-copy-ratios $out_denoised/${sample}.standardizedCR.tsv \
    --denoised-copy-ratios $out_denoised/${sample}.denoisedCR.tsv"

run3="$GATK PlotDenoisedCopyRatios \
    --standardized-copy-ratios $out_denoised/${sample}.standardizedCR.tsv \
    --denoised-copy-ratios $out_denoised/${sample}.denoisedCR.tsv \
    --sequence-dictionary $REFGENOME_DICT \
    --minimum-contig-length 50000000 \
    --output $out_denoised \
    --output-prefix ${sample}"

rm -f ${logs}/GATK_Denoise_${sample}.o.log ${logs}/GATK_Denoise_${sample}.e.log   
echo bsub -M 5G -w $wait_cmd -o ${logs}/GATK_Denoise_${sample}.o.log -e ${logs}/GATK_Denoise_${sample}.e.log  -J GATK_RCden_${sample} "\"$run2 && $run3\"" >> $out/cmd.logs

hold_jobs="$hold_jobs && ended(GATK_RCden_${sample})"


#-----------------------------------------------------
# 3. GATK CollectAllelicCounts from bam
#-----------------------------------------------------
out_AC=$out/AllelicCounts
mkdir -p $out_AC

# Waiting bsub command
wait_cmd="'$hold_jobs'"

run4="$GATK CollectAllelicCounts \
    -L $GNOMAD \
    -I $BAM \
    -R $REFGENOME \
    -O $out_AC/${sample}.AllelicCounts.tsv"

rm -f ${logs}/GATK_CollectAllelicCounts.${sample}.o.log ${logs}/GATK_CollectAllelicCounts.${sample}.e.log
echo bsub -M 30G -w $wait_cmd -o ${logs}/GATK_CollectAllelicCounts.${sample}.o.log -e ${logs}/GATK_CollectAllelicCounts.${sample}.e.log  -J GATK_AC_${sample} "\"$run4\"" >> $out/cmd.logs

hold_jobs="$hold_jobs && ended(GATK_AC_${sample})"


#-----------------------------------------------------
# 4. Segmentation
#-----------------------------------------------------

# Call segments
out_segments=$out/segments
mkdir -p $out_segments

# Waiting bsub command
wait_cmd="'$hold_jobs'"

# To get smoother results
# --number-of-changepoints-penalty-factor 5.0
run5="$GATK ModelSegments \
    --denoised-copy-ratios $out_denoised/${sample}.denoisedCR.tsv \
    --minimum-total-allele-count-normal 20 \
    --number-of-changepoints-penalty-factor 5.0 \
    --allelic-counts $out_AC/${sample}.AllelicCounts.tsv \
    --output $out_segments \
    --output-prefix ${sample}"

run6="$GATK CallCopyRatioSegments \
     --input $out_segments/${sample}.cr.seg \
     --output $out_segments/${sample}.clean.called.seg"

run7="$GATK PlotModeledSegments \
    --denoised-copy-ratios $out_denoised/${sample}.denoisedCR.tsv \
    --allelic-counts $out_segments/${sample}.hets.tsv \
    --segments $out_segments/${sample}.modelFinal.seg \
    --sequence-dictionary $REFGENOME_DICT \
    --minimum-contig-length 50000000 \
    --output $out_segments \
    --output-prefix ${sample}"

rm -f ${logs}/GATK_Segment.${sample}.o.log ${logs}/GATK_Segment.${sample}.e.log
echo bsub -M 5G -w $wait_cmd -o ${logs}/GATK_Segment.${sample}.o.log -e ${logs}/GATK_Segment.${sample}.e.log  -J GATK_Segments_${sample} "\"$run5 && $run6 && $run7\"" >> $out/cmd.logs

#-----------------
# 5. PureCN
# Step 1
# Necessary before running step 2: Preparing intervals for PureCN
#-----------------

# Export script from PureCN
export PURECN=/nfs/research/icortes/fmuyas/anaconda3/envs/PureCN_v2.0.1/lib/R/library/PureCN/extdata
R_LIBS=/nfs/research/icortes/fmuyas/anaconda3/envs/PureCN_v2.0.1/lib/R/library

# Create PureCN output folder
out2=$out/PureCN_results
mkdir -p $out2

# Create interval files
out_ref=$out2/reference_files
mkdir -p $out_ref

grep -v '^@' $INTERVALS | awk -F'\t' -v OFS='\t' '{print $1,$2,$3}' > $out_ref/original_gatk_intervals.bed

run8="Rscript $PURECN/IntervalFile.R --infile $out_ref/original_gatk_intervals.bed \
    --fasta $REFGENOME --outfile $out_ref/PureCN_intervals.${REFGENOME_STRING}.txt \
    --offtarget --genome $REFGENOME_STRING \
    --export $out_ref/PureCN_intervals.${REFGENOME_STRING}.optimized.bed"

# Job to submit
rm -rf $logs/Intervals.o.log $logs/Intervals.e.log
echo bsub -M 5G -o $logs/Intervals.o.log -e $logs/Intervals.e.log -J PureCN_intervals.${sample} "\"$run8\"" >> $out/cmd.logs

#-----------------------------------------------------
# 6. PureCN
# Step 2
# PureCN calling using the GATK CNV calls
#-----------------------------------------------------
wait_cmd="'ended(GATK_Segments_${sample}) && ended(PureCN_intervals.${sample})'"

# Check if we know purity
if [[ $PURITY == 'Unknown' ]]; then
    run9="Rscript $PURECN/PureCN.R --out $out2/${sample}  \
        --sampleid ${sample} \
        --tumor $out_hdf5/${sample}.hdf5 \
        --log-ratio-file $out_denoised/${sample}.denoisedCR.tsv \
        --seg-file $out_segments/${sample}.modelFinal.seg \
        --vcf $out_AC/${sample}.AllelicCounts.tsv \
        --sex ? \
        --min-af 0.1 \
        --intervals $out_ref/PureCN_intervals.${REFGENOME_STRING}.txt \
        --genome ${REFGENOME_STRING} \
        --fun-segmentation Hclust \
        --force --post-optimize --seed 123"
else

    # Get the purity range
    min_purity=$(python -c "print(round($PURITY-0.1,3))")
    max_purity=$(python -c "print(round($PURITY+0.1,3))")

    run9="Rscript $PURECN/PureCN.R --out $out2/${sample}  \
        --sampleid ${sample} \
        --tumor $out_hdf5/${sample}.hdf5 \
        --log-ratio-file $out_denoised/${sample}.denoisedCR.tsv \
        --seg-file $out_segments/${sample}.modelFinal.seg \
        --vcf $out_AC/${sample}.AllelicCounts.tsv \
        --sex ? \
        --min-purity ${min_purity} \
        --max-purity ${max_purity} \
        --min-af 0.1 \
        --intervals $out_ref/PureCN_intervals.${REFGENOME_STRING}.txt \
        --genome ${REFGENOME_STRING} \
        --fun-segmentation Hclust \
        --force --post-optimize --seed 123"
fi

# Job to submit
rm -rf $logs/PureCN.o.log $logs/PureCN.e.log 
echo bsub -M 10G -w $wait_cmd -o $logs/PureCN.o.log -e $logs/PureCN.e.log -J PureCN_calling.${sample} "\"$run9\"" >> $out/cmd.logs

# Run command
cat $out/cmd.logs | sh 