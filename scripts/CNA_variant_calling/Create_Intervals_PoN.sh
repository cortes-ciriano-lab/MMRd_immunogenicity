#!/bin/bash

#############################################
# Pipeline to run GATK with WES parameters
# Panel of Normals
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
parser.add_argument('-normals','--normals', required=True, help='List of normal bam files')
parser.add_argument('-bed','--bed', required=True, help='WES bed file')
parser.add_argument('-bedout','--bedout', required=True, help='Regions to be ignored')
parser.add_argument('-mappability','--mappability', required=True, help='Mappability file')
parser.add_argument('-r','--refgenome', required=True, help='Reference genome - fasta')
parser.add_argument('-o', '--out_folder', default='$PWD', help='Out folder (Full path)')
parser.add_argument('-p', '--prefix', default='PanelOfNormals', help='Prefix out files')
EOF


echo "#-----------------"
echo "# Script to buid the PanelOfNormals for GATK CNV pipeline"
echo "#-----------------"
echo

echo Parameters:
echo  List of normal bam files: "$NORMALS"
echo  WES bed file: "$BED"
echo  Regions to ignore: "$BEDOUT"
echo  Mappability file: "$MAPPABILITY"
echo  Reference genome: "$REFGENOME"
echo  Out folder: "$OUT_FOLDER"
echo  Prefix of the out files: "$PREFIX"
echo 


hold_jobs=""

#------------------------
# 1. Interval calculation
#------------------------

echo "#-----------------"
echo "# Interval calculation"
echo "#-----------------"
echo

out=$OUT_FOLDER
logs=$out/logs
mkdir -p $out
mkdir -p $logs

## GENERATE BINS
run1="$GATK PreprocessIntervals \
          -R $REFGENOME \
          -L $BED \
          --bin-length 0 \
          --padding 250 \
          --exclude-intervals $BEDOUT \
          -imr OVERLAPPING_ONLY \
          -O $out/${PREFIX}.preprocessed_intervals.interval_list"

## GENERATE ANNOTATED INTERVALS
run2="$GATK AnnotateIntervals \
          -R $REFGENOME \
          -L $out/${PREFIX}.preprocessed_intervals.interval_list \
          --mappability-track $MAPPABILITY \
          --interval-merging-rule OVERLAPPING_ONLY \
          -O $out/${PREFIX}.annotated_intervals.tsv"


# Send command to cluster
echo bsub -M 5G -o ${logs}/GATK_Annotated_intervals.o.log -e ${logs}/GATK_Annotated_intervals.e.log  -J ${PREFIX}.GATK_intervals "\"$run1 && $run2\"" > $out/cmd.logs
hold_jobs="ended(${PREFIX}.GATK_intervals)"

#-----------------------------------------------------
# 2. GATK read count for each bam
#-----------------------------------------------------

echo "#-----------------"
echo "# Read count"
echo "#-----------------"
echo

out_hdf5=$out/hdf5_samples
mkdir -p $out_hdf5

# Wait bsub command
wait_cmd="'$hold_jobs'"

hdf5_samples=""
for BAM in $(cat $NORMALS); do
       
    sample=$(basename $BAM | sed 's/.bam$//g')

    run3="$GATK CollectReadCounts \
        -L $out/${PREFIX}.preprocessed_intervals.interval_list \
        -R $REFGENOME \
        -imr OVERLAPPING_ONLY \
        --minimum-mapping-quality 30 \
        -I $BAM \
        -O $out_hdf5/${sample}.hdf5"

	# Get samples files for next step
	#hdf5_samples="$hdf5_samples -I $out_hdf5/${sample}.hdf5"
    touch $out_hdf5/${sample}.hdf5

    rm -f ${logs}/GATK_RC_${sample}.o.log ${logs}/GATK_RC_${sample}.e.log
    echo bsub -M 20G -w $wait_cmd -o ${logs}/GATK_RC_${sample}.o.log -e ${logs}/GATK_RC_${sample}.e.log  -J GATK_RC_${sample} "\"$run3\"" >> $out/cmd.logs
    # Collapse jobs to hold on for next step of bsubs
    hold_jobs="$hold_jobs && ended(GATK_RC_${sample})"
done

#---------------------------
# 3. PoN for all
#---------------------------

echo "#-----------------"
echo "# Building CNV PanelOfNormals"
echo "#-----------------"
echo

export LD_PRELOAD=/nfs/research/icortes/fmuyas/anaconda3/envs/openblas_env/lib/libopenblas.so 

# Create a list of samples
hdf5_samples=""
for hdf5_file in $(ls -d $out_hdf5/*hdf5);do
    hdf5_samples="$hdf5_samples -I $hdf5_file"
done

out_pon=$out/pon
mkdir -p $out_pon

# Merge commands
run4="$GATK CreateReadCountPanelOfNormals \
    --minimum-interval-median-percentile 5.0 \
    -O $out_pon/${PREFIX}.cnvponC.pon.hdf5"
run4="$run4 $hdf5_samples"

# Waiting bsub command
wait_cmd="'$hold_jobs'"

# Command for cluster
rm -f ${logs}/GATK_CNV_PoN.o.log ${logs}/GATK_CNV_PoN.e.log
echo bsub -M 20G -w $wait_cmd -o ${logs}/GATK_CNV_PoN.o.log -e ${logs}/GATK_CNV_PoN.e.log  -J GATK_CNV_PoN "\"$run4\"" >> $out/cmd.logs


# run all commands
cat $out/cmd.logs | sh
