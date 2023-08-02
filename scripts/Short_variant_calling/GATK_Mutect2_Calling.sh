#!/bin/bash

#############################################
# Pipeline to run GATK Mutect2
# Panel of Normals
#############################################

## Path to tools
# Path to installed version of GATK
GATK=/path/to/GATK4/gatk-4.2.0.0/gatk

####https://github.com/nhoffman/argparse-bash
##wget https://raw.githubusercontent.com/nhoffman/argparse-bash/master/argparse.bash
##chmod +x argparse.bash
ARGPARSE=/path/to/argparse_bash/argparse.bash
source $ARGPARSE || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-t','--tumor_bam', required=True, help='Tumor bam file')
parser.add_argument('-n','--normal_bam', required=True, help='Normal bam file')
parser.add_argument('-t_id','--tumor_id', required=True, help='Tumor id (must match the one found in the bam file)')
parser.add_argument('-n_id','--normal_id', required=True, help='Normal id (must match the one found in the bam file)')
parser.add_argument('-bed','--bed', required=True, help='WES bed file')
parser.add_argument('-r','--refgenome', required=True, help='Reference genome - fasta')
parser.add_argument('-pon', '--pon', required=True, help='Panel of Normals')
parser.add_argument('-o', '--out_folder', default='$PWD', help='Out folder (Full path)')
parser.add_argument('-g', '--gnomad', default='Unknown', help='Gnomad vcf with population allele frequencies')
EOF


echo "#-----------------"
echo "# Script to buid the PanelOfNormals for GATK CNV pipeline"
echo "#-----------------"
echo

echo Parameters:
echo  Tumor bam file: "$TUMOR_BAM"
echo  Normal bam file: "$NORMAL_BAM"
echo  Tumor ID : "$TUMOR_ID"
echo  Normal ID: "$NORMAL_ID"
echo  WES bed file: "$BED"
echo  Panel of Normals: "$PON"
echo  Reference genome: "$REFGENOME"
echo  Gnomad vcf file: "$GNOMAD"
echo  Out folder: "$OUT_FOLDER"
echo 


# Out folder
out=$OUT_FOLDER
mkdir -p $out

# Logs folder
logs=$out/logs
mkdir -p ${logs}

# Restart commands
rm -rf $out/cmd.logs
touch $out/cmd.logs

## Now it must be passed as arguments to avoid problems
# # Getting IDs
# TUMOR_ID=$(basename $TUMOR_BAM | sed 's/.bam//g')
# NORMAL_ID=$(basename $NORMAL_BAM | sed 's/.bam//g')
full_id=${TUMOR_ID}.vs.${NORMAL_ID}

#----------------
# Mutect2 variant calling
# PoN must be calculated before hand
#----------------

# 1. Step1. Simple variant calling using PoN (previously calculated)
if [ "$GNOMAD" = "Unknown" ]; then
    run1="$GATK Mutect2 -R $REFGENOME \
     -L $BED \
     -I $TUMOR_BAM \
     -I $NORMAL_BAM \
     -tumor $TUMOR_ID \
     -normal $NORMAL_ID \
     -pon $PON \
     --f1r2-tar-gz $out/f1r2.tar.gz \
     -O $out/$full_id.unfiltered.vcf.gz"
else
    run1="$GATK Mutect2 -R $REFGENOME \
     -L $BED \
     -I $TUMOR_BAM \
     -I $NORMAL_BAM \
     -tumor $TUMOR_ID \
     -normal $NORMAL_ID \
     -germline-resource $GNOMAD \
     -pon $PON \
     --f1r2-tar-gz $out/f1r2.tar.gz \
     -O $out/$full_id.unfiltered.vcf.gz"
fi

# 2. LearnReadOrientationModelrun
run2="$GATK LearnReadOrientationModel -I $out/f1r2.tar.gz \
 -O $out/read-orientation-model.tar.gz"

# 3. Run GetPileupSummaries to summarize read support for a set number of known variant sites
# 4. Estimate contamination with CalculateContamination.
# In mouse we cannot to run this step, as there is no population study providing allele frequencies for common variants

# 5.Filter somatic calls
run3="$GATK FilterMutectCalls -V $out/$full_id.unfiltered.vcf.gz \
 -R $REFGENOME \
 -L $BED \
 --ob-priors $out/read-orientation-model.tar.gz \
 -O $out/$full_id.filtered.vcf.gz"

# 6. Filter for PASS calls
run4="zgrep '#\|PASS' $out/$full_id.filtered.vcf.gz > $out/$full_id.filtered.pass.vcf"

## Submit command to cluster
# Cluster commands
rm -f $logs/Mutect2.$full_id.o.log $logs/Mutect2.$full_id.e.log
echo bsub -M 20G -n 4 -o $logs/Mutect2.$full_id.o.log -e $logs/Mutect2.$full_id.e.log  -J Mutect2.$full_id "\"$run1 && $run2 && $run3 && $run4\"" > $out/cmd.logs
cat $out/cmd.logs | sh



