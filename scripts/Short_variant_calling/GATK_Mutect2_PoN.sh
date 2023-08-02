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
parser.add_argument('-normals','--normals', required=True, help='List of normal bam files')
parser.add_argument('-bed','--bed', required=True, help='WES bed file')
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
echo  Reference genome: "$REFGENOME"
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

#----------------------
# Run PoN step 1
# Vcf calculation
#----------------------

out2=$out/vcfs
mkdir -p $out2

for normal_bam in $(cat $NORMALS); do
    
    #-------------------------------------------
    #   Get parameters and path to files
    #-------------------------------------------
    
    # Normal id
    sample=$(basename $normal_bam | sed 's/.bam$//g')
    
    run1="$GATK Mutect2 -R $REFGENOME -I $normal_bam --max-mnp-distance 0 -O $out2/${sample}.vcf.gz -L $BED"
    
    touch $out2/${sample}.vcf.gz

    rm -f ${logs}/PoN1.${sample}.o.log ${logs}/PoN1.${sample}.e.log
    echo bsub -M 15G -n 1 -o ${logs}/PoN1.${sample}.o.log -e ${logs}/PoN1.${sample}.e.log -J PoN1.${sample} "\"$run1\"" >> $out/cmd.logs

    # Collapse jobs to hold on for next step of bsubs
    hold_jobs="$hold_jobs && ended(PoN1.${sample})"

done


#----------------------
# Run PoN step 2
#----------------------
# Wait bsub command
wait_cmd="'$hold_jobs'"

# Intervals file
out3=$out/interval-files-folder
mkdir -p $out3

run2="$GATK SplitIntervals \
   -R $REFGENOME \
   -L $BED \
   --scatter-count 1 \
   -O $out3"

# Path to normal vcfs plus "-V" in front for GATK run 
samples=''
for normal in $(ls $out2/*gz);do
    samples="$samples -V $normal"
done


# Main command
rm -rf $out/pon_db
run3="$GATK GenomicsDBImport -R $REFGENOME -L $BED --merge-input-intervals true --genomicsdb-workspace-path $out/pon_db"

# Add samples
run3="$run3 $samples"

# Run in the cluster
echo bsub -w $wait_cmd -M 10G -n 1 -o ${logs}/PoN2.o.log -e ${logs}/PoN2.e.log -J PoN2.${PREFIX} "\"$run2 & $run3\"" >> $out/cmd.logs
hold_jobs="$hold_jobs && ended(PoN2.${PREFIX})"


#----------------------
# Run PoN step 3
#----------------------
# Wait bsub command
wait_cmd="'$hold_jobs'"

# Final out folder
out4=$out/final_PoN
mkdir -p $out4

# Main command
run4="$GATK CreateSomaticPanelOfNormals -R $REFGENOME -V gendb://pon_db -O $out4/${PREFIX}.pon.vcf.gz"

# Run in the cluster
echo bsub -w $wait_cmd -M 10G -n 1 -o ${logs}/PoN3.o.log -e ${logs}/PoN3.e.log -J PoN3.${PREFIX} "\"cd $out && $run4\"" >> $out/cmd.logs

sed -i "s/' && /'/g" $out/cmd.logs

# Submit jobs
cat $out/cmd.logs | sh


