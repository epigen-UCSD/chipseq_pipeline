#/bin/bash

# input
samples=$1
PBS_ARRAYID=$2
FASTQDIR=$3
WORKDIR=$4
PBS_NP=1

source activate aquas_chipseq


# select libs
samplenames=(`cat $samples`)
INPREFIX=${samplenames[${PBS_ARRAYID}*3]} #index start from 0
genome=${samplenames[${PBS_ARRAYID}*3+1]}
ctl=${samplenames[${PBS_ARRAYID}*3+2]} #control or treatment
fastq1="${FASTQDIR}/${INPREFIX}_R1.fastq.gz"
fastq2="${FASTQDIR}/${INPREFIX}_R2.fastq.gz"
[[ ! -f $fastq1 ]] && {echo "fastq1 not found"; exit 0;}
[[ ! -f $fastq2 ]] && {echo "fastq2 not found"; exit 0;}


# parameters
type="histone"
true_rep="true"
no_pseudo_rep="true"
outdir="${WORKDIR}${INPREFIX}"; mkdir -p $outdir

############################################################
# Step 1
############################################################
fastq_input="-fastq1_1 ${fastq1} -fastq1_2 ${fastq2}"
[[ "$ctl" = true ]] && fastq_input="-ctl_fastq1_1 ${fastq1} -ctl_fastq1_2 ${fastq2}"

#final_stage="filt_bam"
final_stage="xcor"

# run pipeline
bds /projects/ps-epigen/software/chipseq_pipeline/chipseq.bds \
    -type $type \
    -final_stage $final_stage \
    -true_rep $true_rep \
    -no_pseudo_rep $no_pseudo_rep \
    -pe_no_trim_fastq \
    -nth $PBS_NP \
    -out_dir $outdir \
    -species $genome \
    $fastq_input

wait

############################################################
# Step 2
############################################################

if [ "$ctl" = true ]
then
    runFastQC_screen.sh $INPREFIX
    results_transfer.sh $INPREFIX
    exit 0 # no step 2 for ctl
fi

[[ ! -z ${no_peak+x} ]] &&  exit 0 # no step 2 if -v no_peak

ctl_id=$(cat $samples | grep true -n |cut -d':' -f1)
ctl_id=$[$ctl_id-1] #0-based index
ctl_name=${samplenames[${ctl_id}*3]}
input_dir="${WORKDIR}${ctl_name}/align";
ctl_prefix="$input_dir/ctl1/${samplenames[${ctl_id}*3]}" #index start from 0
tr_prefix="$outdir/align/rep1/${samplenames[${PBS_ARRAYID}*3]}"
tag_suffix="_R1.PE2SE.nodup.bam"
tag_input="-filt_bam ${tr_prefix}${tag_suffix} -ctl_filt_bam ${ctl_prefix}${tag_suffix}"

# run pipeline
bds /projects/ps-epigen/software/chipseq_pipeline/chipseq.bds \
    -type $type \
    -true_rep $true_rep \
    -no_pseudo_rep $no_pseudo_rep \
    -nth $PBS_NP \
    -out_dir $outdir \
    -species $genome \
    -pe \
    $tag_input

wait

############################################################
# runFastQC & fastq_screen
############################################################

#runFastQC_screen.sh $INPREFIX
#results_transfer.sh $INPREFIX
