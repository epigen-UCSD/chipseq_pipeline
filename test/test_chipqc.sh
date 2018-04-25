#!/bin/bash
source activate aquas_chipseq
samples=$1
PBS_ARRAYID=$2
FASTQDIR=$3
WORKDIR=$4
PBS_NP=1

samplenames=(`cat $samples`)

prefix_basename=${samplenames[${PBS_ARRAYID}*3]} #index start from 0
o_dir="${WORKDIR}${prefix_basename}_chip"
w_dir="${WORKDIR}${prefix_basename}_chip/qc/ctl1"
species_chipqc=${samplenames[${PBS_ARRAYID}*3+1]}
chrsz=/home/zhc268/data/GENOME/${species_chipqc}/${species_chipqc}.chrom.sizes
ref_fa=/home/zhc268/data/GENOME/${species_chipqc}/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
tss_enrich=/projects/ps-epigen/GENOME/hg38/ataqc/hg38_gencode_tss_unique.bed.gz
pbc_log=${o_dir}/qc/ctl1/${prefix_basename}_R1.PE2SE.nodup.pbc.qc
fq1="${FASTQDIR}/${prefix_basename}_R1.fastq.gz"
fq2="${FASTQDIR}/${prefix_basename}_R2.fastq.gz"
param_fastq="--fastq1 ${fq1} --fastq2 ${fq2}"
bam=${o_dir}/align/ctl1/${prefix_basename}_R1.PE2SE.bam
dup_log=${o_dir}/qc/ctl1/${prefix_basename}_R1.PE2SE.dup.qc
filt_bam=${o_dir}/align/ctl1/${prefix_basename}_R1.PE2SE.nodup.bam
bed=${o_dir}/align/ctl1/${prefix_basename}_R1.PE2SE.nodup.tagAlign.gz
blacklist="/projects/ps-epigen/GENOME/hg38/hg38.blacklist.bed.gz"
prom=/projects/ps-epigen/GENOME/hg38/ataqc/reg2map_honeybadger2_dnase_prom_p2.hg19_to_hg38.bed.gz
enh=/projects/ps-epigen/GENOME/hg38/ataqc/reg2map_honeybadger2_dnase_enh_p2.hg19_to_hg38.bed.gz
param_blacklist="--blacklist $blacklist"
param_prom="--blacklist $prom"
param_enh="--blacklist $enh"


../chipqc/run_chipqc.py   --workdir $o_dir \
		    --outdir $w_dir \
		    --outprefix $prefix_basename \
		    --genome $species_chipqc \
		    --chromsizes $chrsz \
		    --ref $ref_fa \
		    --tss $tss_enrich \
		    --pbc $pbc_log\
		    $param_fastq \
		    --alignedbam $bam \
		    --coordsortbam $bam \
		    --duplog $dup_log \
		    --finalbam $filt_bam \
		    --finalbed $bed \
		    $param_blacklist $param_prom $param_enh



