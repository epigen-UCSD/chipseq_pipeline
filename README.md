Transcription Factor and Histone ChIP-Seq processing pipeline
==============
* This pipeline is customized/developed on top of [kundajelab/chipseq_pipeline@dec296](https://github.com/kundajelab/chipseq_pipeline/commit/dec2962f51951269dc2d2517f343ac991ceb36cc)
* Legend [README](./README_legend.md)

# Changes logs 
1. added [chipqc](./chipqc) module 
2. enabled xcor for ctrl 
3. enabled pileup track for both ctrl and treatments (ongoing) 

# Usage 

``` shell
$bds ~/data/software/chipseq_pipeline/chipseq.bds

== chipseq pipeline settings
	-type <string>                        : Type of ChIP-Seq pipeline. TF or histone (default: TF).
	-histone <bool>                       : (LEGACY PARAM) Histone ChIP-Seq. Equivalent to '-type histone'
	-final_stage <string>                 : Final stage for pipeline (bam, filt_bam, tag, xcor and peak).
	-sig_trk_for_pooled_rep_only <bool>   : Generate MACS2 signal tracks for pooled replicate only.
	-aligner <string>                     : Aligner to map raw reads in FASTQs (default: bwa).
	-subsample_xcor <string>              : # reads to be subsampled for cross corr. analysis (default: 15M).
	-subsample_chip <string>              : # reads to subsample exp. replicate. Subsampled tagalign will be used for steps downstream (default: 0; no subsampling).
	-subsample <string>                   : (LEGACY PARAM) # reads to subsample exp. replicate. Subsampled tagalign will be used for steps downstream (default: 0; no subsampling).
	-subsample_ctl <string>               : # reads to subsample control if non-zero (recommended: 40M or lower).
	-anon_filt_bam <bool>                 : Generate an annomymized filtered bam. This will not affect tasks downsteam.
	-peak_caller <string>                 : Peak caller for IDR analysis (default: spp for TF ChIP-seq and macs2 for Histone ChIP-seq ).
	-ctl_depth_ratio <real>               : Cut-off ratio of two control tagaligns for pooling (default: 1.2).
	-fraglen <int>                        : (LEGACY PARAM) Manually specify fragment length (default: 0, 0 means parsing fraglen from cross-corr analysis log file).
	-idr_rank <string>                    : Scoring column in narrow peak files for IDR. If not specified, signal.value for SPP peaks (TF) and p.value for MACS2 peaks (histone) are used.
	-idr_thresh <real>                    : IDR threshold : -log_10(score) (default: 0.05).
	-use_pooled_ctl <bool>                : Force to use pooled control (ignoring criteria to choose control for each replicate).
	-true_rep <bool>                      : Call peaks on true replicates only.
	-no_pseudo_rep <bool>                 : Do not call peaks on self pseudo replicates.
	-no_xcor <bool>                       : Disable cross-correlation analysis.
	-no_gpeak_filt <bool>                 : Disable gapped peak filtering through narrow peak (for histone ChIP-Seq only).
	-no_browser_tracks <bool>             : Disable generation of genome browser tracks (workaround for bzip2 shared library issue).
	-pe_xcor_only <bool>                  : (PE ONLY) Align R1 of paired end fastqs only and for cross-correlation analysis. All other analyses and QCs will be disabled.
	-pe_no_trim_fastq <bool>              : (PE ONLY) No fastq trimming and use PE tagAlign for cross-correlation analysis.
	-no_jsd <bool>                        : Disable JSD plot generation
	-no_chipqc <bool>                     : Disable advanced chipseq qc
== configuration file settings
	-c <string>                           : Configuration file path.
	-env <string>                         : Environment file path.
== parallelization settings
	-no_par <bool>                        : Serialize all tasks (individual tasks can still use multiple threads up to '-nth').
	-nth <int>                            : Maximum # threads for a pipeline. (default: 8).
== cluster/system/resource settings
	-wt <string>                          : Walltime for all single-threaded tasks (example: 8:10:00, 3h, 3600, default: 5h50m, 5:50:00).
	-memory <string>                      : Maximum memory for all single-threaded tasks (equivalent to '-mem', example: 4.5G, 1024M, default: 7G).
	-use_system <string>                  : Force to use a system (equivalent to 'bds -s [SYSTEM_NAME] ...', any system defined in bds.config can be used).
	-nice <int>                           : Set process priority for all tasks (default: 0; -20 (highest) ~ 19 (lowest) ).
	-retrial <int>                        : # of Retrial for failed tasks (default: 0).
	-q <string>                           : Submit tasks to a specified cluster queue.
	-unlimited_mem_wt <bool>              : Use unlimited max. memory and walltime.
	-java_tmp_dir <string>                : Java temporary directory. (change it when you get 'Disk quota exceeded' error in Java, default: ${TMPDIR}).
== shell environment settings
	-mod <string>                         : Modules separated by ; (example: "bowtie/2.2.4; bwa/0.7.7; picard-tools/1.92").
	-shcmd <string>                       : Shell commands separated by ;. Shell var. must be written as ${VAR} not as $VAR (example: "export PATH=${PATH}:/usr/test; VAR=test").
	-addpath <string>                     : Path separated by ; or : to be PREPENDED to \$PATH (example: "/bin/test:${HOME}/utils").
	-conda_env <string>                   : Anaconda Python (or Miniconda) environment name for all softwares including Python2.
	-conda_env_py3 <string>               : Anaconda Python (or Miniconda) environment name for Python3.
	-conda_bin_dir <string>               : Anaconda Python (or Miniconda) bin directory.
	-cluster_task_min_len <int>           : Minimum length for a cluster job in seconds (dealing with NFS delayed write, default: 60).
	-cluster_task_delay <int>             : Constant delay for every job in seconds (dealing with NFS delayed write, default: 0).
== output/title settings
	-out_dir <string>                     : Output directory (default: out).
	-title <string>                       : Prefix for HTML report and outputs without given prefix.
== species settings
	-species <string>                     : Species. need to specify '-species_file' too if you have not installed genome database with 'install_genome_data.sh'.
	-species_file <string>                : Species file path.
	-species_browser <string>             : Species name in WashU genome browser.
	-ref_fa <string>                      : Reference genome sequence fasta.
	-chrsz <string>                       : Chromosome sizes file path (use fetchChromSizes from UCSC tools).
	-blacklist <string>                   : Blacklist bed.
	-seq_dir <string>                     : Reference genome sequence directory path (where chr*.fa exist).
== ENCODE accession settings
	-ENCODE_accession <string>            : ENCODE experiment accession ID (or dataset).
	-ENCODE_award_rfa <string>            : ENCODE award RFA (e.g. ENCODE3).
	-ENCODE_assay_category <string>       : ENCODE assay category.
	-ENCODE_assay_title <string>          : ENCODE assay title.
	-ENCODE_award <string>                : ENCODE award (e.g. /awards/U41HG007000/).
	-ENCODE_lab <string>                  : Lab (e.g. /labs/anshul-kundaje/)
	-ENCODE_assembly <string>             : hg19, GRCh38, mm9, mm10.
	-ENCODE_alias_prefix <string>         : Alias prefix, Alias = alias_prefix: + filename + alias_suffix
	-ENCODE_alias_suffix <string>         : Alias suffix, Alias = alias_prefix: + filename + alias_suffix
== report settings
	-url_base <string>                    : URL base for output directory.
	-viz_genome_coord <string>            : WashU genome browser genome coordinate (e.g. chr7:27117661-27153380).
== fastq input definition :
        Single-ended : For replicate '-fastq[REP_ID]', For control '-ctl_fastq[REP_ID]'
        Paired end : For replicate '-fastq[REP_ID]_[PAIR_ID]', For control '-ctl_fastq[REP_ID]_[PAIR_ID]'
== bam input (raw or filtered) definition :
        Raw bam : For replicate '-bam[REP_ID]', For control '-ctl_bam[REP_ID]'.
        Filtered bam : For replicate '-filt_bam[REP_ID]', For control '-ctl_filt_bam[REP_ID]'.
== tagalign input definition :
        For replicate '-tag[REP_ID]', For control '-ctl_tag[REP_ID]'.
== narrow peak input definition : 
        For true replicates, use '-peak1' and '-peak2',
        For pooled replicates, use '-peak_pooled',
        For two PR (self-pseudo-replicates), use '-peak[REP_ID]_pr1' and '-peak[REP_ID]_pr2'
        For two PPR (pooled pseudo-replicates), use '-peak_ppr1' and '-peak_ppr2'
== input endedness settings (SE or PE) :
	-se <bool>                            : Singled-ended data set. To specify it for each replicate, '-se[REP_ID]' for exp. reps, '-ctl_se[CTL_ID]' for control.
	-pe <bool>                            : Paired end data set. To specify it for each replicate, '-pe[REP_ID]' for exp. reps, '-ctl_pe[CTL_ID]' for controls.
== align bwa settings (requirements: -bwa_idx)
	-param_bwa_aln <string>               : Parameters for bwa aln (default: "-q 5 -l 32 -k 2").
	-bwa_idx <string>                     : BWA index (full path prefix of *.bwt file) .
	-wt_bwa <string>                      : Walltime for bwa (default: 47, 47:00:00).
	-mem_bwa <string>                     : Max. memory for bwa (default: 12G).
== fastq trimmer settings
	-trim_bp <int>                        : Number of basepairs after trimming fastqs (default: 50).
== align multimapping settings
	-multimapping <int>                   : # alignments reported for multimapping (default: 0).
== postalign bam settings
	-mapq_thresh <int>                    : Threshold for low MAPQ reads removal (default: 30).
	-rm_chr_from_tag <string>             : Perl style reg-ex to exclude reads from tag-aligns. (example: 'other|ribo|mito|_', '_', default: blank)
	-no_dup_removal <bool>                : No dupe removal when filtering raw bam.
	-wt_dedup <string>                    : Walltime for post-alignment filtering (default: 23h, 24:00:00).
	-mem_dedup <string>                   : Max. memory for post-alignment filtering (default: 12G).
	-dup_marker <string>                  : Dup marker for filtering mapped reaads in BAMs: picard or sambamba (default: picard).
	-use_sambamba_markdup <bool>          : Use sambamba markdup instead of Picard MarkDuplicates (default: false).
== postalign bed/tagalign settings
	-mem_shuf <string>                    : Max. memory for UNIX shuf (default: 12G).
	-fraglen0 <bool>                      : (LEGACY PARAM) Set predefined fragment length as zero for cross corr. analysis (add -speak=0 to run_spp.R).
	-speak_xcor <int>                     : Set user-defined cross-corr. peak strandshift (-speak= in run_spp.R). Use -1 to disable (default: -1).
	-extra_param_xcor <string>            : Set extra parameters for run_spp.R (cross-corr. analysis only).
	-mem_xcor <string>                    : Max. memory for cross-corr. analysis (default: 15G).
== postalign bed/tagalign settings
	-mem_shuf <string>                    : Max. memory for UNIX shuf (default: 12G).
	-fraglen0 <bool>                      : (LEGACY PARAM) Set predefined fragment length as zero for cross corr. analysis (add -speak=0 to run_spp.R).
	-speak_xcor <int>                     : Set user-defined cross-corr. peak strandshift (-speak= in run_spp.R). Use -1 to disable (default: -1).
	-extra_param_xcor <string>            : Set extra parameters for run_spp.R (cross-corr. analysis only).
	-mem_xcor <string>                    : Max. memory for cross-corr. analysis (default: 15G).
== callpeak spp settings
	-cap_num_peak_spp <string>            : Cap number of peaks (-npeak= in run_spp.R) (default: 300000).
	-max_ppsize_spp <string>              : R stack size (R parameter --max-ppsize=; between 5000 and 5000000) for SPP.
	-speak_spp <int>                      : User-defined cross-corr. peak strandshift (-speak= in run_spp.R). Use -1 to get from upstream cross-corr. analysis (default: -1).
	-extra_param_spp <string>             : Extra parameters for SPP (run_spp.R, peak calling only).
	-wt_spp <string>                      : Walltime for spp (default: 47h, 47:00:00).
	-mem_spp <string>                     : Max. memory for spp (default: 12G).
== callpeak gem settings
	-npeak_gem <int>                      : Threshold on # of peaks for GEM (default: 300000).
	-k_min_gem <int>                      : Minimum length of k-mers (--k_min in GEM, default: 6).
	-k_max_gem <int>                      : Maximum length of k-mers (--k_max in GEM, default: 13).
	-q_val_thresh_gem <real>              : Q-value threshold (--q in GEM, default: 0).
	-read_dist_gem <string>               : Read distribution txt file for GEM (default: $script_dir/etc/Read_Distribution_default.txt).
	-extra_param_gem <string>             : Extra parameters for GEM.
	-wt_gem <string>                      : Walltime for GEM (default: 47h, 47:00:00).
	-mem_gem <string>                     : Max. memory for GEM (default: 15G).
== callpeak PeakSeq settings
	-target_fdr_peakseq <real>            : Target FDR for PeakSeq (default: 0.05).
	-n_sim_peakseq <int>                  : Number of simulations for PeakSeq (default: 10).
	-enrich_mapped_fraglen_peakseq <int>  : Enrichment mapped fragment length for PeakSeq. Use -1 to get from upstream cross-corr. analysis (default: -1).
	-min_interpeak_dist_peakseq <int>     : Minimum interpeak distance for PeakSeq. Use -1 to get from upstream cross-corr. analysis (default: -1).
	-mappability_map_peakseq <string>     : Mappability map file for PeakSeq (http://archive.gersteinlab.org/proj/PeakSeq/Mappability_Map).
	-max_qval_peakseq <real>              : Maximum Q-value for PeakSeq (default: 0.1).
	-bckgrnd_model_peakseq <string>       : Background model for PeakSeq (default: Simulated).
	-extra_param_peakseq <string>         : Extra parameters for PeakSeq.
	-wt_peakseq <string>                  : Walltime for PeakSeq (default: 47h, 47:00:00).
	-mem_peakseq <string>                 : Max. memory for PeakSeq (default: 12G).
== callpeak macs2 settings (requirements: -chrsz -gensz)
	-gensz <string>                       : Genome size; hs for human, mm for mouse.
	-wt_macs2 <string>                    : Walltime for MACS2 (default: 23h, 23:00:00).
	-mem_macs2 <string>                   : Max. memory for MACS2 (default: 15G).
	-pval_thresh_macs2 <real>             : --pvalue for macs2 callpeak (https://github.com/taoliu/MACS#-p--pvalue, default: 0.01)
	-keep_dup_macs2 <string>              : --keep-dup for macs2 callpeak (https://github.com/taoliu/MACS#--keep-dup, default: all).
	-extsize_macs2 <int>                  : --extsize for macs2 callpeak (https://github.com/taoliu/MACS#--extsize). Use -1 to get from upstream cross-corr. analysis (default: -1).
	-shift_macs2 <int>                    : --shift for macs2 callpeak (https://github.com/taoliu/MACS#--shift, default: 0).
	-cap_num_peak_macs2 <string>          : Cap number of peaks by taking top N peaks for MACS2 (default: 500K).
	-extra_param_macs2 <string>           : Extra parameters for macs2 callpeak.
== callpeak naive overlap settings
	-nonamecheck <bool>                   : bedtools intersect -nonamecheck (bedtools>=2.24.0, use this if you get bedtools intersect naming convenction warnings/errors).
== IDR settings
	-idr_suffix <bool>                    : Append IDR threshold to IDR output directory.
== CHIPQC settings
	-tss_enrich <string>                  : TSS enrichment bed for chipqc.
	-prom <string>                        : Promoter bed (promoter region file) for chipqc.
	-enh <string>                         : Enhancer bed (enhancer region file) for chipqc.
	-mem_chipqc <string>                  : Max. memory for CHIPQC (default: 20G).
	-wt_chipqc <string>                   : Walltime for CHIPQC (default: 47h, 47:00:00).

```

