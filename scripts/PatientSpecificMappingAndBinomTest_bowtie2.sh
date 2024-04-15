#!/bin/bash
#PatientSpecificMappingAndBinomTest_bowtie2.sh

################################################################################
# The script below was used to map fastqs to the hg38 genome for comparison to
#     a mapping to the personalized vg genome.
#
#     The script below is tailored to our local compute cluster's specific
#     architecture, and will need to be rewritten for your own compute environment's
#     specific needs if you would like to reproduce its results.  In particular,
#     the various sourcings will likely need to be removed, as well as the log_msg
#     commands.  One will need to provide or define their own TEMP_DIR.
#
# Expected to be in the user's PATH:
#
#     bcftools version 1.9
#     samtools
#     bowtie2 version 2.2.5
#     cutadapt
#     MarkDuplicates.jar from picard-tools version 1.88
#     bedtools
#
# Required Environment Variables:
#     LIBRARY             - Name of library to be processed
#           GENOME              - Name of the genome to use through processing
#     REFERENCE_DIR       - Location of the patient-specific genome
#     FASTQ_FILE         - The FASTQ of interest
#     OUTPUT_DIR          - Change default output directory
#     THE_TARGET            - A target against which to perform the bcftools pileup.
#                           generated via a bcftools query and tabix.
#
# Optional Environment Variables:
#     CONTROL             - Flags that the library is a control library
#                           will skip running if a down-sampled version
#                           of the library exists
#     ADAPTER_SEQUENCE    - Alternative adapter sequence for cutadapt
#
################################################################################



export SOURCE_DIR="/gpfs/gpfs2/sdi/etc"
export RPATH="/gpfs/gpfs2/software/R-3.4.3"
export PIPELINE_VENV_PYTHON3="/gpfs/gpfs2/sdi/software/pipeline-venv-python3/bin/activate"

### Mandatory sourcing of bashrc for necessary environment variables. ###
if [ -e $SOURCE_DIR/bashrc ]; then
    . $SOURCE_DIR/bashrc
else echo "[fatal] - Could not find myerslab bashrc file. Exiting"; exit 1; fi

### Mandatory sourcing of functions to get helper functions (like call_cmd). ###
if [ -e $SOURCE_DIR/functions ]; then
    . $SOURCE_DIR/functions
else echo "[fatal] - Could not find functions file. Exiting"; exit 1; fi

if [ -z "$ADAPTER_SEQUENCE" ]; then
    ADAPTER_SEQUENCE="/gpfs/gpfs1/myerslab/reference/sequence-adapters/adapters.fa"
fi

### Verify we are not running on the head node. ###
if [ -z "$LSB_JOBID" ]; then echo "Please run on a compute node. Exiting"; exit 1; fi
if [ -z "$LIBRARY" ];            then empty_param_quit "LIBRARY"; fi
if [ -z "$GENOME" ];             then empty_param_quit "GENOME"; fi


### Set up temp directory variable ###
if [ -z "$TEMP_DIR" ]; then TEMP_DIR=$(get_temp_dir); fi
if [ ! -d "$TEMP_DIR" ]; then
    if ! mkdir -p "$TEMP_DIR"; then
      echo "Could not create output dir: $TEMP_DIR. Exiting";
      exit 1;
    fi
fi

### Activate Python3 vevn ###
if [ -e $PIPELINE_VENV_PYTHON3 ]; then
    . $PIPELINE_VENV_PYTHON3
else echo "[fatal] - Could not find python3 activate file[$PIPELINE_VENV_PYTHON3]. Exiting"; exit 1; fi

export PATH=$PATH:/gpfs/gpfs2/software/samtools-1.8/bin
export LOGFILE_NAME=$OUTPUT_DIR/$LIBRARY.log



NTHREADS=8

log_msg info "Beginning Fastq to TagAlign"
log_msg info "    LIBRARY:              $LIBRARY"
log_msg info "    GENOME:               $GENOME"
log_msg info "    OUTPUT_DIR:           $OUTPUT_DIR"

log_msg info "Consolidating fastq files..."
# Set up output files
FASTQ_FILE_1=$TEMP_DIR/$LIBRARY.fastq.gz
SAI_FILE_1=$TEMP_DIR/$LIBRARY.sai
RAW_BAM_PREFIX=$LIBRARY.raw.srt
RAW_BAM_FILE=$OUTPUT_DIR/$RAW_BAM_PREFIX.bam
RAW_BAM_FILE_MAPSTATS=$OUTPUT_DIR/$RAW_BAM_PREFIX.flagstat.qc

# Run cutadapt on each of the fastqs
for fastq in $FASTQ_FILE ; do
  base_fastq_name="$(basename $fastq .fastq.gz)"
  flowcell_id=$(echo "$base_fastq_name" | cut -d_ -f1)

  CUTADAPT_INPUT="$fastq"
  CUTADAPT_OUTPUT="$TEMP_DIR/${base_fastq_name}_cutadapt.fastq.gz"
  CUTADAPT_STATISTICS="$OUTPUT_DIR/${base_fastq_name}_cutadapt_report.out"
  CMD="cutadapt -a file:$ADAPTER_SEQUENCE -j 8 -m 40 -o $CUTADAPT_OUTPUT $CUTADAPT_INPUT > $CUTADAPT_STATISTICS"
  run_cmd "$CMD" "$CUTADAPT_OUTPUT"
done

TRIMMED_FASTQS="$(ls ${TEMP_DIR}/*cutadapt.fastq.gz)"

CMD="cat $TRIMMED_FASTQS > $FASTQ_FILE_1"
run_cmd "$CMD" "$FASTQ_FILE_1"

#BWA_SOFTWARE=$(get_software_dir bwa-0.7.12/bwa)
BOWTIE2_SOFTWARE=$(get_software_dir bowtie2-2.2.5/bowtie2)

#ALN_SAM=$OUTPUT_DIR/${LIBRARY}.raw.srt.sam
ALN_SAM=$OUTPUT_DIR/${LIBRARY}.raw.sam

#log_msg info "Aligning with bwa aln..."
log_msg info "Aligning with bowtie2..."
#CMD="$BWA_SOFTWARE aln -q 5 -l 20 -k 3 -t $NTHREADS $REFERENCE_DIR/$GENOME.fa $FASTQ_FILE_1 > $SAI_FILE_1"
CMD="$BOWTIE2_SOFTWARE --very-sensitive -x $REFERENCE_DIR/$GENOME.fa -U $FASTQ_FILE_1 > $ALN_SAM"
run_cmd "$CMD" "$ALN_SAM"



#Sort the sam file using samtools sort, done by chromosomal location.
SRT_SAM=$OUTPUT_DIR/${LIBRARY}.raw.srt.sam
#log_msg info "sorting sam file by read name..."
log_msg info "sorting sam file by mapping location..."
CMD="samtools sort --threads $NTHREADS -O SAM -o $SRT_SAM -T $OUTPUT_DIR/$RAW_BAM_PREFIX $ALN_SAM"
run_cmd "$CMD" "$SRT_SAM"

rm $ALN_SAM

#Getting some mapping statistics from samtools
CMD="samtools flagstat $SRT_SAM > $RAW_BAM_FILE_MAPSTATS"
run_cmd "$CMD" "$RAW_BAM_FILE_MAPSTATS"


#Here we filter by quality score (def:30) in preparation for removing PCR duplicates.
#We simultaneously convert to bam format.
log_msg info "Converting raw.srt.sam to filt.srt.bam"
# Set up filter bam output files
FILT_BAM_PREFIX=$LIBRARY.filt.srt
FILT_BAM_FILE=$OUTPUT_DIR/$FILT_BAM_PREFIX.bam
DUP_FILE_QC=$OUTPUT_DIR/$FILT_BAM_PREFIX.dup.qc
MAPQ_THRESH=30

CMD="samtools view -F 1548 -q $MAPQ_THRESH -b $SRT_SAM > $FILT_BAM_FILE"
run_cmd "$CMD" "$FILT_BAM_FILE"




######################
#This marks PCR duplicates, based on position and orientation of reads.
#If two reads share the same start position and the same orientation, the read
#with the lower FASTQ quality score is removed.
######################
#BRYAN NOTE:  I changed the parameter -Xmx4G to -Xmx16G to avoid errors due to
#insufficient memory.
TMP_FILT_BAM_FILE=$TEMP_DIR/$FILT_BAM_PREFIX.dupmark.bam
PICARD_JAR=$(get_software_dir picard-tools-1.88/MarkDuplicates.jar)
CMD="java -Xmx16G -jar $PICARD_JAR INPUT=$FILT_BAM_FILE OUTPUT=$TMP_FILT_BAM_FILE METRICS_FILE=$DUP_FILE_QC VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false"
run_cmd "$CMD" "$TMP_FILT_BAM_FILE $DUP_FILE_QC"
mv "$TMP_FILT_BAM_FILE" "$FILT_BAM_FILE"



#This section filters the sam file to remove unmapped reads and PCR dups, indexes it,
#and extracts quality metrics from the sam.
# Set up final bam output files
FINAL_BAM_PREFIX=$LIBRARY.filt.nodup.srt
FINAL_BAM_FILE=$OUTPUT_DIR/$FINAL_BAM_PREFIX.bam
#FINAL_BAM_INDEX_FILE=$OUTPUT_DIR/$FINAL_BAM_PREFIX.bai
FINAL_BAM_FILE_MAPSTATS=$OUTPUT_DIR/$FINAL_BAM_PREFIX.flagstat.qc

CMD="samtools view -F 1548 -b $FILT_BAM_FILE > $FINAL_BAM_FILE"
run_cmd "$CMD" "$FINAL_BAM_FILE"
#CMD="samtools index $FINAL_BAM_FILE $FINAL_BAM_INDEX_FILE"
#run_cmd "$CMD" "$FINAL_BAM_INDEX_FILE"
CMD="samtools index $FINAL_BAM_FILE"
run_cmd "$CMD" "$FINAL_BAM_FILE.bai"
CMD="samtools flagstat $FINAL_BAM_FILE > $FINAL_BAM_FILE_MAPSTATS"
run_cmd "$CMD" "$FINAL_BAM_FILE_MAPSTATS"


#Get stats on library complexity...
# Set up library complexity output file
PBC_FILE_QC=$OUTPUT_DIR/$FINAL_BAM_PREFIX.pbc.qc

#This converts the bam file to a bed file, remove the reads mapping to the mitochondrial genome, sort the bed file, count the number of unique entries,
#then run an awk command that I don't understand.  However, the output is a quality control file, specifically for the BAM file.
BEDTOOLS=$(get_software_dir bedtools2-2.20.0/bin/bedtools)
CMD="$BEDTOOLS bamtobed -i $FILT_BAM_FILE | awk 'BEGIN{OFS=\"\t\"}{print \$1,\$2,\$3,\$6}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} (\$1==1){m1=m1+1} (\$1==2){m2=m2+1} {m0=m0+1} {mt=mt+1} END{printf \"%d\t%d\t%d\t%d\t%f\t%f\t%f\n\",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > $PBC_FILE_QC"
run_cmd "$CMD" "$PBC_FILE_QC"


#Convert final BAM to a SAM and sort by read name
FINAL_SAM_FILE=$OUTPUT_DIR/$FINAL_BAM_PREFIX.sam
CMD="samtools view -h $FINAL_BAM_FILE | samtools sort -n --threads $NTHREADS -O SAM -o $FINAL_SAM_FILE -T $OUTPUT_DIR/$RAW_BAM_PREFIX"
run_cmd "$CMD" "$FINAL_SAM_FILE"


module load g/bcftools/1.9

OUTPUT_VCF=$OUTPUT_DIR/${FINAL_BAM_PREFIX}.pileup.vcf

bcftools mpileup -T ${THE_TARGET} -d 100000 -f ${REFERENCE_DIR}/${GENOME}.fa ${FINAL_BAM_FILE} | bcftools call -T ${THE_TARGET} -m -Ov -o ${OUTPUT_VCF}





deactivate
