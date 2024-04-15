#!/bin/bash
#DonorSpecificMapping.sh


################################################################################
#This script is used to map reads to the donor-specific graph genomes, and then
#     extract the reads which map to heterozygous locations for further counting
#     and significance assessment.
#
#     The script below is tailored to our local compute cluster's specific
#     architecture, and will need to be rewritten for your own compute environment's
#     specific needs if you would like to reproduce its results.  In particular,
#     the various sourcings will likely need to be removed, as well as the log_msg
#     commands.  One will need to provide or define their own TEMP_DIR. You will
#     also need to load packags to ensure that tools are available in your PATH,
#     and change the R script path to be called correctly.
#
#The script expects that the following are available:
# R version 4.1.0, but other versions are acceptable so long as they contain
#    the GenomicRanges package.
# python 3.6.6 is available.
# vg toolkit version 1.2.0 is available.
# piccard 2.21.6 is available.
#     Note that the MarkedDuplicates.jar must be in your PATH.
#
#This program takes the following arguments:
# LIBRARY: Name of the library to be processed.  Should be in the form of
#     TISSUE_FACTOR
# GENOME: Name of the genome to use throughout.
# REFERENCE_DIR: Directory containing the patient specific genome.  The genome
#     should be named identical to the GENOME argument.
# FASTQ_FILE: The fastq file to be mapped to the genome.
# ALL_HET_REGIONS: A FASTA file which contains all regions that contain a
#     non-reference base, whether heterozygous or homozygous. Produced by the
#     Building_fasta_variant_sequences.sh script.
#
#This program has the following optional arguments:
# ADAPTER_SEQUENCE: adapter sequence for cutadapt. Provided as AdapterSequence.fa
#
################################################################################

module load cluster/python/3.6.6

module list

NTHREADS=8

TEMP_DIR=${OUTPUT_DIR}/${LIBRARY}
mkdir ${TEMP_DIR}

# Set up output files
FASTQ_FILE_1=$OUTPUT_DIR/$LIBRARY.fastq.gz
FASTQ_FILE_2=$OUTPUT_DIR/$LIBRARY.fastq.gz
SAI_FILE_1=$OUTPUT_DIR/$LIBRARY.sai
RAW_BAM_PREFIX=$LIBRARY.raw.srt
RAW_BAM_FILE=$OUTPUT_DIR/$RAW_BAM_PREFIX.bam
RAW_BAM_FILE_MAPSTATS=$OUTPUT_DIR/$RAW_BAM_PREFIX.flagstat.qc


module load cluster/cutadapt/2.8
# Run cutadapt on each of the fastqs
for fastq in $FASTQ_FILE ; do
  base_fastq_name="$(basename $fastq .fastq.gz)"
  flowcell_id=$(echo "$base_fastq_name" | cut -d_ -f1)

  CUTADAPT_INPUT="$fastq"
  CUTADAPT_OUTPUT="$OUTPUT_DIR/${base_fastq_name}_cutadapt.fastq.gz"
  CUTADAPT_STATISTICS="$OUTPUT_DIR/${LIBRARY}_cutadapt_report.out"
  cutadapt -a file:$ADAPTER_SEQUENCE -j 8 -m 40 -o $CUTADAPT_OUTPUT $CUTADAPT_INPUT > $CUTADAPT_STATISTICS
done

TRIMMED_FASTQS="$(ls ${OUTPUT_DIR}/*cutadapt.fastq.gz)"

cat $TRIMMED_FASTQS > $FASTQ_FILE_1
cp $FASTQ_FILE_1 $FASTQ_FILE_2



module load cluster/picard/2.21.6
module load cluster/vg/1.20.0
module load cluster/samtools/1.8


ALN_GAM=$OUTPUT_DIR/${LIBRARY}.raw.gam
ALN_BAM=$OUTPUT_DIR/${LIBRARY}.raw.bam
ALN_SAM=$OUTPUT_DIR/${LIBRARY}.raw.sam
vg map -t ${NTHREADS} -A -K -M 3 -f ${FASTQ_FILE_2} -x ${REFERENCE_DIR}/${GENOME}.xg -g ${REFERENCE_DIR}/${GENOME}.gcsa -1 ${REFERENCE_DIR}/${GENOME}.gbwt > $ALN_GAM

vg surject -t $NTHREADS -b -x $REFERENCE_DIR/${GENOME}.xg -b $ALN_GAM > $ALN_BAM

#vg surject -t 4 -x /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/VariantGenome_Graph/hg38_with_variants.xg -s /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/MappedReads_vg/AMY_BCL11A_1224.raw.gam > /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/MappedReads_vg/AMY_BCL11A_1224.raw.sam


samtools view -h -o $ALN_SAM $ALN_BAM


#Sort the sam file using samtools sort, done by chromosomal location.
SRT_SAM=$OUTPUT_DIR/${LIBRARY}.raw.srt.sam
samtools sort --threads $NTHREADS -O SAM -o $SRT_SAM -T $OUTPUT_DIR/$RAW_BAM_PREFIX $ALN_SAM
rm $ALN_SAM



#Getting some mapping statistics from samtools
samtools flagstat $SRT_SAM > $RAW_BAM_FILE_MAPSTATS


#Here we filter by quality score (def:30) in preparation for removing PCR duplicates.
#We simultaneously convert to bam format.
# Set up filter bam output files
FILT_BAM_PREFIX=$LIBRARY.filt.srt
FILT_BAM_FILE=$OUTPUT_DIR/$FILT_BAM_PREFIX.bam
DUP_FILE_QC=$OUTPUT_DIR/$FILT_BAM_PREFIX.dup.qc
MAPQ_THRESH=30

samtools view -q $MAPQ_THRESH -b $SRT_SAM > $FILT_BAM_FILE




######################
#This marks PCR duplicates, based on position and orientation of reads.
#If two reads share the same start position and the same orientation, the read
#with the lower FASTQ quality score is removed.
######################
#BRYAN NOTE:  I changed the parameter -Xmx4G to -Xmx16G to avoid errors due to
#insufficient memory.
TMP_FILT_BAM_FILE=$TEMP_DIR/$FILT_BAM_PREFIX.dupmark.bam
PICARD_JAR=MarkDuplicates.jar
java -Xmx16G -jar $PICARD_JAR INPUT=$FILT_BAM_FILE OUTPUT=$TMP_FILT_BAM_FILE METRICS_FILE=$DUP_FILE_QC VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
mv "$TMP_FILT_BAM_FILE" "$FILT_BAM_FILE"




#This section filters the sam file to remove unmapped reads and PCR dups, indexes it,
#and extracts quality metrics from the sam.
# Set up final bam output files
FINAL_BAM_PREFIX=$LIBRARY.filt.nodup.srt
FINAL_BAM_FILE=$OUTPUT_DIR/$FINAL_BAM_PREFIX.bam
#FINAL_BAM_INDEX_FILE=$OUTPUT_DIR/$FINAL_BAM_PREFIX.bai
FINAL_BAM_FILE_MAPSTATS=$OUTPUT_DIR/$FINAL_BAM_PREFIX.flagstat.qc

samtools view -F 1548 -b $FILT_BAM_FILE > $FINAL_BAM_FILE
#CMD="samtools index $FINAL_BAM_FILE $FINAL_BAM_INDEX_FILE"
#run_cmd "$CMD" "$FINAL_BAM_INDEX_FILE"
samtools index $FINAL_BAM_FILE
samtools flagstat $FINAL_BAM_FILE > $FINAL_BAM_FILE_MAPSTATS


#Get stats on library complexity...
# Set up library complexity output file
PBC_FILE_QC=$OUTPUT_DIR/$FINAL_BAM_PREFIX.pbc.qc



#Convert final BAM to a SAM and sort by read name
FINAL_SAM_FILE=$OUTPUT_DIR/$FINAL_BAM_PREFIX.sam
samtools view -h $FINAL_BAM_FILE | samtools sort -n --threads $NTHREADS -O SAM -o $FINAL_SAM_FILE -T $OUTPUT_DIR/$RAW_BAM_PREFIX

FINAL_HET_SAM_FILE=${OUTPUT_DIR}/${LIBRARY}.filt.nodup.srt.hetsOnly.sam

module unload cluster/python/3.6.6
module unload cluster/picard/2.21.6
module unload cluster/vg/1.20.0
module unload cluster/samtools/1.8

module load cluster/R/4.1.0

env
Rscript Parsing_Sam_Reads_In_Variant_Regions.R ${FINAL_SAM_FILE} ${FINAL_HET_SAM_FILE} ${ALL_HET_REGIONS}

rm ${RAW_BAM_FILE} ${ALN_GAM} ${ALN_BAM} ${SRT_SAM} ${FILT_BAM_FILE} ${FINAL_BAM_FILE_MAPSTATS} ${FINAL_SAM_FILE} ${FASTQ_FILE_2}
#rm -rf ${TEMP_DIR}
