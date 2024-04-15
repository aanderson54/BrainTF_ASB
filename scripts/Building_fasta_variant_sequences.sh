#!/bin/bash
#Building_fasta_variant_sequences.sh

################################################################################
# The script below was used to create a fasta file of all the variant regions in
#     a donor using phased variants in vcf format. It converts that vcf into a
#     table of the passed variants, and then identifies variants within some
#     designated range of one another (HANGOVER_LEN).
#
#     The script below is tailored to our local compute cluster's specific
#     architecture, and will need to be rewritten for your own compute environment's
#     specific needs if you would like to reproduce its results.  In particular,
#     the various sourcings will likely need to be removed, as well as the log_msg
#     commands.  One will need to provide or define their own TEMP_DIR.
#
# Expected to be in the user's PATH:
#
#     R 3.4.3 (or another R version with necessary packages)
#
#
# Required Environment Variables:
#     VCF            - A file in VCF format corresponding to variants found in
#                      the donor of interest. Provided via the following DOI:
#                      https://doi.org/10.7303/syn4921369
#     OUTDIR         - Full path to the directory to which all outputs will be
#                      written.  NOTE that if this directory exists, it will be
#                      overwritten by this pipeline.
#
# Optional Environment Variables:
#     HANGOVER_LEN   - Number of basepairs to be included past each variant.
#                      Min is 101, default is 101.
#     GENOME_FILE    - Path to the reference genome fasta file.  Should match
#                      genome used for deriving VCF.
#                      Default: /gpfs/gpfs1/home/jlawlor/longranger_jacob/refdata-GRCh38-2.1.0/fasta/genome.fa
#
################################################################################

export SOURCE_DIR="/gpfs/gpfs2/sdi/etc"
export PIPELINE_VENV_PYTHON3="/gpfs/gpfs2/sdi/software/pipeline-venv-python3/bin/activate"

### Mandatory sourcing of bashrc for necessary environment variables. ###
if [ -e $SOURCE_DIR/bashrc ]; then
    . $SOURCE_DIR/bashrc
else echo "[fatal] - Could not find myerslab bashrc file. Exiting"; exit 1; fi

### Mandatory sourcing of functions to get helper functions (like call_cmd). ###
if [ -e $SOURCE_DIR/functions ]; then
    . $SOURCE_DIR/functions
else echo "[fatal] - Could not find functions file. Exiting"; exit 1; fi

####Make Outdir
if [ -e $OUTDIR ]; then
    echo "OUTDIR already exists."
    #rm -rf $OUTDIR
else
  echo "Making OUTDIR"
  mkdir $OUTDIR
  mkdir ${OUTDIR}/VariantGenome
fi


if [ -z ${HANGOVER_LEN+x} ]; then
    echo "Hangover length not set.  Setting to 101 bp"
    HANGOVER_LEN="101"
else
    echo "Hangover length set to $HANGOVER_LEN"
fi

if [ -z ${GENOME_FILE+x} ]; then
    echo "Genome File not set.  Setting to /gpfs/gpfs1/home/jlawlor/longranger_jacob/refdata-GRCh38-2.1.0/fasta/genome.fa"
    GENOME_FILE=/gpfs/gpfs1/home/jlawlor/longranger_jacob/refdata-GRCh38-2.1.0/fasta/genome.fa
else
    echo "Genome file set to $GENOME_FILE"
fi


### Activate Python3 vevn ###
if [ -e $PIPELINE_VENV_PYTHON3 ]; then
    . $PIPELINE_VENV_PYTHON3
else echo "[fatal] - Could not find python3 activate file[$PIPELINE_VENV_PYTHON3]. Exiting"; exit 1; fi

export RPATH="/gpfs/gpfs2/software/R-3.4.3"
export R_LIBS_SITE=/gpfs/gpfs2/sdi/software/R-site-packages-3.4.3
export R_SCRIPT_BIN=/gpfs/gpfs2/software/R-3.4.3/bin/
export LOGFILE_NAME=$OUTDIR/MakingDonorGenome.log




### First, we need to build the donor genome fasta file in multiple steps. ###
VCF_TABLE=${OUTDIR}/phased_variants_PASSED.vcf.table

log_msg info "Parsing VCF into a table..."
CMD="${RPATH}/bin/Rscript /gpfs/gpfs1/home/bmoyers/Scripts/BrainTF_AlleleSpecificBinding/VCF_to_VCFTable.R ${VCF} ${VCF_TABLE}"
run_cmd "$CMD" "$VCF_TABLE"

log_msg info "Creating all variant sequences for each chromosome..."
CMD="${RPATH}/bin/Rscript /gpfs/gpfs1/home/bmoyers/Scripts/BrainTF_AlleleSpecificBinding/Construct_fasta_variant_sequences.R ${VCF_TABLE} ${OUTDIR} ${GENOME_FILE} ${HANGOVER_LEN}"
run_cmd "$CMD" "${OUTDIR}/AllVar_Sequences_chr1.fasta"


log_msg info "Combining all variant sequences into one file..."
ALL_VAR_SEQS=${OUTDIR}/AllVar_Sequences.fasta
ALL_PROBLEM_SEQS=${OUTDIR}/ProblemSequences.fasta
ALL_PROBLEM_VARS=${OUTDIR}/ProblemVariants.txt
CMD="cat ${OUTDIR}/AllVar_Sequences_chr*.fasta > ${ALL_VAR_SEQS}"
run_cmd "$CMD" "${ALL_VAR_SEQS}"
rm ${OUTDIR}/AllVar_Sequences_chr*.fasta

CMD="cat ${OUTDIR}/ProblemSequences_chr* > ${ALL_PROBLEM_SEQS}"
run_cmd "$CMD" "${ALL_PROBLEM_SEQS}"
rm ${OUTDIR}/ProblemSequences_chr*

CMD="cat ${OUTDIR}/ProblemVariants_chr* > ${ALL_PROBLEM_VARS}"
run_cmd "$CMD" "${ALL_PROBLEM_VARS}"
rm ${OUTDIR}/ProblemVariants_chr*



deactivate
