#!/bin/bash

# Activate general conda
source /usr/local/miniforge/etc/profile.d/conda.sh

conda activate airlom

# Main presentation
cat <<-'PRESENTATION'

######################################################
#                                                    #
#   AIRLoM (Adaptive Immune-Receptor Locus Mapper)   #
#                                                    #
######################################################

PRESENTATION

# Get start date in seconds.
start=$SECONDS

# Load functions.
source ./SCRIPTS/MASTER/functions.sh

# If no argument is given, show error and help.
[ $# -eq 0 ] && echo -e "\nNo arguments given.\n" && usage

# Give default values to every fasta file to be used to create databases
IGH_V_CDNA="./SCRIPTS/SEQ_DB/RAW/IGHV_CDS_7bat_01_11_23.fna"
IGH_C_CDNA="./SCRIPTS/SEQ_DB/RAW/IGHC_F_cDNA.fna"
IGH_V_AA="./SCRIPTS/SEQ_DB/RAW/IGHV_CDS_7bat_01_11_23.faa"
IGH_J_CDNA="./SCRIPTS/SEQ_DB/RAW/IGHJ_7bat.fna"
RSS_V="./SCRIPTS/SEQ_DB/RAW/IGHV_RSS3_7bat.fna"
RSS_D3="./SCRIPTS/SEQ_DB/RAW/IGHD_RSS3_7bat.fna"
RSS_D5="./SCRIPTS/SEQ_DB/RAW/IGHD_RSS5_7bat.fna"
RSS_J="./SCRIPTS/SEQ_DB/RAW/IGHJ_RSS5_7bat.fna"
SIG_P="./SCRIPTS/SEQ_DB/RAW/L1_exon_7bat.fna"
#NOTE: OPTION IN GETOPTS MUST ALSO HAVE A FLAG TO CATCH THE ARGUMENT
while getopts ":hs:g:a:m:i:j:k:l:p:w:x:y:z:" flag
do
        case "${flag}" in
        s | specie) SPECIE=${OPTARG} ;;
                g | genome) GENOME=${OPTARG} ;;
                a | annotation) ANNOTATION=${OPTARG} ;;
                i | igh_v_cdna) IGH_V_CDNA=${OPTARG} ;;
                j | igh_c_cdna) IGH_C_CDNA=${OPTARG} ;;
                k | igh_v_aa) IGH_V_AA=${OPTARG} ;;
                l | igh_j_cdna) IGH_J_CDNA=${OPTARG} ;;
                p | signal_p)  SIG_P=${OPTARG} ;;
                w | rss_v) RSS_V=${OPTARG} ;;
                x | rss_d3) RSS_D3=${OPTARG} ;;
                y | rss_d5) RSS_D5=${OPTARG} ;;
                z | rss_j) RSS_J=${OPTARG} ;;
                h | help) usage ;;
                :) echo -e "\nMissing argument submitted ($OPTARG).\n" && usage && exit 1 ;;
                \?) echo -e "\nInvalid option submitted.\n" && usage && exit 1 ;;
        esac
done

# ----- CHECK ARGUMENTS -----
check_arguments

## ----- FORMAT NAMES -----
format_name 

# ----- RUN CD-HIT -----
run_cdhit 

## ----- REASSING VARIABLES -----
reassign_vars

#at <<-'JUMP' >/dev/null 2>&1

# ----- MOVE GENOME AND MAKE BLAST DATABASE -----
move_genome
make_blast_database

# ----- BLAST -----
blast_v_cdna
blast_v_aa 
blast_j_cdna 
blast_c_cdna

conda deactivate

# ----- REFORMAT M6 TO GFF -----
conda activate R

reformat_blast_v_cdna
reformat_blast_v_aa
reformat_blast_j_cdna
reformat_blast_c_cdna

conda deactivate

# ----- EXTRACT SCAFFOLDS NAMES-----
extract_scaffolds

# ----- EXTRACT SCAFFOLDS SEQUENCE-----
conda activate airlom

extract_scaffolds_seq 

# ----- EXONERATE -----
# Make exonerate directory
mkdir ./RESULTS/$SPECIE/EXONERATE 

exonerate_v_cdna
exonerate_v_aa
exonerate_j_cdna
exonerate_c_cdna

# ----- EXTRACT VULGAR FROM EXONERATE -----
extract_exonerate_vulgar

# ----- EXTRACT GFF FROM EXONERATE -----
extract_exonerate_gff

# ----- CONVERT VULGAR TO TABLE -----
conda deactivate
conda activate R
vulgar_to_table

# Remove raw exonerate
rm ./RESULTS/$SPECIE/EXONERATE/exonerate*$SHORT_GS

# ----- HMMER -----
conda deactivate
conda activate airlom
run_hmmer

# Make exonerate filtered directory
mkdir ./RESULTS/$SPECIE/EXONERATE_FILTERED

# ----- FILTRAE EXONERATE TO SEPARATE EXONS AND GENES, PLUS AND MINUS -----
exonerate_filtration

conda deactivate

# ----- REDUCE EXONERATE FILTERED FILES -----
mkdir ./RESULTS/$SPECIE/REDUCTION

conda activate R
exonerate_reduction

# ----- OVERLAP ANALYSIS -----
mkdir ./RESULTS/$SPECIE/OVERLAP
exonerate_overlap_v_cdna
exonerate_overlap_v_aa

# ----- CORRECT J COORDINATES -----
mkdir -p ./RESULTS/$SPECIE/J_RSS_CORRECTED
correct_j_plus
correct_j_minus

#JUMP

conda activate R

# ----- DETECT D SEGMENTS -----
mkdir -p ./RESULTS/$SPECIE/D_SEGMENTS/
detect_d

# ----- CORRECT V COORDINATES -----
mkdir -p ./RESULTS/$SPECIE/V_RSS_CORRECTED
correct_v_minus
correct_v_plus

# ----- MAKE RESULTS DIRECTORY -----
merge_gffs

# ----- Make BED from final GFF -----
gff_to_bed

# ----- FINAL MESSAGE -----
end=$SECONDS
runtime=$((end-start))
printf "### ANALYSIS COMPLETED IN $runtime SECONDS\n"
printf "### EXITING\n"

