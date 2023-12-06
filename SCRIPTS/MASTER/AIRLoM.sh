#!/bin/bash

# Get start date in seconds
start=$SECONDS

# Load functions
source ./SCRIPTS/MASTER/fun.sh

#Si no hay argumentos, mostrar error y la ayuda
[ $# -eq 0 ] && echo -e "\nNo arguments given.\n" && usage

#OJO EN GETOPTS DEBEN ESTAR TAMBIEN LA ETIQUETA QUE RECIBIRA EL ARGUMENTO
while getopts ":hs:g:a:m:i:j:k:l:w:x:y:z:" flag
do
        case "${flag}" in
        s | specie) SPECIE=${OPTARG} ;;
		g | genome) GENOME=${OPTARG} ;;
		a | annotation) ANNOTATION=${OPTARG} ;;
		i | igh_v_cdna) IGH_V_CDNA=${OPTARG} ;;
		j | igh_c_cdna) IGH_C_CDNA=${OPTARG} ;;
		k | igh_v_aa) IGH_V_AA=${OPTARG} ;;
		l | igh_j_cdna) IGH_J_CDNA=${OPTARG} ;;
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

# ----- RUN CD-HIT -----
#run_cdhit 

# ----- REASSING VARIABLES -----
reassign_vars

# ----- FORMAT NAMES -----
format_name 

# ----- MOVE GENOME AND MAKE BLAST DATABASE -----
#move_genome
#make_blast_database

# ----- BLAST -----
#blast_v_cdna 
#blast_v_aa 
#blast_j_cdna 
#blast_c_cdna

# ----- REFORMAT M6 TO GFF -----
#reformat_blast_v_cdna
#reformat_blast_v_aa
#reformat_blast_j_cdna
#reformat_blast_c_cdna

# ----- EXTRACT SCAFFOLDS NAMES-----
#extract_scaffolds

# ----- EXTRACT SCAFFOLDS SEQUENCE-----
#extract_scaffolds_seq 

# Make exonerate directory
#mkdir ./RESULTS/$SPECIE/EXONERATE 
# ----- EXONERATE -----
#exonerate_v_cdna
#exonerate_v_aa
#exonerate_j_cdna
#exonerate_c_cdna

# ----- EXTRACT VULGAR FROM EXONERATE -----
#extract_exonerate_vulgar

# ----- EXTRACT GFF FROM EXONERATE -----
#extract_exonerate_gff

# ----- CONVERT VULGAR TO TABLE -----
#vulgar_to_table

# Remove raw exonerate
#rm ./RESULTS/$SPECIE/EXONERATE/exonerate*$SHORT_GS

# ----- HMMER -----
#run_hmmer

# Make exonerate filtered directory
#mkdir ./RESULTS/$SPECIE/EXONERATE_FILTERED

# ----- FILTRAE EXONERATE TO SEPARATE EXONS AND GENES, PLUS AND MINUS -----
#exonerate_filtration

# ----- REDUCE EXONERATE FILTERED FILES -----
#mkdir ./RESULTS/$SPECIE/REDUCTION
#exonerate_reduction

# ----- OVERLAP ANALYSIS -----
#mkdir ./RESULTS/$SPECIE/OVERLAP
#exonerate_overlap_v_cdna
#exonerate_overlap_v_aa

# ----- CORRECT J COORDINATES -----
#mkdir -p ./RESULTS/$SPECIE/J_RSS_CORRECTED
#correct_j_plus
#correct_j_minus
#correct_v_plus



# ----- FINAL MESSAGE -----
end=$SECONDS
runtime=$((end-start))
printf "### ANALYSIS COMPLETED IN $runtime SECONDS\n"
printf "### EXITING\n"
