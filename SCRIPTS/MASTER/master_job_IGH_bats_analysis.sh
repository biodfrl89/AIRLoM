#!/bin/bash

#$ -S /bin/bash
#$ -N prueba_union_scripts
#$ -cwd
#$ -pe mpi 4
#$ -l h_vmem=1G
#$ -l qname=all.q
#$ -e error
#$ -o output

#Este programa debe ejecutarse en una carpeta que contenga el genoma de referencia en formato fna.gz y su correspondiente gff.gz
#El programa creara los directorios correspondientes y movera los archivos a sus destinos antes de empezar el procesamiento.

#Obtener el nombre del programa solamente para la función usage
#PROGRAM=$(basename $0)

# The script must be executed in ANALISIS/

# Get start date in seconds
start=$SECONDS

# Documentation
usage(){
cat <<-'EOF'

DESCRIPTION:

SYNTAX:
 ./master_job_IGH_bats_analysis.sh [-s|-specie GENUS_SPECIE] [-g|--genome GENOME_FILE]
[-a|--annotation ANNOTATION_FILE] [-i|--igh_v_cdna FASTA_FILE]
[-j|--igh_c_cdna FASTA_FILE] [-k|--igh_v_aa FASTA_FILE] [-w|--rss_v FASTA_FILE]
[-x|--rss_d3 FASTA_FILE] [-y|--rss_d5 FASTA_FILE] [-z|--rss_j FASTA_FILE]

OPTIONS:
-s, --specie		Name of the specie to be analized. Format: GENUS_SPECIE. Mandatory.
-g, --genome		The file that has the compressed genome to be analized. Must be a .fna.gz file. Mandatory.
-a, --annotation	The annotations file of the .fna.gz. It also has to be compressed. Must be a .gff.gz file. Optional.
-i, --igh_v_cdna	A fasta file with the cDNA sequences of IGHV. Mandatory.
-j, --igh_c_cdna	A fasta file with the cDNA sequences of IGHC. Mandatory.
-k, --igh_v_aa		A fasta file with the aminoacid sequences of IGHV. Mandatory.
-w, --rss_v		A fasta file with the RSS sequences for the V genes. Mandatory.
-x, --rss_d3		A fasta file with the RSS sequences for the 3' of D genes. Mandatory.
-y, --rss_d5		A fasta file with the RSS sequences for the 5' of D genes. Mandatory.
-z, --rss_j		A fasta file with the RSS sequences for the J genes. Mandatory

EOF

exit 0

}

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

## GENERAL

#Si la variable SPECIE contiene algun espacio, mostrar el error y la ayuda.
[[ $SPECIE == *[[:space:]]* ]] && echo -e "\nIncorrect species name format. Check help.\n" && exit 1

#Si la variable GENOME esta vacía, mostrar el error y la ayuda
[[ -z $GENOME ]] && echo -e "\nNo genome submitted.\n" && exit 1

#Checar que exista el archivo .fna.gz
[[ ! -f $GENOME ]] && echo -e "\nIncorrect .fna.gz sugmitted. Check help.\n" && exit 1

## IGHV FILES

# Check if VAR is empty.
[[ -z $IGH_V_CDNA ]] && echo -e "\nNo IGHV cDNA fasta file submitted.\n" && exit 1

# Check that file exists.
[[ ! -f $IGH_V_CDNA ]] && echo -e "\nIGHV cDNA fasta file not found.\n" && exit 1

# Check if VAR is empty.
[[ -z $IGH_C_CDNA ]] && echo -e "\nNo IGHC cDNA fasta file submitted.\n" && exit 1

# Check that file exists.
[[ ! -f $IGH_C_CDNA ]] && echo -e "\nIGHC cDNA fasta file not found.\n" && exit 1

# Check if VAR is empty.
[[ -z $IGH_J_CDNA ]] && echo -e "\nNo IGHJ cDNA fasta file submitted.\n" && exit 1

# Check that file exists.
[[ ! -f $IGH_J_CDNA ]] && echo -e "\nIGHJ cDNA fasta file not found.\n" && exit 1

# Check if VAR is empty.
[[ -z $IGH_V_AA ]] && echo -e "\nNo IGHV AA fasta file submitted.\n" && exit 1

# Check that file exists.
[[ ! -f $IGH_V_AA ]] && echo -e "\nIGHV AA fasta file not found.\n" && exit 1

## RSS FILES

# Check if VAR is empty.
[[ -z $RSS_V ]] && echo -e "\nNo RSS_V fasta file submitted.\n" && exit 1

# Check that file exists.
[[ ! -f $RSS_V ]] && echo -e "\nRSS_V fasta file not found.\n" && exit 1

# Check if VAR is empty.
[[ -z $RSS_D3 ]] && echo -e "\nNo RSS_D3 fasta file submitted.\n" && exit 1

# Check that file exists.
[[ ! -f $RSS_D3 ]] && echo -e "\nRSS_D3 fasta file not found.\n" && exit 1

# Check if VAR is empty.
[[ -z $RSS_D5 ]] && echo -e "\nNo RSS_D5 fasta file submitted.\n" && exit 1

# Check that file exists.
[[ ! -f $RSS_D5 ]] && echo -e "\nRSS_D5 fasta file not found.\n" && exit 1

# Check if VAR is empty.
[[ -z $RSS_J ]] && echo -e "\nNo RSS_J fasta file submitted.\n" && exit 1

# Check that file exists.
[[ ! -f $RSS_J ]] && echo -e "\nRSS_J fasta file not found.\n" && exit 1

## OPTIONALS

#Si se suministra la variable annotation, checar que exista el archivo:
if [[ -n $ANNOTATION ]] 
then 
	[[ ! -f $ANNOTATION ]] && echo -e "\nIncorrect .gff.gz sugmitted\n" && exit 1
fi

# ----- CHECKPOINT -----

#Si todo está en orden, se muestran las variables a utilizar
cat<<-EOF

##########################################

Inmunoglobuling Detection Tool v0.3

#########################################

The following arguments were correctly submited

Specie to analize: $SPECIE
Genome file: $GENOME
Annotation file: $ANNOTATION
IGHV cDNA file: $IGH_V_CDNA
IGHV aa file: $IGH_V_AA
IGHJ cDNA file:$IGH_J_CDNA
IGHC cDNA file:	$IGH_C_CDNA
RSS_V file: $RSS_V
RSS_D5 file: $RSS_D5
RSS_D3 file: $RSS_D3
RSS_J file: $RSS_J

Press (y|Y) to begin the analysis or (n|N) to exit.

EOF

cat <<-'EOF' > /dev/null 2>&1

# CHECK KEY PRESSED
while true; do
	# Wait for the user to press a key
	read -s -n 1 key

	# Check which key was pressed
	case $key in
    		y|Y) printf "You pressed 'y'. Continuing...\n\n" && break ;;
    		n|N) printf "You pressed 'n'. Exiting...\n\n" && exit 1 ;;
    		*) printf "Invalid input. Please press 'y' or 'n'.\n\n" ;;
	esac
done
EOF

# -----CD-HIT-----

# Run CD-HIT to cluster redundancy
printf "### Running CD-HIT to cluster raw database sequences redundancy\n"
conda run -n roary cd-hit-est -i $IGH_V_CDNA -o cdhit_IGHV_CDS.fna -c 0.95 -n 5 -M 2000 -T 4 >/dev/null 2>&1
conda run -n roary cd-hit -i $IGH_V_AA -o cdhit_IGHV_CDS.faa -c 0.95 -n 5 -M 2000 -T 4 >/dev/null 2>&1
conda run -n roary cd-hit-est -i $IGH_J_CDNA -o cdhit_IGHJ.fna -c 0.95 -n 5 -M 2000 -T 4 >/dev/null 2>&1
conda run -n roary cd-hit-est -i $IGH_C_CDNA -o cdhit_IGHC.fna -c 0.95 -n 5 -M 2000 -T 4 >/dev/null 2>&1
conda run -n roary cd-hit-est -i $RSS_V -o cdhit_IGHV_RSS3.fna -c 0.95 -n 5 -M 2000 -T 4 >/dev/null 2>&1
conda run -n roary cd-hit-est -i $RSS_D5 -o cdhit_IGHD_RSS5.fna -c 0.95 -n 5 -M 2000 -T 4 >/dev/null 2>&1
conda run -n roary cd-hit-est -i $RSS_D3 -o cdhit_IGHD_RSS3.fna -c 0.95 -n 5 -M 2000 -T 4 >/dev/null 2>&1
conda run -n roary cd-hit-est -i $RSS_J -o cdhit_IGHJ_RSS5.fna -c 0.95 -n 5 -M 2000 -T 4 >/dev/null 2>&1

# Run CD-HIT to cluster redundancy
printf "### Running CLUSTAL alignment in RSS clustered sequences\n"
conda run -n roary clustalw2 -align -quiet -infile=cdhit_IGHV_RSS3.fna -OUTFILE=clustal_IGHV_RSS3.fna -OUTPUT=FASTA >/dev/null 2>&1
conda run -n roary clustalw2 -align -quiet -infile=cdhit_IGHD_RSS5.fna -OUTFILE=clustal_IGHD_RSS5.fna -OUTPUT=FASTA >/dev/null 2>&1
conda run -n roary clustalw2 -align -quiet -infile=cdhit_IGHD_RSS3.fna -OUTFILE=clustal_IGHD_RSS3.fna -OUTPUT=FASTA >/dev/null 2>&1
conda run -n roary clustalw2 -align -quiet -infile=cdhit_IGHJ_RSS5.fna -OUTFILE=clustal_IGHJ_RSS5.fna -OUTPUT=FASTA >/dev/null 2>&1

printf "### Moving database files\n"
rm cdhit*.dnd
rm cdhit*RSS?.fna
mv clustal*.fna ./SCRIPTS/SEQ_DB/RSS
mv cdhit*.fna cdhit*.faa ./SCRIPTS/SEQ_DB/FASTAS
mv cdhit*.clstr ./SCRIPTS/SEQ_DB/FASTAS/CLUSTERS

# Reassign variables
printf "### Reassigning variables\n"
IGH_V_CDNA="./SCRIPTS/SEQ_DB/FASTAS/cdhit_IGHV_CDS.fna"
IGH_V_AA="./SCRIPTS/SEQ_DB/FASTAS/cdhit_IGHV_CDS.faa"
IGH_J_CDNA="./SCRIPTS/SEQ_DB/FASTAS/cdhit_IGHJ.fna"
IGH_C_CDNA="./SCRIPTS/SEQ_DB/FASTAS/cdhit_IGHC.fna"
RSS_V="./SCRIPTS/SEQ_DB/RSS/clustal_IGHV_RSS3.fna"
RSS_D5="./SCRIPTS/SEQ_DB/RSS/clustal_IGHD_RSS5.fna"
RSS_D3="./SCRIPTS/SEQ_DB/RSS/clustal_IGHD_RSS3.fna"
RSS_J="./SCRIPTS/SEQ_DB/RSS/clustal_IGHJ_RSS5.fna"

# -----SPECIE-----

# Create short specie name
printf "### Creating short specie name\n"
SPECIE_L=${SPECIE,,} #To lower all
SHORT_GENUS=${SPECIE_L:0:2}
SPECIE_NAME=$(echo $SPECIE_L | cut -d "_" -f2)
SHORT_SPECIE=${SPECIE_NAME:0:2}
SHORT_GS=${SHORT_GENUS^}${SHORT_SPECIE} #Capitalize genus.

# -----BLAST-----

# Create specie directory with reference genome directory. Move genome to reference genome folder and locate inside.
mkdir -p ./RESULTS/$SPECIE/REFERENCE_GENOME && \
mv *.gz ./RESULTS/$SPECIE/REFERENCE_GENOME/ && \
cd ./RESULTS/$SPECIE/REFERENCE_GENOME/

# Decompress genome and pipe it to makeblastdb
printf "### Making BLAST Database\n"

gunzip -dc $GENOME | makeblastdb -in - -dbtype nucl -out BLAST_DB -title $GENOME >/dev/null 2>&1

# Move behind and make blast directory and database directory.
cd .. && mkdir -p BLAST/DB 

# Move to DB directory and move all BLAST file to that location.
cd ./BLAST/DB/ && mv ../../REFERENCE_GENOME/BLAST* .

# Return to original path.
cd ../../../..

# V region with cDNA fasta
printf "### Running BLAST for V regions using cDNA\n"
conda run -n murcielagos tblastx \
-query $IGH_V_CDNA \
-db ./RESULTS/$SPECIE/BLAST/DB/BLAST_DB \
-outfmt 6 \
-evalue 1e-20 \
-num_threads 12 \
-max_target_seqs 100 \
-out ./RESULTS/$SPECIE/BLAST/tblastx_IGHV_cDNA_vs_$SHORT_GS.m6

# V region with aa fasta
printf "### Running BLAST for V regions using AA\n"
conda run -n murcielagos tblastn \
-query $IGH_V_AA \
-db ./RESULTS/$SPECIE/BLAST/DB/BLAST_DB \
-outfmt 6 \
-evalue 1e-20 \
-num_threads 12 \
-max_target_seqs 100 \
-out ./RESULTS/$SPECIE/BLAST/tblastn_IGHV_aa_vs_$SHORT_GS.m6

# J region with cDNA fasta
printf "### Running BLAST for J regions using cDNA\n"
conda run -n murcielagos tblastx \
-query $IGH_J_CDNA \
-db ./RESULTS/$SPECIE/BLAST/DB/BLAST_DB \
-outfmt 6 \
-evalue 1e-4 \
-num_threads 12 \
-max_target_seqs 100 \
-out ./RESULTS/$SPECIE/BLAST/tblastx_IGHJ_cDNA_vs_$SHORT_GS.m6

# C region with cDNA fasta
printf "### Running BLAST for C regions using cDNA\n"
conda run -n murcielagos tblastx \
-query $IGH_C_CDNA \
-db ./RESULTS/$SPECIE/BLAST/DB/BLAST_DB \
-outfmt 6 \
-evalue 1e-20 \
-num_threads 12 \
-max_target_seqs 100 \
-out ./RESULTS/$SPECIE/BLAST/tblastx_IGHC_cDNA_vs_$SHORT_GS.m6

# Use subscript to format blast m6 to gff

printf "### Reformating blast m6 tables of IGHV_cDNA to gff\n"
#conda run -n murcielagos \
#python3 ./SCRIPTS/SUBSCRIPTS/blast_m6_to_gff.py --file ./RESULTS/$SPECIE/BLAST/tblastx_IGHV_cDNA_vs_$SHORT_GS.m6 --source tblastx --bitscore 0
Rscript ./SCRIPTS/SUBSCRIPTS/blast_m6_to_gff.R --file ./RESULTS/$SPECIE/BLAST/tblastx_IGHV_cDNA_vs_$SHORT_GS.m6 --source tblastx --bitscore 0

printf "### Reformating blast m6 tables of IGHV_aa to gff\n"
#conda run -n murcielagos \
#python3 ./SCRIPTS/SUBSCRIPTS/blast_m6_to_gff.py --file ./RESULTS/$SPECIE/BLAST/tblastn_IGHV_aa_vs_$SHORT_GS.m6 --source tblastn --bitscore 0
Rscript ./SCRIPTS/SUBSCRIPTS/blast_m6_to_gff.R --file ./RESULTS/$SPECIE/BLAST/tblastn_IGHV_aa_vs_$SHORT_GS.m6 --source tblastn --bitscore 0

printf "### Reformating blast m6 tables of IGHJ_cDNA to gff\n"
#conda run -n murcielagos \
#python3 ./SCRIPTS/SUBSCRIPTS/blast_m6_to_gff.py --file ./RESULTS/$SPECIE/BLAST/tblastx_IGHJ_cDNA_vs_$SHORT_GS.m6 --source tblastx --bitscore 0
Rscript ./SCRIPTS/SUBSCRIPTS/blast_m6_to_gff.R --file ./RESULTS/$SPECIE/BLAST/tblastx_IGHJ_cDNA_vs_$SHORT_GS.m6 --source tblastx --bitscore 0

printf "### Reformating blast m6 tables of IGHC_cDNA to gff\n"
#conda run -n murcielagos \
#python3 ./SCRIPTS/SUBSCRIPTS/blast_m6_to_gff.py --file ./RESULTS/$SPECIE/BLAST/tblastx_IGHC_cDNA_vs_$SHORT_GS.m6 --source tblastx --bitscore 0
Rscript ./SCRIPTS/SUBSCRIPTS/blast_m6_to_gff.R --file ./RESULTS/$SPECIE/BLAST/tblastx_IGHC_cDNA_vs_$SHORT_GS.m6 --source tblastx --bitscore 0

# Obtain scaffolds ID in each m6 file. Merge the result, and get scaffolds with match obtained in every BLAST. 
printf "### Getting scaffolds ID and extract their sequences\n"
cut -f2 ./RESULTS/$SPECIE/BLAST/tblastx_IGHV_cDNA_vs_$SHORT_GS.m6 | sort | uniq >./RESULTS/$SPECIE/BLAST/scaffolds_tblastx_IGHV_cDNA_vs_$SPECIE
cut -f2 ./RESULTS/$SPECIE/BLAST/tblastn_IGHV_aa_vs_$SHORT_GS.m6 | sort | uniq >./RESULTS/$SPECIE/BLAST/scaffolds_tblastn_IGHV_vs_$SPECIE
cut -f2 ./RESULTS/$SPECIE/BLAST/tblastx_IGHJ_cDNA_vs_$SHORT_GS.m6 | sort | uniq >./RESULTS/$SPECIE/BLAST/scaffolds_tblastx_IGHJ_cDNA_vs_$SPECIE
cut -f2 ./RESULTS/$SPECIE/BLAST/tblastx_IGHC_cDNA_vs_$SHORT_GS.m6 | sort | uniq >./RESULTS/$SPECIE/BLAST/scaffolds_tblastx_IGHC_cDNA_vs_$SPECIE
cat ./RESULTS/$SPECIE/BLAST/scaffolds_t* | sort | uniq -c | awk ' $1 >= 3 ' | tr -s ' ' | cut -d' ' -f3 >./RESULTS/$SPECIE/BLAST/scaffolds_to_extract

# Decompress genome to prepare for sequence extraction
printf "### Decompressing genome\n"
pigz -dk -p8 ./RESULTS/$SPECIE/REFERENCE_GENOME/$GENOME

# Extract the desired scaffold sequences
printf "### Getting sequences of the matching scaffolds\n"
samtools faidx ./RESULTS/$SPECIE/REFERENCE_GENOME/$(echo $GENOME | sed 's/\.gz$//') $(cat ./RESULTS/$SPECIE/BLAST/scaffolds_to_extract) >./RESULTS/$SPECIE/BLAST/scaffolds_extracted

# Index sequences
printf "### Indexing obtained sequences\n"
samtools faidx ./RESULTS/$SPECIE/BLAST/scaffolds_extracted

# Remove decompressed genome
printf "### Removing decompressed genome\n"
rm $(echo "./RESULTS/$SPECIE/REFERENCE_GENOME/$GENOME" | sed 's/\.gz$//')

#-----EXONERATE-----

# Make exonerate directory
mkdir ./RESULTS/$SPECIE/EXONERATE 

# V region with cDNA fasta
printf "### Running EXONERATE for V regions using cDNA\n"
conda run -n murcielagos exonerate \
--cores 4 \
--model est2genome \
--showtargetgff TRUE \
--verbose 0 \
--showalignment no \
--showvulgar yes \
--bestn 30  \
--softmasktarget TRUE \
--maxintron 1000 \
--score 1000 \
-q $IGH_V_CDNA \
-t ./RESULTS/$SPECIE/BLAST/scaffolds_extracted >./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_cDNA_vs_$SHORT_GS

# V region with AA fasta
printf "### Running EXONERATE for V regions using AA \n"
conda run -n murcielagos exonerate \
--cores 2 \
--model protein2genome \
--showtargetgff TRUE \
--verbose 0 \
--showalignment no \
--showvulgar yes \
--bestn 30  \
--score 300 \
--softmasktarget TRUE \
--maxintron 1000 \
-q $IGH_V_AA \
-t ./RESULTS/$SPECIE/BLAST/scaffolds_extracted >./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_vs_$SHORT_GS

# J region with cDNA fasta
printf "### Running EXONERATE for J regions using cDNA\n"

##### Agregar score despues de inspeccio manua
conda run -n murcielagos exonerate \
--cores 4 \
--model est2genome \
--showtargetgff TRUE \
--verbose 0 \
--showalignment no \
--showvulgar yes \
--bestn 30  \
--softmasktarget TRUE \
--maxintron 1000 \
-q $IGH_J_CDNA \
-t ./RESULTS/$SPECIE/BLAST/scaffolds_extracted >./RESULTS/$SPECIE/EXONERATE/exonerate_IGHJ_cDNA_vs_$SHORT_GS

# C region with cDNA fasta
printf "### Running EXONERATE for C regions using cDNA\n"

######### Cambiar intron de la region constante
conda run -n murcielagos exonerate \
--cores 4 \
--model est2genome \
--showtargetgff TRUE \
--verbose 0 \
--showalignment no \
--showvulgar yes \
--bestn 30  \
--softmasktarget TRUE \
--maxintron 1000 \
--score 500 \
-q $IGH_C_CDNA \
-t ./RESULTS/$SPECIE/BLAST/scaffolds_extracted >./RESULTS/$SPECIE/EXONERATE/exonerate_IGHC_cDNA_vs_$SHORT_GS

# Extract the vulgar anotation from from exonerate file for further procesing
printf "### Extracting vulgar anotation from exonerate output file\n"
grep "vulgar" ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_cDNA_vs_$SHORT_GS | sed 's/vulgar: //' >./RESULTS/$SPECIE/EXONERATE/vulgar_IGHV_cDNA_vs_$SHORT_GS
grep "vulgar" ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_vs_$SHORT_GS | sed 's/vulgar: //' >./RESULTS/$SPECIE/EXONERATE/vulgar_IGHV_vs_$SHORT_GS
grep "vulgar" ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHC_cDNA_vs_$SHORT_GS | sed 's/vulgar: //' >./RESULTS/$SPECIE/EXONERATE/vulgar_IGHC_cDNA_vs_$SHORT_GS
grep "vulgar" ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHJ_cDNA_vs_$SHORT_GS | sed 's/vulgar: //' >./RESULTS/$SPECIE/EXONERATE/vulgar_IGHJ_cDNA_vs_$SHORT_GS

# Clean exonerate output files to leave a gff formated file
printf "### Making exonerate gff output file\n"
grep "#" -v  ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_cDNA_vs_$SHORT_GS | grep "vulgar" -v >./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_cDNA_vs_${SHORT_GS}.gff
grep "#" -v  ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_vs_$SHORT_GS | grep "vulgar" -v >./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_vs_${SHORT_GS}.gff
grep "#" -v  ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHC_cDNA_vs_$SHORT_GS | grep "vulgar" -v >./RESULTS/$SPECIE/EXONERATE/exonerate_IGHC_cDNA_vs_${SHORT_GS}.gff
grep "#" -v  ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHJ_cDNA_vs_$SHORT_GS | grep "vulgar" -v >./RESULTS/$SPECIE/EXONERATE/exonerate_IGHJ_cDNA_vs_${SHORT_GS}.gff

# Create temp vulgar files with semicolon separator for reading in the next R script 
awk  -F " " ' function SUB(F) {sub("$", ";", $F)}; {SUB(1); SUB(2) SUB(3) SUB(4) SUB(5) SUB(6) SUB(7) SUB(8) SUB(9)} 1 ' \
 <./RESULTS/$SPECIE/EXONERATE/vulgar_IGHV_cDNA_vs_$SHORT_GS >./RESULTS/$SPECIE/EXONERATE/temp_vulgar_IGHV_cDNA_vs_$SHORT_GS
 awk  -F " " ' function SUB(F) {sub("$", ";", $F)}; {SUB(1); SUB(2) SUB(3) SUB(4) SUB(5) SUB(6) SUB(7) SUB(8) SUB(9)} 1 ' \
 <./RESULTS/$SPECIE/EXONERATE/vulgar_IGHV_vs_$SHORT_GS >./RESULTS/$SPECIE/EXONERATE/temp_vulgar_IGHV_vs_$SHORT_GS
awk  -F " " ' function SUB(F) {sub("$", ";", $F)}; {SUB(1); SUB(2) SUB(3) SUB(4) SUB(5) SUB(6) SUB(7) SUB(8) SUB(9)} 1 ' \
 <./RESULTS/$SPECIE/EXONERATE/vulgar_IGHJ_cDNA_vs_$SHORT_GS >./RESULTS/$SPECIE/EXONERATE/temp_vulgar_IGHJ_cDNA_vs_$SHORT_GS
awk  -F " " ' function SUB(F) {sub("$", ";", $F)}; {SUB(1); SUB(2) SUB(3) SUB(4) SUB(5) SUB(6) SUB(7) SUB(8) SUB(9)} 1 ' \
 <./RESULTS/$SPECIE/EXONERATE/vulgar_IGHC_cDNA_vs_$SHORT_GS >./RESULTS/$SPECIE/EXONERATE/temp_vulgar_IGHC_cDNA_vs_$SHORT_GS

# Reformat temp files to count the elements in the vulgar format
printf "### Making exonerate vulgar table files\n"
Rscript ./SCRIPTS/SUBSCRIPTS/vulgar_to_table.R --file ./RESULTS/$SPECIE/EXONERATE/temp_vulgar_IGHV_cDNA_vs_$SHORT_GS
Rscript ./SCRIPTS/SUBSCRIPTS/vulgar_to_table.R --file ./RESULTS/$SPECIE/EXONERATE/temp_vulgar_IGHV_vs_$SHORT_GS
Rscript ./SCRIPTS/SUBSCRIPTS/vulgar_to_table.R --file ./RESULTS/$SPECIE/EXONERATE/temp_vulgar_IGHC_cDNA_vs_$SHORT_GS
Rscript ./SCRIPTS/SUBSCRIPTS/vulgar_to_table.R --file ./RESULTS/$SPECIE/EXONERATE/temp_vulgar_IGHJ_cDNA_vs_$SHORT_GS

# Remove temp files
rm ./RESULTS/$SPECIE/EXONERATE/temp_*

#-----HMMER-----

# Make HMMER directory
mkdir -p ./RESULTS/$SPECIE/HMMER/DB

# Build hmmer profiles
printf "### Building hmmer databases ###\n"
hmmbuild --dna ./RESULTS/$SPECIE/HMMER/DB/hmmprofile_RSS_IGHV.hmm $RSS_V >/dev/null 2>&1
hmmbuild --dna ./RESULTS/$SPECIE/HMMER/DB/hmmprofile_RSS_IGHD_5.hmm $RSS_D5 >/dev/null 2>&1
hmmbuild --dna ./RESULTS/$SPECIE/HMMER/DB/hmmprofile_RSS_IGHD_3.hmm $RSS_D3 >/dev/null 2>&1
hmmbuild --dna ./RESULTS/$SPECIE/HMMER/DB/hmmprofile_RSS_IGHJ.hmm $RSS_J >/dev/null 2>&1

RSS_DIRS="RSS_IGHV RSS_IGHD_5 RSS_IGHD_3 RSS_IGHJ"

# Prepare profile database for hmmscan
printf "### Pressing hmmer databases for hmmerscan\n"
for dir in $RSS_DIRS
do
	hmmpress ./RESULTS/$SPECIE/HMMER/DB/hmmprofile_$dir.hmm >/dev/null 2>&1
done

# Run nhmmerscan analysis.
for dir in $RSS_DIRS
do
	# Make output directory
	mkdir ./RESULTS/$SPECIE/HMMER/$dir

	# Run nhmmerscan
	printf "### Running nhmmerscan for $dir \n"
    nhmmscan \
    -o ./RESULTS/$SPECIE/HMMER/$dir/nhmmer_$dir \
    --tblout ./RESULTS/$SPECIE/HMMER/$dir/nhmmer_$dir.tbl \
    --cpu 4 \
    --nobias \
    -E 10 \
    ./RESULTS/$SPECIE/HMMER/DB/hmmprofile_$dir.hmm \
    ./RESULTS/$SPECIE/BLAST/scaffolds_extracted

	# Reformat output file to gff 
	conda run -n murcielagos python3 ./SCRIPTS/SUBSCRIPTS/hmmer_tbl_to_gff.py --file ./RESULTS/$SPECIE/HMMER/$dir/nhmmer_$dir.tbl --source hmmerscan --bitscore 0 

	# Rename output file
	mv ./RESULTS/$SPECIE/HMMER/$dir/nhmmer_$dir.gff ./RESULTS/$SPECIE/HMMER/$dir/nhmmer_${dir}_${SHORT_GS}.gff 
done

#-----MOVE RESULTS-----
# Copy all gff of interest to GFF_RESULTS
# Make GFF_RESULTS directory
printf "### Making gff results directory \n"
mkdir ./RESULTS/$SPECIE/GFF_RESULTS

printf "### Copying gff files into gff results directory \n"
# Copy BLAST gff
cp ./RESULTS/$SPECIE/BLAST/*.gff ./RESULTS/$SPECIE/GFF_RESULTS
# Copy EXONERATE
cp ./RESULTS/$SPECIE/EXONERATE/*.gff ./RESULTS/$SPECIE/GFF_RESULTS
# Copy RSS
for dir in $RSS_DIRS
do
	cp ./RESULTS/$SPECIE/HMMER/$dir/*.gff ./RESULTS/$SPECIE/GFF_RESULTS
done

#-----REDUCE-----

# Filter exonerate results to get exons hits.
cat ./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHV_cDNA_vs_${SHORT_GS}.gff | awk ' $3 ~ "exon"  ' >./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHV_cDNA_vs_${SHORT_GS}_exons.gff
cat ./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHV_vs_${SHORT_GS}.gff | awk ' $3 ~ "exon"  ' >./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHV_vs_${SHORT_GS}_exons.gff 
cat ./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHC_cDNA_vs_${SHORT_GS}.gff | awk ' $3 ~ "exon"  ' >./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHC_cDNA_vs_${SHORT_GS}_exons.gff
cat ./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHJ_cDNA_vs_${SHORT_GS}.gff | awk ' $3 ~ "exon"  ' >./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHJ_cDNA_vs_${SHORT_GS}_exons.gff

# Filter exonerate results to get gene hits
cat ./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHV_cDNA_vs_${SHORT_GS}.gff | awk ' $3 ~ "gene"  ' >./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHV_cDNA_vs_${SHORT_GS}_gene.gff
cat ./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHV_vs_${SHORT_GS}.gff | awk ' $3 ~ "gene"  ' >./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHV_vs_${SHORT_GS}_gene.gff 
cat ./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHC_cDNA_vs_${SHORT_GS}.gff | awk ' $3 ~ "gene"  ' >./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHC_cDNA_vs_${SHORT_GS}_gene.gff
cat ./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHJ_cDNA_vs_${SHORT_GS}.gff | awk ' $3 ~ "gene"  ' >./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHJ_cDNA_vs_${SHORT_GS}_gene.gff

# Make directory to perform the overlap reduction.
printf "### Making directory for overlap reduction analysis \n"
mkdir ./RESULTS/$SPECIE/OVERLAP_REDUCTION

# Make reduction
printf "### Runing Rscript for overlap reduction over exonerate exons \n"
Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHV_cDNA_vs_${SHORT_GS}_exons.gff >/dev/null 2>&1
Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHV_vs_${SHORT_GS}_exons.gff >/dev/null 2>&1
Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHC_cDNA_vs_${SHORT_GS}_exons.gff >/dev/null 2>&1
Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHJ_cDNA_vs_${SHORT_GS}_exons.gff >/dev/null 2>&1

printf "### Runing Rscript for overlap reduction over exonerate genes \n"
Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHV_cDNA_vs_${SHORT_GS}_gene.gff >/dev/null 2>&1
Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHV_vs_${SHORT_GS}_gene.gff >/dev/null 2>&1
Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHC_cDNA_vs_${SHORT_GS}_gene.gff >/dev/null 2>&1
Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/GFF_RESULTS/exonerate_IGHJ_cDNA_vs_${SHORT_GS}_gene.gff >/dev/null 2>&1

# Move produced files into the overlap reduction directory
mv ./RESULTS/$SPECIE/GFF_RESULTS/reduced_exonerate* ./RESULTS/$SPECIE/OVERLAP_REDUCTION/


#-----REDUCE-----

# Make overlap analysis for IGHV using cDNA
printf "### Runing overlap analysis for IGHV using cDNA \n"
Rscript ./SCRIPTS/SUBSCRIPTS/predict_ighv_by_overlaps.R \
--query ./RESULTS/$SPECIE/OVERLAP_REDUCTION/reduced_exonerate_IGHV_cDNA_vs_${SHORT_GS}_exons.gff \
--subject ./RESULTS/$SPECIE/OVERLAP_REDUCTION/reduced_exonerate_IGHV_cDNA_vs_${SHORT_GS}_gene.gff >/dev/null 2>&1

mv ./gene_prediction.gff ./RESULTS/$SPECIE/OVERLAP_REDUCTION/overlap_gene_prediction_IGHV_cDNA_vs_${SHORT_GS}.gff 
mv ./exon_prediction.gff ./RESULTS/$SPECIE/OVERLAP_REDUCTION/overlap_exon_prediction_IGHV_cDNA_vs_${SHORT_GS}.gff 

# Make overlap analysis for IGHV using amino acids
printf "### Runing overlap analysis for IGHV using aa \n"
Rscript  ./SCRIPTS/SUBSCRIPTS/predict_ighv_by_overlaps.R \
--query ./RESULTS/$SPECIE/OVERLAP_REDUCTION/reduced_exonerate_IGHV_vs_${SHORT_GS}_exons.gff \
--subject ./RESULTS/$SPECIE/OVERLAP_REDUCTION/reduced_exonerate_IGHV_vs_${SHORT_GS}_gene.gff >/dev/null 2>&1

mv ./gene_prediction.gff ./RESULTS/$SPECIE/OVERLAP_REDUCTION/overlap_gene_prediction_IGHV_vs_${SHORT_GS}.gff 
mv ./exon_prediction.gff ./RESULTS/$SPECIE/OVERLAP_REDUCTION/overlap_exon_prediction_IGHV_vs_${SHORT_GS}.gff 

# Make overlap analysis for IGHJ using cDNA
#printf "### Runing overlap analysis for IGHJ using cDNA \n"
#Rscript  ./SCRIPTS/SUBSCRIPTS/predict_ighv_by_overlaps.R \
#--query ./RESULTS/$SPECIE/OVERLAP_REDUCTION/reduced_exonerate_IGHJ_cDNA_vs_${SHORT_GS}_exons.gff \
#--subject ./RESULTS/$SPECIE/OVERLAP_REDUCTION/reduced_exonerate_IGHJ_cDNA_vs_${SHORT_GS}_gene.gff >/dev/null 2>&1

#mv ./gene_prediction.gff ./RESULTS/$SPECIE/OVERLAP_REDUCTION/overlap_gene_prediction_IGHJ_vs_${SHORT_GS}.gff 
#mv ./exon_prediction.gff ./RESULTS/$SPECIE/OVERLAP_REDUCTION/overlap_exon_prediction_IGHJ_vs_${SHORT_GS}.gff 

# Copying GFF file to GFF_RESULTS
printf "### Copying Overlap GFF files to gff results directory \n"
cp ./RESULTS/$SPECIE/OVERLAP_REDUCTION/overlap_* ./RESULTS/$SPECIE/GFF_RESULTS

#-----MINIPROT-----

# Make miniprot results directory
printf "### Making miniprot directory \n"
mkdir ./RESULTS/$SPECIE/MINIPROT

# Index genome for miniprot
printf "### Indexing scaffold for miniprot analysis \n"
conda run -n murcielagos \
miniprot -t8 -d ./RESULTS/$SPECIE/BLAST/scaffolds_extracted.mpi ./RESULTS/$SPECIE/BLAST/scaffolds_extracted >/dev/null 2>&1

# Run miniprot
printf "### Runing miniprot \n"
conda run -n murcielagos \
miniprot -t 8 --aln --gff -G 1000 ./RESULTS/$SPECIE/BLAST/scaffolds_extracted.mpi $IGH_V_AA >./RESULTS/$SPECIE/MINIPROT/miniprot_IGHV_vs_${SHORT_GS}.paf 

# Identify complete ORFs and saving to table
printf "### Making miniprot target results table \n"
cat ./RESULTS/$SPECIE/MINIPROT/miniprot_IGHV_vs_${SHORT_GS}.paf | grep "##AQA" | sed 's/^##AQA\t//' | sed 's/ //g'| awk '{ if ($1 ~ "-") {print $0,"\t", "FALSE"} else {print $0, "\t", "TRUE"}  } ' >./RESULTS/$SPECIE/MINIPROT/miniprot_IGHV_vs_${SHORT_GS}_table_orfs.tsv

#TEST
# cat miniprot_IGHV_vs_Arja.paf |  grep "##ATA" | sed 's/^##ATA\t//' | sed 's/ //g'| awk '{ if ($1 ~ "*") {print $0,"\t", "FALSE"} else {print $0, "\t", "TRUE"}  } '
# cat miniprot_IGHV_vs_Arja.paf |  grep "##ATA" | sed 's/^##ATA\t//' | sed 's/ //g'| sed 's/\.//g' | uniq | less

# Get GFF with miniprot results
printf "### Making GFF of miniprot results \n"
cat ./RESULTS/$SPECIE/MINIPROT/miniprot_IGHV_vs_${SHORT_GS}.paf | grep -v "##" >./RESULTS/$SPECIE/MINIPROT/miniprot_IGHV_vs_${SHORT_GS}.gff

# Copying GFF file to GFF_RESULTS
printf "### Copying miniprot GFF file to gff results directory \n"
cp ./RESULTS/$SPECIE/MINIPROT/miniprot_IGHV_vs_${SHORT_GS}.gff ./RESULTS/$SPECIE/GFF_RESULTS

# Remove indexed genome
printf "### Removing indexed scaffolds \n"
rm ./RESULTS/$SPECIE/BLAST/scaffolds_extracted.mpi

# Make reduction of miniprot gff
printf "### Reducing miniprot gff results\n"
Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/MINIPROT/miniprot_IGHV_vs_${SHORT_GS}.gff >/dev/null 2>&1

#-----MERGE RSS AND VDJ SEGMENTS-----

# Make results directory
printf "### Making V_RSS_CORRECTED_GFF directory \n"
mkdir -p ./RESULTS/$SPECIE/V_RSS_CORRECTED_GFF

# Make overlap search between V segments and RSS and correct joining coordinates 
printf "### Making overlap search between V segments and V RSS, and correct joining coordinates.\n"
Rscript ./SCRIPTS/SUBSCRIPTS/locate_nearby_rss.R \
-n ./RESULTS/$SPECIE/HMMER/RSS_IGHV/nhmmer_RSS_IGHV_${SHORT_GS}.gff \
-t ./RESULTS/$SPECIE/HMMER/RSS_IGHV/nhmmer_RSS_IGHV.tbl \
-v ./RESULTS/$SPECIE/OVERLAP_REDUCTION/overlap_gene_prediction_IGHV_vs_${SHORT_GS}.gff \
-r 50 \
-m V_segments >/dev/null 2>&1

mv v_rss_analysis.gff ./RESULTS/$SPECIE/V_RSS_CORRECTED_GFF/IGHV_RSS_analysis_${SHORT_GS}.gff

# Make results directory
printf "### Making J_RSS_CORRECTED_GFF directory \n"
mkdir -p ./RESULTS/$SPECIE/J_RSS_CORRECTED_GFF

# Make overlap search between J segments and RSS and correct joining coordinates 
# REMEMBER FOR J SEGMENTS USE reduced_exonerate_IGHJ_cDNA_vs_${SHORT_GS}_exons
printf "### Making overlap search between V segments and V RSS, and correct joining coordinates.\n"
Rscript ./SCRIPTS/SUBSCRIPTS/locate_nearby_rss.R \
-m J_segments \
-n ./RESULTS/$SPECIE/HMMER/RSS_IGHJ/nhmmer_RSS_IGHJ_${SHORT_GS}.gff \
-t ./RESULTS/$SPECIE/HMMER/RSS_IGHJ/nhmmer_RSS_IGHJ.tbl \
-v ./RESULTS/$SPECIE/OVERLAP_REDUCTION/reduced_exonerate_IGHJ_cDNA_vs_${SHORT_GS}_exons.gff \
-r 50 >/dev/null 2>&1

mv j_rss_analysis.gff ./RESULTS/$SPECIE/J_RSS_CORRECTED_GFF/IGHJ_RSS_analysis_${SHORT_GS}.gff

#-----DETECT D SEGMENTS SEQUENCE USING D RSS-----
mkdir -p ./RESULTS/$SPECIE/D_SEGMENTS/PREDICTION

printf "### Predict D segments using the D-RSS sequences.\n"

Rscript ./SCRIPTS/SUBSCRIPTS/predict_ighd_by_rss.R \
-f ./RESULTS/$SPECIE/GFF_RESULTS/nhmmer_RSS_IGHD_5_${SHORT_GS}.gff \
-t ./RESULTS/$SPECIE/GFF_RESULTS/nhmmer_RSS_IGHD_3_${SHORT_GS}.gff >/dev/null 2>&1

mv RSS_D_IGH_D_segments.gff ./RESULTS/$SPECIE/D_SEGMENTS/PREDICTION/D_segments_and_RSS_${SHORT_GS}.gff

# Clean gff
cd ./RESULTS/$SPECIE/D_SEGMENTS/PREDICTION/
grep -v "#" D_segments_and_RSS_${SHORT_GS}.gff >temp.gff && \
rm D_segments_and_RSS_${SHORT_GS}.gff && \
mv temp.gff D_segments_and_RSS_${SHORT_GS}.gff
cd ../../../..

# -----MAKE FINAL DIRECTORY WITH CONDENSED RESULTS----- 
printf "### Making final condensed gff.\n"

mkdir -p ./RESULTS/$SPECIE/CONDENSED_RESULTS

cat ./RESULTS/$SPECIE/V_RSS_CORRECTED_GFF/IGHV_RSS_analysis_${SHORT_GS}.gff \
./RESULTS/$SPECIE/J_RSS_CORRECTED_GFF/IGHJ_RSS_analysis_${SHORT_GS}.gff \
./RESULTS/$SPECIE/D_SEGMENTS/PREDICTION/D_segments_and_RSS_${SHORT_GS}.gff \
./RESULTS/$SPECIE/OVERLAP_REDUCTION/overlap_gene_prediction_IGHV_cDNA_vs_${SHORT_GS}.gff \
./RESULTS/$SPECIE/OVERLAP_REDUCTION/overlap_exon_prediction_IGHV_cDNA_vs_${SHORT_GS}.gff \
./RESULTS/$SPECIE/EXONERATE/exonerate_IGHC_cDNA_vs_${SHORT_GS}.gff | \
grep -v "#" >./RESULTS/$SPECIE/CONDENSED_RESULTS/${SHORT_GS}_cond.gff

#-----END-----

end=$SECONDS
runtime=$((end-start))

# Ending Program
printf "### ANALYSIS COMPLETED IN $runtime SECONDS\n"
printf "### EXITING\n"
exit 0

# stockholm to fasta
#fgrep -e "#" -e "//" -v IGHV_RSS3_7bat.sto | awk '{gsub(/-/, "", $2)} {print $1, "\t", $2}' | awk ' {print ">"$1"\n"$2 '} >IGHV_RSS3_7bat.fna
#fgrep -e "#" -e "//" -v IGHJ_7bat.sto | awk '{gsub(/-/, "", $2)} {print $1, "\t", $2}' | awk ' {print ">"$1"\n"$2 '} >IGHJ_7bat.fna

# Usar P_discolor porque esta en un solo scaffold

# Para generar tablas de conteo de tamaño de clusters
# grep "^>" cdhit_IGHV_RSS3.clstr >cdhit_IGHV_RSS3.clstr.clusters
# awk '!/^>/{count++}/^>/{print count; count = 0} END {print count}' cdhit_IGHV_RSS3.clstr | sed '/^$/d' >cdhit_IGHV_RSS3.clstr.count
# paste cdhit_IGHV_RSS3.clstr.clusters cdhit_IGHV_RSS3.clstr.count > cdhit_IGHV_RSS3.clstr.tbl
# rm cdhit_IGHV_RSS3.clstr.clusters cdhit_IGHV_RSS3.clstr.count
