#!/bin/bash

############# USAGE #################
usage () {
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
exit 1
}

############# CHECK ARGUMENTS #################
check_arguments () {
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
}

############# CDHIT #################

run_cdhit () {
    printf "### Running CD-HIT to cluster raw database sequences redundancy\n"
    cd-hit-est -i $IGH_V_CDNA -o cdhit_IGHV_CDS.fna -c 0.95 -n 5 -M 2000 -T 4 >/dev/null 2>&1
    cd-hit -i $IGH_V_AA -o cdhit_IGHV_CDS.faa -c 0.95 -n 5 -M 2000 -T 4 >/dev/null 2>&1
    cd-hit-est -i $IGH_J_CDNA -o cdhit_IGHJ.fna -c 0.95 -n 5 -M 2000 -T 4 >/dev/null 2>&1
    cd-hit-est -i $IGH_C_CDNA -o cdhit_IGHC.fna -c 0.95 -n 5 -M 2000 -T 4 >/dev/null 2>&1
    cd-hit-est -i $RSS_V -o cdhit_IGHV_RSS3.fna -c 0.95 -n 5 -M 2000 -T 4 >/dev/null 2>&1
    cd-hit-est -i $RSS_D5 -o cdhit_IGHD_RSS5.fna -c 0.95 -n 5 -M 2000 -T 4 >/dev/null 2>&1
    cd-hit-est -i $RSS_D3 -o cdhit_IGHD_RSS3.fna -c 0.95 -n 5 -M 2000 -T 4 >/dev/null 2>&1
    cd-hit-est -i $RSS_J -o cdhit_IGHJ_RSS5.fna -c 0.95 -n 5 -M 2000 -T 4 >/dev/null 2>&1

    printf "### Running CLUSTAL alignment in RSS clustered sequences\n"
    clustalw2 -align -quiet -infile=cdhit_IGHV_RSS3.fna -OUTFILE=clustal_IGHV_RSS3.fna -OUTPUT=FASTA >/dev/null 2>&1
    clustalw2 -align -quiet -infile=cdhit_IGHD_RSS5.fna -OUTFILE=clustal_IGHD_RSS5.fna -OUTPUT=FASTA >/dev/null 2>&1
    clustalw2 -align -quiet -infile=cdhit_IGHD_RSS3.fna -OUTFILE=clustal_IGHD_RSS3.fna -OUTPUT=FASTA >/dev/null 2>&1
    clustalw2 -align -quiet -infile=cdhit_IGHJ_RSS5.fna -OUTFILE=clustal_IGHJ_RSS5.fna -OUTPUT=FASTA >/dev/null 2>&1

    printf "### Moving database files\n"
    rm cdhit*.dnd
    rm cdhit*RSS?.fna
    mv clustal*.fna ./SCRIPTS/SEQ_DB/RSS
    mv cdhit*.fna cdhit*.faa ./SCRIPTS/SEQ_DB/FASTAS
    mv cdhit*.clstr ./SCRIPTS/SEQ_DB/FASTAS/CLUSTERS
}

############# REASSINGMENT #################

reassign_vars () {
    printf "### Reassigning variables\n"
    IGH_V_CDNA="./SCRIPTS/SEQ_DB/FASTAS/cdhit_IGHV_CDS.fna"
    IGH_V_AA="./SCRIPTS/SEQ_DB/FASTAS/cdhit_IGHV_CDS.faa"
    IGH_J_CDNA="./SCRIPTS/SEQ_DB/FASTAS/cdhit_IGHJ.fna"
    IGH_C_CDNA="./SCRIPTS/SEQ_DB/FASTAS/cdhit_IGHC.fna"
    RSS_V="./SCRIPTS/SEQ_DB/RSS/clustal_IGHV_RSS3.fna"
    RSS_D5="./SCRIPTS/SEQ_DB/RSS/clustal_IGHD_RSS5.fna"
    RSS_D3="./SCRIPTS/SEQ_DB/RSS/clustal_IGHD_RSS3.fna"
    RSS_J="./SCRIPTS/SEQ_DB/RSS/clustal_IGHJ_RSS5.fna"
}

############# CREATE SPECIES SHORT NAME #################
format_name () {
    printf "### Creating short specie name\n"
    SPECIE_L=${SPECIE,,} #To lower all
    SHORT_GENUS=${SPECIE_L:0:2}
    SPECIE_NAME=$(echo $SPECIE_L | cut -d "_" -f2)
    SHORT_SPECIE=${SPECIE_NAME:0:2}
    SHORT_GS=${SHORT_GENUS^}${SHORT_SPECIE} #Capitalize genus.
}

############# MOVE GENOME AND MAKE BLAST DATABASE #################

# Create specie directory with reference genome directory. Move genome to reference genome folder and locate inside.
move_genome () {
    mkdir -p ./RESULTS/$SPECIE/REFERENCE_GENOME && \
    cp $GENOME ./RESULTS/$SPECIE/REFERENCE_GENOME/ && \
    cd ./RESULTS/$SPECIE/REFERENCE_GENOME/
}

make_blast_database () {
    # Decompress genome and pipe it to makeblastdb
    printf "### Making BLAST Database\n"

    gunzip -dc $GENOME | makeblastdb -in - -dbtype nucl -out BLAST_DB -title $GENOME >/dev/null 2>&1

    # Move behind and make blast directory and database directory.
    cd .. && mkdir -p BLAST/DB 

    # Move to DB directory and move all BLAST file to that location.
    cd ./BLAST/DB/ && mv ../../REFERENCE_GENOME/BLAST* .

    # Return to original path.
    cd ../../../..
}

############# BLAST FUNCTIONS #################

blast_v_cdna () {
    # V region with cDNA fasta
    printf "### Running BLAST for V regions using cDNA\n"
    tblastx \
    -query $IGH_V_CDNA \
    -db ./RESULTS/$SPECIE/BLAST/DB/BLAST_DB \
    -outfmt 6 \
    -evalue 1e-20 \
    -num_threads 12 \
    -max_target_seqs 100 \
    -out ./RESULTS/$SPECIE/BLAST/tblastx_IGHV_cDNA_vs_$SHORT_GS.m6
}

blast_v_aa () {
    # V region with aa fasta
    printf "### Running BLAST for V regions using AA\n"
    tblastn \
    -query $IGH_V_AA \
    -db ./RESULTS/$SPECIE/BLAST/DB/BLAST_DB \
    -outfmt 6 \
    -evalue 1e-20 \
    -num_threads 12 \
    -max_target_seqs 100 \
    -out ./RESULTS/$SPECIE/BLAST/tblastn_IGHV_aa_vs_$SHORT_GS.m6
}

blast_j_cdna () {
    # J region with cDNA fasta
    printf "### Running BLAST for J regions using cDNA\n"
    tblastx \
    -query $IGH_J_CDNA \
    -db ./RESULTS/$SPECIE/BLAST/DB/BLAST_DB \
    -outfmt 6 \
    -evalue 1e-4 \
    -num_threads 12 \
    -max_target_seqs 100 \
    -out ./RESULTS/$SPECIE/BLAST/tblastx_IGHJ_cDNA_vs_$SHORT_GS.m6
}

blast_c_cdna () {
    # C region with cDNA fasta
    printf "### Running BLAST for C regions using cDNA\n"
    tblastx \
    -query $IGH_C_CDNA \
    -db ./RESULTS/$SPECIE/BLAST/DB/BLAST_DB \
    -outfmt 6 \
    -evalue 1e-20 \
    -num_threads 12 \
    -max_target_seqs 100 \
    -out ./RESULTS/$SPECIE/BLAST/tblastx_IGHC_cDNA_vs_$SHORT_GS.m6
}

############# REFORMAT M6 TO GFF #################

reformat_blast_v_cdna () {
    printf "### Reformating blast m6 tables of IGHV_cDNA to gff\n"
    #conda activate R
    Rscript ./SCRIPTS/SUBSCRIPTS/blast_m6_to_gff.R --file ./RESULTS/$SPECIE/BLAST/tblastx_IGHV_cDNA_vs_$SHORT_GS.m6 --source tblastx --bitscore 0
    rm ./RESULTS/$SPECIE/BLAST/tblastx_IGHV_cDNA_vs_$SHORT_GS.m6
}

reformat_blast_v_aa () {
    printf "### Reformating blast m6 tables of IGHV_aa to gff\n"
    Rscript ./SCRIPTS/SUBSCRIPTS/blast_m6_to_gff.R --file ./RESULTS/$SPECIE/BLAST/tblastn_IGHV_aa_vs_$SHORT_GS.m6 --source tblastn --bitscore 0
    rm ./RESULTS/$SPECIE/BLAST/tblastn_IGHV_aa_vs_$SHORT_GS.m6
}

reformat_blast_j_cdna () {
    printf "### Reformating blast m6 tables of IGHJ_cDNA to gff\n"
    Rscript ./SCRIPTS/SUBSCRIPTS/blast_m6_to_gff.R --file ./RESULTS/$SPECIE/BLAST/tblastx_IGHJ_cDNA_vs_$SHORT_GS.m6 --source tblastx --bitscore 0
    rm ./RESULTS/$SPECIE/BLAST/tblastx_IGHJ_cDNA_vs_$SHORT_GS.m6
}

reformat_blast_c_cdna (){
    printf "### Reformating blast m6 tables of IGHC_cDNA to gff\n"
    Rscript ./SCRIPTS/SUBSCRIPTS/blast_m6_to_gff.R --file ./RESULTS/$SPECIE/BLAST/tblastx_IGHC_cDNA_vs_$SHORT_GS.m6 --source tblastx --bitscore 0
    rm ./RESULTS/$SPECIE/BLAST/tblastx_IGHC_cDNA_vs_$SHORT_GS.m6
}

############# EXTRACT SCAFFOLDS NAME #################
extract_scaffolds () {
    printf "### Getting scaffolds ID and extract their sequences\n"
    cut -f1 ./RESULTS/$SPECIE/BLAST/tblastx_IGHV_cDNA_vs_$SHORT_GS.gff | sort | uniq >./RESULTS/$SPECIE/BLAST/scaffolds_tblastx_IGHV_cDNA_vs_$SHORT_GS
    cut -f1 ./RESULTS/$SPECIE/BLAST/tblastn_IGHV_aa_vs_$SHORT_GS.gff | sort | uniq >./RESULTS/$SPECIE/BLAST/scaffolds_tblastn_IGHV_vs_$SHORT_GS
    cut -f1 ./RESULTS/$SPECIE/BLAST/tblastx_IGHJ_cDNA_vs_$SHORT_GS.gff | sort | uniq >./RESULTS/$SPECIE/BLAST/scaffolds_tblastx_IGHJ_cDNA_vs_$SHORT_GS
    cut -f1 ./RESULTS/$SPECIE/BLAST/tblastx_IGHC_cDNA_vs_$SHORT_GS.gff | sort | uniq >./RESULTS/$SPECIE/BLAST/scaffolds_tblastx_IGHC_cDNA_vs_$SHORT_GS
    #cat ./RESULTS/$SPECIE/BLAST/scaffolds_t* | sort | uniq -c | awk ' $1 >= 3 ' | tr -s ' ' | cut -d' ' -f3 >./RESULTS/$SPECIE/BLAST/scaffolds_to_extract_count
    cat ./RESULTS/$SPECIE/BLAST/scaffolds_t* | sort | uniq -c | tr -s ' ' | cut -d' ' -f3 >./RESULTS/$SPECIE/BLAST/scaffolds_to_extract
    rm ./RESULTS/$SPECIE/BLAST/scaffolds_tblast*

}


############# EXTRACT SCAFFOLDS SEQUENCE #################
extract_scaffolds_seq () {
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
}

############# EXONERATE FUNCTIONS #################
exonerate_v_cdna () {
    printf "### Running EXONERATE for V regions using cDNA\n"
    exonerate \
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
}

exonerate_v_aa () {
    printf "### Running EXONERATE for V regions using AA \n"
    exonerate \
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
}

exonerate_j_cdna () {
    printf "### Running EXONERATE for J regions using cDNA\n"
    # Add score after manual inspection
    exonerate \
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
}

exonerate_c_cdna () {
    # C region with cDNA fasta
    printf "### Running EXONERATE for C regions using cDNA\n"

    # Change intron size of C region
    exonerate \
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
}

############# EXTRACT VULGAR FROM EXONERATE #################
extract_exonerate_vulgar () {
    printf "### Extracting vulgar anotation from exonerate output file\n"
    grep "vulgar" ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_cDNA_vs_$SHORT_GS | sed 's/vulgar: //' >./RESULTS/$SPECIE/EXONERATE/vulgar_IGHV_cDNA_vs_$SHORT_GS
    grep "vulgar" ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_vs_$SHORT_GS | sed 's/vulgar: //' >./RESULTS/$SPECIE/EXONERATE/vulgar_IGHV_vs_$SHORT_GS
    grep "vulgar" ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHC_cDNA_vs_$SHORT_GS | sed 's/vulgar: //' >./RESULTS/$SPECIE/EXONERATE/vulgar_IGHC_cDNA_vs_$SHORT_GS
    grep "vulgar" ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHJ_cDNA_vs_$SHORT_GS | sed 's/vulgar: //' >./RESULTS/$SPECIE/EXONERATE/vulgar_IGHJ_cDNA_vs_$SHORT_GS
}

############# EXTRACT GFF FROM EXONERATE #################
extract_exonerate_gff () {
    printf "### Making exonerate gff output file\n"
    grep "#" -v  ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_cDNA_vs_$SHORT_GS | grep "vulgar" -v >./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_cDNA_vs_${SHORT_GS}.gff
    grep "#" -v  ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_vs_$SHORT_GS | grep "vulgar" -v >./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_vs_${SHORT_GS}.gff
    grep "#" -v  ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHC_cDNA_vs_$SHORT_GS | grep "vulgar" -v >./RESULTS/$SPECIE/EXONERATE/exonerate_IGHC_cDNA_vs_${SHORT_GS}.gff
    grep "#" -v  ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHJ_cDNA_vs_$SHORT_GS | grep "vulgar" -v >./RESULTS/$SPECIE/EXONERATE/exonerate_IGHJ_cDNA_vs_${SHORT_GS}.gff
}

############# CONVERT VULGAR TO TABLE #################

vulgar_to_table () {
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
}

############# HMMER #################
run_hmmer () {
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
        python3 ./SCRIPTS/SUBSCRIPTS/hmmer_tbl_to_gff.py --file ./RESULTS/$SPECIE/HMMER/$dir/nhmmer_$dir.tbl --source hmmerscan --bitscore 0 

        # Rename output file
        mv ./RESULTS/$SPECIE/HMMER/$dir/nhmmer_$dir.gff ./RESULTS/$SPECIE/HMMER/$dir/nhmmer_${dir}_${SHORT_GS}.gff 
    done
}

############# EXONERATE FILTRATION BY EXONS AND GENES, PLUS AND MINUS #################
exonerate_filtration () {
    ############# GET PLUS EXONS IN EXONERATE #################
    awk ' $3 ~ "exon" && $7 ~ "+" ' ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_cDNA_vs_${SHORT_GS}.gff >./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHV_cDNA_vs_${SHORT_GS}_plus_exons.gff
    awk ' $3 ~ "exon" && $7 ~ "+" ' ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_vs_${SHORT_GS}.gff >./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHV_vs_${SHORT_GS}_plus_exons.gff 
    awk ' $3 ~ "exon" && $7 ~ "+" ' ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHC_cDNA_vs_${SHORT_GS}.gff >./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHC_cDNA_vs_${SHORT_GS}_plus_exons.gff
    awk ' $3 ~ "exon" && $7 ~ "+" ' ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHJ_cDNA_vs_${SHORT_GS}.gff >./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHJ_cDNA_vs_${SHORT_GS}_plus_exons.gff

    ############# GET MINUS EXONS IN EXONERATE #################
    awk ' $3 ~ "exon" && $7 ~ "-" ' ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_cDNA_vs_${SHORT_GS}.gff >./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHV_cDNA_vs_${SHORT_GS}_minus_exons.gff
    awk ' $3 ~ "exon" && $7 ~ "-" ' ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_vs_${SHORT_GS}.gff >./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHV_vs_${SHORT_GS}_minus_exons.gff 
    awk ' $3 ~ "exon" && $7 ~ "-" ' ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHC_cDNA_vs_${SHORT_GS}.gff >./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHC_cDNA_vs_${SHORT_GS}_minus_exons.gff
    awk ' $3 ~ "exon" && $7 ~ "-" ' ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHJ_cDNA_vs_${SHORT_GS}.gff >./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHJ_cDNA_vs_${SHORT_GS}_minus_exons.gff

    ############# GET PLUS GENES IN EXONERATE #################
    awk ' $3 ~ "gene" && $7 ~ "+" ' ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_cDNA_vs_${SHORT_GS}.gff >./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHV_cDNA_vs_${SHORT_GS}_plus_genes.gff
    awk ' $3 ~ "gene" && $7 ~ "+" ' ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_vs_${SHORT_GS}.gff >./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHV_vs_${SHORT_GS}_plus_genes.gff 
    awk ' $3 ~ "gene" && $7 ~ "+" ' ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHC_cDNA_vs_${SHORT_GS}.gff >./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHC_cDNA_vs_${SHORT_GS}_plus_genes.gff
    awk ' $3 ~ "gene" && $7 ~ "+" ' ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHJ_cDNA_vs_${SHORT_GS}.gff >./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHJ_cDNA_vs_${SHORT_GS}_plus_genes.gff

    ############# GET MINUS GENES IN EXONERATE #################
    awk ' $3 ~ "gene" && $7 ~ "-" ' ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_cDNA_vs_${SHORT_GS}.gff >./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHV_cDNA_vs_${SHORT_GS}_minus_genes.gff
    awk ' $3 ~ "gene" && $7 ~ "-" ' ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_vs_${SHORT_GS}.gff >./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHV_vs_${SHORT_GS}_minus_genes.gff 
    awk ' $3 ~ "gene" && $7 ~ "-" ' ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHC_cDNA_vs_${SHORT_GS}.gff >./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHC_cDNA_vs_${SHORT_GS}_minus_genes.gff
    awk ' $3 ~ "gene" && $7 ~ "-" ' ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHJ_cDNA_vs_${SHORT_GS}.gff >./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHJ_cDNA_vs_${SHORT_GS}_minus_genes.gff
}

############# REDUCE EXONERATE FILTERED FILES #################
exonerate_reduction () {
    printf "### Runing Rscript for overlap reduction over exonerate plus exons \n"
    Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHV_cDNA_vs_${SHORT_GS}_plus_exons.gff >/dev/null 2>&1
    Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHV_vs_${SHORT_GS}_plus_exons.gff >/dev/null 2>&1
    Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHC_cDNA_vs_${SHORT_GS}_plus_exons.gff >/dev/null 2>&1
    Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHJ_cDNA_vs_${SHORT_GS}_plus_exons.gff >/dev/null 2>&1

    printf "### Runing Rscript for overlap reduction over exonerate minus exons \n"
    Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHV_cDNA_vs_${SHORT_GS}_minus_exons.gff >/dev/null 2>&1
    Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHV_vs_${SHORT_GS}_minus_exons.gff >/dev/null 2>&1
    Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHC_cDNA_vs_${SHORT_GS}_minus_exons.gff >/dev/null 2>&1
    Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHJ_cDNA_vs_${SHORT_GS}_minus_exons.gff >/dev/null 2>&1

    printf "### Runing Rscript for overlap reduction over exonerate plus genes \n"
    Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHV_cDNA_vs_${SHORT_GS}_plus_genes.gff >/dev/null 2>&1
    Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHV_vs_${SHORT_GS}_plus_genes.gff >/dev/null 2>&1
    Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHC_cDNA_vs_${SHORT_GS}_plus_genes.gff >/dev/null 2>&1
    Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHJ_cDNA_vs_${SHORT_GS}_plus_genes.gff >/dev/null 2>&1

    printf "### Runing Rscript for overlap reduction over exonerate minus genes \n"
    Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHV_cDNA_vs_${SHORT_GS}_minus_genes.gff >/dev/null 2>&1
    Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHV_vs_${SHORT_GS}_minus_genes.gff >/dev/null 2>&1
    Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHC_cDNA_vs_${SHORT_GS}_minus_genes.gff >/dev/null 2>&1
    Rscript ./SCRIPTS/SUBSCRIPTS/gff_disambiguation.R --file ./RESULTS/$SPECIE/EXONERATE_FILTERED/exonerate_IGHJ_cDNA_vs_${SHORT_GS}_minus_genes.gff >/dev/null 2>&1

    mv reduced_exonerate* ./RESULTS/$SPECIE/REDUCTION
}

############# OVERLAP ANALYSIS ON V SEGMENTS USING CDNA #################
exonerate_overlap_v_cdna () {
    printf "### Runing overlap analysis for IGHV using cDNA - plus \n"
    Rscript ./SCRIPTS/SUBSCRIPTS/predict_ighv_by_overlaps.R \
    --query ./RESULTS/$SPECIE/REDUCTION/reduced_exonerate_IGHV_cDNA_vs_${SHORT_GS}_plus_exons.gff \
    --subject ./RESULTS/$SPECIE/REDUCTION/reduced_exonerate_IGHV_cDNA_vs_${SHORT_GS}_plus_genes.gff \
    --strand plus >/dev/null 2>&1

    mv ./gene_prediction.gff ./RESULTS/$SPECIE/OVERLAP/overlap_gene_prediction_IGHV_cDNA_vs_${SHORT_GS}_plus.gff 
    mv ./exon_prediction.gff ./RESULTS/$SPECIE/OVERLAP/overlap_exon_prediction_IGHV_cDNA_vs_${SHORT_GS}_plus.gff 

    printf "### Runing overlap analysis for IGHV using cDNA - minus \n"
    Rscript ./SCRIPTS/SUBSCRIPTS/predict_ighv_by_overlaps.R \
    --query ./RESULTS/$SPECIE/REDUCTION/reduced_exonerate_IGHV_cDNA_vs_${SHORT_GS}_minus_exons.gff \
    --subject ./RESULTS/$SPECIE/REDUCTION/reduced_exonerate_IGHV_cDNA_vs_${SHORT_GS}_minus_genes.gff \
    --strand minus >/dev/null 2>&1

    mv ./gene_prediction.gff ./RESULTS/$SPECIE/OVERLAP/overlap_gene_prediction_IGHV_cDNA_vs_${SHORT_GS}_minus.gff 
    mv ./exon_prediction.gff ./RESULTS/$SPECIE/OVERLAP/overlap_exon_prediction_IGHV_cDNA_vs_${SHORT_GS}_minus.gff 
}

############# OVERLAP ANALYSIS ON V SEGMENTS USING AA #################
exonerate_overlap_v_aa () {
    printf "### Runing overlap analysis for IGHV using aa - plus \n"
    Rscript ./SCRIPTS/SUBSCRIPTS/predict_ighv_by_overlaps.R \
    --query ./RESULTS/$SPECIE/REDUCTION/reduced_exonerate_IGHV_vs_${SHORT_GS}_plus_exons.gff \
    --subject ./RESULTS/$SPECIE/REDUCTION/reduced_exonerate_IGHV_vs_${SHORT_GS}_plus_genes.gff \
    --strand plus >/dev/null 2>&1

    mv ./gene_prediction.gff ./RESULTS/$SPECIE/OVERLAP/overlap_gene_prediction_IGHV_vs_${SHORT_GS}_plus.gff 
    mv ./exon_prediction.gff ./RESULTS/$SPECIE/OVERLAP/overlap_exon_prediction_IGHV_vs_${SHORT_GS}_plus.gff 

    printf "### Runing overlap analysis for IGHV using aa - minus \n"
    Rscript ./SCRIPTS/SUBSCRIPTS/predict_ighv_by_overlaps.R \
    --query ./RESULTS/$SPECIE/REDUCTION/reduced_exonerate_IGHV_vs_${SHORT_GS}_minus_exons.gff \
    --subject ./RESULTS/$SPECIE/REDUCTION/reduced_exonerate_IGHV_vs_${SHORT_GS}_minus_genes.gff \
    --strand minus >/dev/null 2>&1

    mv ./gene_prediction.gff ./RESULTS/$SPECIE/OVERLAP/overlap_gene_prediction_IGHV_vs_${SHORT_GS}_minus.gff 
    mv ./exon_prediction.gff ./RESULTS/$SPECIE/OVERLAP/overlap_exon_prediction_IGHV_vs_${SHORT_GS}_minus.gff 
}

############# CORRECT COORDINATES FOR J AND RSS#################
correct_j_plus () {
    printf "### Making overlap search between V segments and V RSS, and correct joining coordinates. Plus strand.\n"
    Rscript ./SCRIPTS/SUBSCRIPTS/locate_nearby_rss_j_plus.R \
    -n ./RESULTS/$SPECIE/HMMER/RSS_IGHJ/nhmmer_RSS_IGHJ_${SHORT_GS}.gff \
    -t ./RESULTS/$SPECIE/HMMER/RSS_IGHJ/nhmmer_RSS_IGHJ.tbl \
    -v ./RESULTS/$SPECIE/REDUCTION/reduced_exonerate_IGHJ_cDNA_vs_${SHORT_GS}_plus_genes.gff \
    -r 30 >/dev/null 2>&1

    test -f "j_rss_plus_analysis.gff" && mv j_rss_plus_analysis.gff ./RESULTS/$SPECIE/J_RSS_CORRECTED/J_RSS_plus_analysis_${SHORT_GS}.gff
}

correct_j_minus () {
    printf "### Making overlap search between V segments and V RSS, and correct joining coordinates. Minus strand.\n"
    Rscript ./SCRIPTS/SUBSCRIPTS/locate_nearby_rss_j_minus.R \
    -n ./RESULTS/$SPECIE/HMMER/RSS_IGHJ/nhmmer_RSS_IGHJ_${SHORT_GS}.gff \
    -t ./RESULTS/$SPECIE/HMMER/RSS_IGHJ/nhmmer_RSS_IGHJ.tbl \
    -v ./RESULTS/$SPECIE/REDUCTION/reduced_exonerate_IGHJ_cDNA_vs_${SHORT_GS}_minus_genes.gff \
    -r 30 >/dev/null 2>&1

    test -f "j_rss_minus_analysis.gff" && mv j_rss_minus_analysis.gff ./RESULTS/$SPECIE/J_RSS_CORRECTED/J_RSS_minus_analysis_${SHORT_GS}.gff
}

############# DETECT D SEGMENTS #################

detect_d () {
    printf "### Predict D segments using the D-RSS sequences.\n"

    Rscript ./SCRIPTS/SUBSCRIPTS/predict_ighd_by_rss.R \
    -f ./RESULTS/$SPECIE/HMMER/RSS_IGHD_5/nhmmer_RSS_IGHD_5_${SHORT_GS}.gff \
    -t ./RESULTS/$SPECIE/HMMER/RSS_IGHD_3/nhmmer_RSS_IGHD_3_${SHORT_GS}.gff >/dev/null 2>&1

    mv RSS_D_IGH_D_segments.gff ./RESULTS/$SPECIE/D_SEGMENTS/D_RSS_analysis_${SHORT_GS}.gff
}

############# CORRECT COORDINATES FOR V AND RSS #################
correct_v_plus () {
    printf "### Making overlap search between V segments and V RSS, and correct joining coordinates. Plus strand.\n"
    Rscript ./SCRIPTS/SUBSCRIPTS/locate_nearby_rss_v_plus.R \
    -n ./RESULTS/$SPECIE/HMMER/RSS_IGHV/nhmmer_RSS_IGHV_${SHORT_GS}.gff \
    -t ./RESULTS/$SPECIE/HMMER/RSS_IGHV/nhmmer_RSS_IGHV.tbl \
    -v ./RESULTS/$SPECIE/OVERLAP/overlap_gene_prediction_IGHV_vs_${SHORT_GS}_plus.gff \
    -r 30 >/dev/null 2>&1

    test -f "v_rss_plus_analysis.gff" && mv v_rss_plus_analysis.gff ./RESULTS/$SPECIE/V_RSS_CORRECTED/V_RSS_plus_analysis_${SHORT_GS}.gff
}

correct_v_minus () {
    printf "### Making overlap search between V segments and V RSS, and correct joining coordinates. Minus strand.\n"
    Rscript ./SCRIPTS/SUBSCRIPTS/locate_nearby_rss_v_minus.R \
    -n ./RESULTS/$SPECIE/HMMER/RSS_IGHV/nhmmer_RSS_IGHV_${SHORT_GS}.gff \
    -t ./RESULTS/$SPECIE/HMMER/RSS_IGHV/nhmmer_RSS_IGHV.tbl \
    -v ./RESULTS/$SPECIE/OVERLAP/overlap_gene_prediction_IGHV_vs_${SHORT_GS}_minus.gff \
    -r 30 >/dev/null 2>&1

    test -f "v_rss_minus_analysis.gff" && mv v_rss_minus_analysis.gff ./RESULTS/$SPECIE/V_RSS_CORRECTED/V_RSS_minus_analysis_${SHORT_GS}.gff
}

############# MERGE GFFS #################

merge_gffs () {
    # Check for directory existence
    if [ ! -d "RESULTS/$SPECIE/GFF" ]; then
        mkdir RESULTS/$SPECIE/GFF
    fi

    # Create empty final temp gff file
    touch RESULTS/$SPECIE/GFF/temp.gff

    printf "### Merging results.\n"

    [[ -f RESULTS/$SPECIE/D_SEGMENTS/D_RSS_analysis_${SHORT_GS}.gff ]] && cat RESULTS/$SPECIE/D_SEGMENTS/D_RSS_analysis_${SHORT_GS}.gff >>temp.gff
    [[ -f RESULTS/$SPECIE/V_RSS_CORRECTED/V_RSS_plus_analysis_${SHORT_GS}.gff ]] && cat RESULTS/$SPECIE/V_RSS_CORRECTED/V_RSS_plus_analysis_${SHORT_GS}.gff >>temp.gff
    [[ -f RESULTS/$SPECIE/V_RSS_CORRECTED/V_RSS_minus_analysis_${SHORT_GS}.gff ]] && cat RESULTS/$SPECIE/V_RSS_CORRECTED/V_RSS_minus_analysis_${SHORT_GS}.gff >>temp.gff
    [[ -f RESULTS/$SPECIE/J_RSS_CORRECTED/J_RSSq_plus_analysis_${SHORT_GS}.gff ]] && cat RESULTS/$SPECIE/J_RSS_CORRECTED/J_RSS_plus_analysis_${SHORT_GS}.gff >>temp.gff
    [[ -f RESULTS/$SPECIE/J_RSS_CORRECTED/J_RSS_minus_analysis_${SHORT_GS}.gff ]] && cat RESULTS/$SPECIE/J_RSS_CORRECTED/J_RSS_minus_analysis_${SHORT_GS}.gff >>temp.gff
    [[ -f RESULTS/$SPECIE/REDUCTION/reduced_exonerate_IGHC_cDNA_vs_${SHORT_GS}_minus_genes.gff ]] && cat RESULTS/$SPECIE/REDUCTION/reduced_exonerate_IGHC_cDNA_vs_${SHORT_GS}_minus_genes.gff >>temp.gff
    [[ -f RESULTS/$SPECIE/REDUCTION/reduced_exonerate_IGHC_cDNA_vs_${SHORT_GS}_minus_exons.gff ]] && cat RESULTS/$SPECIE/REDUCTION/reduced_exonerate_IGHC_cDNA_vs_${SHORT_GS}_minus_exons.gff >>temp.gff
    [[ -f RESULTS/$SPECIE/REDUCTION/reduced_exonerate_IGHC_cDNA_vs_${SHORT_GS}_plus_genes.gff ]] && cat RESULTS/$SPECIE/REDUCTION/reduced_exonerate_IGHC_cDNA_vs_${SHORT_GS}_plus_genes.gff >>temp.gff
    [[ -f RESULTS/$SPECIE/REDUCTION/reduced_exonerate_IGHC_cDNA_vs_${SHORT_GS}_plus_exons.gff ]] && cat RESULTS/$SPECIE/REDUCTION/reduced_exonerate_IGHC_cDNA_vs_${SHORT_GS}_plus_exons.gff >>temp.gff
    [[ -f RESULTS/$SPECIE/OVERLAP/overlap_exon_prediction_IGHV_vs_${SHORT_GS}_plus.gff ]] && cat RESULTS/$SPECIE/OVERLAP/overlap_exon_prediction_IGHV_vs_${SHORT_GS}_plus.gff >>temp.gff
    [[ -f RESULTS/$SPECIE/OVERLAP/overlap_exon_prediction_IGHV_vs_${SHORT_GS}_minus.gff ]] && cat RESULTS/$SPECIE/OVERLAP/overlap_exon_prediction_IGHV_vs_${SHORT_GS}_minus.gff >>temp.gff

    # Clean gff file and give final name
    grep -v "#" temp.gff >RESULTS/$SPECIE/GFF/final_gff_${SHORT_GS}.gff

    # Remove temp
    rm RESULTS/$SPECIE/GFF/temp.gff
    rm temp.gff
}

############# MAKE BED FROM FINAL BED #################
gff_to_bed () {
    printf "### Making BED file from final GFF\n"
    Rscript ./SCRIPTS/SUBSCRIPTS/gff_to_bed.R \
    -f ./RESULTS/$SPECIE/GFF/final_gff_${SHORT_GS}.gff

    echo 'track name="Bed_results" description="Bed anotation for manual inspection" visibility=2 itemRgb="On"' | cat - RESULTS/$SPECIE/GFF/final_gff_${SHORT_GS}.bed >RESULTS/$SPECIE/GFF/final_bed_${SHORT_GS}.bed
    rm RESULTS/$SPECIE/GFF/final_gff_${SHORT_GS}.bed
}

