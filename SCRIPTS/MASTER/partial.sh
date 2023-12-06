#!/bin/bash

SPECIE="DESMODUS_ROTUNDUS"
SHORT_GS="Dero"

#echo ./RESULTS/$SPECIE/EXONERATE 
#echo $SHORT_GS 
#echo ./RESULTS/$SPECIE/EXONERATE/exonerate_IGHV_cDNA_vs_$SHORT_GS

cat <<-'EOF' > /dev/null 2>&1



Rscript ./SCRIPTS/SUBSCRIPTS/predict_ighd_by_rss.R \
-f ./RESULTS/$SPECIE/GFF_RESULTS/nhmmer_RSS_IGHD_5_${SHORT_GS}.gff \
-t ./RESULTS/$SPECIE/GFF_RESULTS/nhmmer_RSS_IGHD_3_${SHORT_GS}.gff >/dev/null 2>&1

mv RSS_D_IGH_D_segments.gff ./RESULTS/$SPECIE/D_SEGMENTS/PREDICTION/D_segments_and_RSS_${SHORT_GS}.gff

mkdir -p ./RESULTS/$SPECIE/CONDENSED_RESULTS
EOF

Rscript ./SCRIPTS/SUBSCRIPTS/locate_nearby_rss.R \
-n ./RESULTS/$SPECIE/HMMER/RSS_IGHV/nhmmer_RSS_IGHV_${SHORT_GS}.gff \
-t ./RESULTS/$SPECIE/HMMER/RSS_IGHV/nhmmer_RSS_IGHV.tbl \
-v ./RESULTS/$SPECIE/OVERLAP_REDUCTION/overlap_gene_prediction_IGHV_vs_${SHORT_GS}.gff \
-r 50 \
-m V_segments >/dev/null 2>&1

mv v_rss_analysis.gff ./RESULTS/$SPECIE/V_RSS_CORRECTED_GFF/IGHV_RSS_analysis_${SHORT_GS}.gff


cat ./RESULTS/$SPECIE/V_RSS_CORRECTED_GFF/IGHV_RSS_analysis_${SHORT_GS}.gff \
./RESULTS/$SPECIE/J_RSS_CORRECTED_GFF/IGHJ_RSS_analysis_${SHORT_GS}.gff \
./RESULTS/$SPECIE/D_SEGMENTS/PREDICTION/D_segments_and_RSS_${SHORT_GS}.gff \
./RESULTS/$SPECIE/OVERLAP_REDUCTION/overlap_gene_prediction_IGHV_cDNA_vs_${SHORT_GS}.gff \
./RESULTS/$SPECIE/OVERLAP_REDUCTION/overlap_exon_prediction_IGHV_cDNA_vs_${SHORT_GS}.gff \
./RESULTS/$SPECIE/EXONERATE/exonerate_IGHC_cDNA_vs_${SHORT_GS}.gff | \
grep -v "#" >./RESULTS/$SPECIE/CONDENSED_RESULTS/${SHORT_GS}_cond.gff








#Rscript ./SCRIPTS/SUBSCRIPTS/blast_m6_to_gff.R --file ./RESULTS/$SPECIE/BLAST/tblastn_IGHV_aa_vs_$SHORT_GS.m6 --source tblastn --bitscore 0


# mv ighv_rss_corrected_coordinates.gff ./RESULTS/$SPECIE/V_RSS_CORRECTED_GFF/IGHV_RSS_merged.gff