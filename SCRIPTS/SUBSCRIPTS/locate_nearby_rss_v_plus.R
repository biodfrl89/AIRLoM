# Name: locate_nearby_rss_v_plus.R
# Description: This script is designed to detect which rss are detected in the vecinity of a exon or gene sequence.
# This is achieved by increasing the range of exon or gene both at the 5' and 3' ends of V segments, and then trying to make an overlap,
# with the detected RSS sequences from nhmmer fro V segments.
# The coordinates of overlapping V segment and RSS are mutualy corrected.
# Author: David Felipe Rendón Luna 
# Date: November-2023

# CHECK LIBRARIES AND ARGUMENTS -------------------------------------------

# Check for optparse library to load arguments from command line
if(suppressMessages(!require("optparse"))) {
  stop("optparse was not found. Exiting.")
}

# Check for aditional libraries
if(nzchar(system.file(package = "rtracklayer"))) {
  cat("-rtracklayer library found.\n")
} else {
  stop("rtracklayer was not found. Exiting.\n")
}

if(nzchar(system.file(package = "GenomicRanges"))) {
  cat("-GenomicRanges library found.\n")
} else {
  stop("GenomicRanges was not found. Exiting.\n")
}

# Load parser
library("optparse")
# Create list of arguments and help asociadted to receive.
opt_list = list(make_option(opt_str = c("-n", "--nhmmer"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "nhmmer gff file corresponding to the founded RSS of the desired region (V-D-J)", 
                            metavar = "[FILENAME]"),
                make_option(opt_str = c("-t", "--tbl"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "nhmmer tbl file corresponding to the founded RSS of the desired region (V-D-J)", 
                            metavar = "[FILENAME]"),
                make_option(opt_str = c("-v", "--overlap"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "Overlap or reduced file of exonerate gene file from the desired region (V-D-J)", 
                            metavar = "[FILENAME]"),
                make_option(opt_str = c("-r", "--range"), 
                            action="store",
                            type = "numeric", 
                            default = NULL, 
                            help = "The range of extension of the V segment in order to detect any RSS signal", 
                            metavar = "[INT]")
                )

# Make the parser
opt_parser = OptionParser(option_list = opt_list)

# Load the arguments in the parser into an object
opt = parse_args(opt_parser)

# CHECK ARGUMENTS -----------------------------------------------------
if (length(opt) == 1) {
  print_help(opt_parser)
  stop("No arguments submitted")
}

if (is.null(opt$nhmmer)) {
  print_help(opt_parser)
  stop("A nhmmer file must be submitted", call.=FALSE)
}

if (is.null(opt$tbl)) {
  print_help(opt_parser)
  stop("A tbl file must be submitted", call.=FALSE)
}

if (is.null(opt$overlap)) {
  print_help(opt_parser)
  stop("An overlap file must be submitted", call.=FALSE)
}

if (is.null(opt$range)) {
  print_help(opt_parser)
  stop("A range must be submitted", call.=FALSE)
}

# SCRIPT -----------------------------------------------------

# Load files from overlap and overlaped exonerate results
gff_nhmmer <- rtracklayer::import.gff(opt$nhmmer)
gff_overlaps <- rtracklayer::import.gff(opt$overlap)
gff_overlaps_raw <- gff_overlaps
hmm_tbl <- read.table(opt$tbl, 
                        col.names = c("target_name", "accession", "query_name", "accession", "hmm_from", "hmm_to", "ali_from", "ali_to", "enf_from", "env_to", 
                                      "modlen", "strand", "E-value", "score", "bias", "description_target"))
detection_range <- opt$range

# Increase the range of the exons 
gff_overlaps@ranges@start <- as.integer(gff_overlaps@ranges@start - detection_range)
gff_overlaps@ranges@width <- as.integer(gff_overlaps@ranges@width + detection_range * 2)

# Make overlaps in both strands
overlaps <- GenomicRanges::findOverlaps(query = gff_nhmmer, subject = gff_overlaps)



# Checkpoint for overlaps and make analysis in plus strand
if (length(overlaps_plus) < 1) {
  print("No overlaps detected for plus strand")
} else {
  # Copy original df of subject to make the selection
  gff_overlaps_plus <- gff_overlaps_raw
  gff_overlaps_plus <- gff_overlaps_plus[which(as.character(gff_overlaps_plus@strand) == "+")]
    
  # Make a loop to cycle trough all the located querys
  for (i in seq_along(overlaps_plus)) {
    # Select from and to using position index
    j_query <- overlaps_plus@from[i]
    j_subject <- overlaps_plus@to[i]
      
    # Calculate RSS real begining for the record in nhhmer query
    rss_real_begin <- hmm_tbl_plus[j_query, "ali_from"] + 1 - hmm_tbl_plus[j_query, "hmm_from"] 
      
    # Replace RSS real begining at the start of RSS
    gff_nhmmer_plus@ranges@start[j_query] <- as.integer(rss_real_begin)
      
    # Calculate the real RSS end
    #rss_real_end <- rss_real_begin + hmm_tbl[j_query,"hmm_to"] + (39 - hmm_tbl[j_query,"hmm_to"]) - 1
    rss_real_end <- rss_real_begin + hmm_tbl[j_query, "modlen"] 
      
    # Replace width for gff_nhmmer
    gff_nhmmer_plus@ranges@width[j_query] <- as.integer(rss_real_end - rss_real_begin) 
      
    # Replace RSS real begining at the end of V segment
    v_start <- gff_overlaps_plus[j_subject]@ranges@start
      
    v_new_width <- abs(rss_real_begin - v_start)
      
    gff_overlaps_plus[j_subject]@ranges@width <- as.integer(v_new_width)
  }

  # Initialice column of metadata for RSS
  gff_nhmmer_plus$ID <- NA
  # Fill metadata
  if (length(gff_nhmmer_plus[overlaps_plus@from]) >= 1) gff_nhmmer_plus[overlaps_plus@from]$ID <- paste0("RSS_V-associated_", seq_along(gff_nhmmer_plus[overlaps_plus@from]))
  if (length(gff_nhmmer_plus[which(is.na(gff_nhmmer_plus$ID))]) >= 1) gff_nhmmer_plus[which(is.na(gff_nhmmer_plus$ID))]$ID <- paste("RSS_V-non-associated_", seq_along(gff_nhmmer_plus[which(is.na(gff_nhmmer_plus$ID))]))
    
  # Select V segments with nearby RSS signals, with corrected coordinates from plus strand
  final_plus <- c(gff_overlaps_plus[overlaps_plus@to], gff_nhmmer_plus[overlaps_plus@from])