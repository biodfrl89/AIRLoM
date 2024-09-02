# Name: blast_m6_to_gff.R
# Description: This scripts reformat the m6 format of the blast output into a gff formated file 
# Author: David Felipe Rendón Luna 
# Date: December-2023
# Update: August-2024 by Ernestina Godoy (elizabeth.godoy@insp.mx)

# CHECK LIBRARIES AND ARGUMENTS -------------------------------------------

# Check for optparse library to load arguments from command line
if(suppressMessages(!require("optparse"))) {
  stop("optparse was not found. Exiting.")
}

# Load parser
library("optparse")
suppressMessages(library("dplyr"))

# Create list of arguments and help asociadted to receive.
opt_list = list(make_option(opt_str = c("-f", "--file"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "The m6 blast file that is going to be formated.", 
                            metavar = "[FILENAME]"),
                make_option(opt_str = c("-b", "--bitscore"), 
                            action="store",
                            type = "numeric", 
                            default = NULL, 
                            help = "The score used to filter the records.", 
                            metavar = "[INTEGER]"),
                make_option(opt_str = c("-s", "--source"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "The name of the program used to obtain the records.", 
                            metavar = "[STRING]")
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

if (is.null(opt$file)) {
  print_help(opt_parser)
  stop("A nhmmer file must be submitted", call.=FALSE)
}

if (is.null(opt$bitscore)) {
  print_help(opt_parser)
  stop("A tbl file must be submitted", call.=FALSE)
}

if (is.null(opt$source)) {
  print_help(opt_parser)
  stop("An overlap file must be submitted", call.=FALSE)
}

# ASSIGN ARGUMENTS -----------------------------------------------------

# TEST
#file <- "./DATA/tblastx_IGHV_cDNA_vs_Phdi.m6"
#bitscore <- 0
#method <- "tblastx"
#outname <- sub(".m6", ".gff", file)

file <- opt$file
bitscore <- opt$bitscore
source <- opt$source
outname <- sub(".m6", ".gff", file)

# SCRIPT

data <- read.delim(file, header = FALSE, 
                   col.names = c("query_seqid", "subject_seqid", "perc_ident", "length", "mismatch", "gap_open", "query_start", "query_end", "subject_start", "subject_end", "e_value", "bit_score"))

# Sort by subject_seqid, subject_start and bit_score (decrease)
data_sort <- data[order(data$subject_seqid,data$subject_start,-data$bit_score, data$e_value),]

# remove redundancy 
data_sort$ss_se <- paste(data_sort$subject_start, data_sort$subject_end, sep=".")

data_sort_uniq <- data_sort %>%
  distinct(ss_se, .keep_all = TRUE)

# Build GFF
df <- data.frame(seqid = data_sort_uniq$subject_seqid,
                 source = source,
                 type = "gene",
                 start = data_sort_uniq$subject_start,
                 end = data_sort_uniq$subject_end,
                 score = ".", strand = ".", phase = ".",
                 attributes = paste0("Query=", data_sort_uniq$query_seqid, ";evalue=", data_sort_uniq$e_value, ";bitscore=", data_sort_uniq$bit_score, ";pident=", data_sort_uniq$perc_ident))

print("Start correcting coordinates and set strand value")

positives <- which(df$end>df$start)

if (length(positives)>0){
  df$strand[positives] <- "+"
  df$strand[-positives] <- "-"
  rec_start <- df$end[-positives]
  rec_end <- df$start[-positives]
  df$start[-positives] <- rec_start
  df$end[-positives] <- rec_end
} else{
  df$strand <- "-"
  rec_start <- df$end
  rec_end <- df$start
  df$start <- rec_start
  df$end <- rec_end
}

write.table(df, outname, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
