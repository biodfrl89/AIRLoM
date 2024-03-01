# Name: blast_m6_to_gff.R
# Description: This scripts reformat the m6 format of the blast output into a gff formated file 
# Author: David Felipe Rendón Luna 
# Date: December-2023

# CHECK LIBRARIES AND ARGUMENTS -------------------------------------------

# Check for optparse library to load arguments from command line
if(suppressMessages(!require("optparse"))) {
  stop("optparse was not found. Exiting.")
}

# Load parser
library("optparse")
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

df <- data.frame(seqid = NA, source = NA, type = NA, start = NA, end = NA, score = NA, strand = NA, phase = NA, attributes = NA) # Begin empty df with gff3 columns

for (i in 1:nrow(data)) {
  if (data[i,]$bit_score > bitscore ) {
  # Make the record
    df[i,] <- c(data[i,]$subject_seqid,
               source,
               "gene",
               data[i,]$subject_start,
               data[i,]$subject_end, 
               ".", ".", ".",
               paste0("Query=", data[i,]$query_seqid))
  }
  # Add positive strand symbol
  if (df[i,]$end > df[i,]$end) {
    df[i,]$strand <- "+"
  } else {
  # Add negative strand symbol and correct coordinates 
    rec_start <- df[i,]$end
    rec_end <- df[i,]$start
    df[i,]$start <- rec_start
    df[i,]$end <- rec_end
    df[i,]$strand <- "-"
  }
}

write.table(df, outname, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)