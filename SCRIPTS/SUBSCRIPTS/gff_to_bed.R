# Name: gff_to_bed.R
# Description: 
# Author: David Felipe Rendón Luna 
# Date: February-2024


# Check for optparse library to load arguments from command line ----------
if(suppressMessages(!require("optparse"))) {
  stop("optparse was not found. Exiting.")
}

# Load parser -------------------------------------------------------------
library("optparse")
# Create list of arguments and help asociadted to receive.
opt_list = list(make_option(opt_str = c("-f", "--file"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "GFF files to reduce the ambiguity", 
                            metavar = "[FILENAME]"))

# Make the parser
opt_parser = OptionParser(option_list = opt_list)

# Load the arguments in the parser into an object
opt = parse_args(opt_parser)

# Check for arguments -----------------------------------------------------
if (length(opt) == 1) {
  print_help(opt_parser)
  stop("No arguments submitted")
}

if (is.null(opt$file)) {
  print_help(opt_parser)
  stop("A file must be submitted", call.=FALSE)
}

# SCRIPT
library(utils)

# Read arguments
FILE <- opt$file
OUTFILE <- gsub(".gff$", ".bed", FILE)

#Extract specie name
obj_grep <- regexpr(pattern = "_.{4}.gff",text = FILE)
name_specie <- substr(FILE, obj_grep + 1, obj_grep + attr(obj_grep, "match.length") - 5)

# Read file
df <- read.delim(FILE, header = FALSE, col.names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "atribute"))
gff <- rtracklayer::import(FILE, "gff")

# Split df attributes by semicolons
splited <- strsplit(df$atribute, ";")

# Create function for cleaning attributes
check_ID <- function(x) {
  if(length(x[grepl("ID=", x)]) >= 1) {
    return(x[grepl("ID=", x)]) }
  else {
    return(x)
  }
}

# Replace attributes
df$atribute <- as.vector(do.call(rbind, lapply(splited, check_ID)))

# Create BED
bed_df <- data.frame(chrom = df$seqname,
                     chromStart = df$start,
                     chromEnd = df$end,
                     name = df$atribute,
                     score = 0, 
                     strand = df$strand,
                     thickStart = df$start,
                     thickEnd = df$end)
bed_df$itemRgb <- with(bed_df,ifelse(strand == "+", "255,0,0", 
                                     ifelse(strand == "-", "0,0,255", "0,255,0" )))

# Select only RSS and exons
bed_df <- bed_df[grepl("RSS|_exon_", bed_df$name),]

# Eliminate "ID=" from name
bed_df$name <- gsub("ID=", "", bed_df$name)

# Make second field
second_field <- c()
c <- 1
for (i in seq_along(bed_df$name)) {
  if (grepl("RSS", bed_df$name[c])) {
    second_field[c] <- "IGHV_RSS3"  
  } else if (grepl("IGHV_exon", bed_df$name[c])) {
    second_field[c] <- "IGHV_exon"  
  } 
  c <- c + 1
}

# Make final attribute for bed column "name"
bed_df$name <- paste0(name_specie, "_", gsub("_plus|_minus", "", bed_df$name), "|", 
                      second_field, "|", 
                      "|",
                      with(bed_df,ifelse(strand == "+", "FOR", ifelse(strand == "-", "REV", "." ))))

#### Correct coordinates for plus and minus
# Select strands
df_plus <- bed_df[bed_df$strand == "+",]
df_minus <- bed_df[bed_df$strand == "-",]

# Correct coordinates
df_plus$chromStart <- df_plus$chromStart - 1
df_plus$thickStart <- df_plus$thickStart - 1
df_minus$chromStart <- df_minus$chromStart - 1
df_minus$thickStart <- df_minus$thickStart - 1

# Make bed again
bed_df <- rbind(df_plus, df_minus)

# Separate genes 
genes <- gff[grepl("IGHV_gene_.*_", gff$ID)]
features <- gff[grepl("RSS_V-associated_.*|IGHV_exon_.*_", gff$ID)]

# Extend genes ranges
genes@ranges@start <- as.integer(genes@ranges@start - 100)
genes@ranges@width <- as.integer(genes@ranges@width + 200)

#Make overlaps
gr <- GenomicRanges::findOverlaps(genes, features)

# Return ranges to original
genes@ranges@start <- as.integer(genes@ranges@start + 100)
genes@ranges@width <- as.integer(genes@ranges@width - 200)

vec_conversion <- c()
for (i in seq_along(gr)) {
  vec_conversion[gr@to[i]] <- gr@from[i]
}

# Make final attribute for bed column "name"
bed_df$name <- paste0(paste0(name_specie, "_IGHV_", sprintf('%0.3d', vec_conversion)), "|", 
                      second_field, "|", 
                      "|",
                      with(bed_df,ifelse(strand == "+", "FOR", ifelse(strand == "-", "REV", "." ))), "|")


bed_df$seg_size <- abs(bed_df$chromStart - bed_df$chromEnd) 

for (i in 1:nrow(bed_df)) {
  if (bed_df$seg_size[i] < 60 & grepl("IGHV_exon", bed_df$name[i])) {
    second_field[i] <- "SP_exon"
  }
}

# Make final attribute for bed column "name"
bed_df$name <- paste0(paste0(name_specie, "_IGHV_", sprintf('%0.3d', vec_conversion)), "|", 
                      second_field, "|", 
                      "|",
                      with(bed_df,ifelse(strand == "+", "FOR", ifelse(strand == "-", "REV", "." ))), "|")

bed_df$seg_size <- NULL

# Save results
write.table(bed_df, file = OUTFILE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

