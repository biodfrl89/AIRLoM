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

# Make final attribute for bed column "name"
bed_df$name <- paste0(name_specie, "_", gsub("_plus|_minus", "", bed_df$name), "|", 
                      "|", 
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

# Save results
write.table(bed_df, file = OUTFILE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

