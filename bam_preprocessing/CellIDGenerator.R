# Read command line args
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2){stop("Incorrect number of arguments check script for info")}

inFile <- args[1] # First argument has to be the input bam
cellIDFile <- args[2] # Second argument is the full path to cell IDs

options(scipen=999)
date <- Sys.Date()

# Get header from all cells
header <- system(paste('samtools view -H', inFile), intern = T, ignore.stderr = T)

# Create Cell IDs for splitting the bam file
cellIDs <- grep('^@CO', header, value = T)
cellIDs <- gsub('.*CB:(.*)', '\\1', cellIDs)
#cellIDs <- data.frame(cellIDs, well = gsub('.*(R[0-9]+\\-C[0-9]+)', '\\1', cellIDs))
cellIDs <- data.frame(cellIDs, well = cellIDs)

# write out the cell IDs
write.table(x = cellIDs, file = cellIDFile, sep = '\t', quote = F, col.names = F, row.names = F)