#' process.cell
#'
#' @param bam path to a bam file
#' @param bedpe path to the accompanying bedfile
#' @param est.ploidy estimated ploidy of the cell (from estimate.ploidy)
#' @param bin.size number of nucleotides per bin
#' @param min.svSize number of bins for the smallest confident structural variant
#'
#' @return
#' @export
#'
#' @examples
process.cell <- function(bam, genome, bedpe = NULL, bin.size = 500000, min.svSize = 1e6, min_length = 50, max_length = 1000, tag_overlap = 9, ext_correction = NULL){
  min.svSize <- min.svSize/bin.size

  reads <- load_cell(bam, binSize = bin.size, genome)
  reads.cor <- Songbird::convert_long(reads)
  # Lets try normalizing the read depth?
  reads.cor$reads <- (reads.cor$reads/sum(reads.cor$reads, na.rm = T))*10000
  num_reads <- c(sum(reads.cor$uncorrected.reads))

  # Get mode CN indices
  reads.cor$ubh_tx <- ubh_segment(reads.cor$reads, reads.cor$use, min_svSize = min.svSize)
  reads.cor$num_reads <- num_reads
  reads.cor$bam_file <- bam
  reads.cor$bedpe_file <- bedpe

  if(is.null(bedpe)){
    reads.cor$est_ploidy <- NA
  }else{
    out <- estimate.ploidy(sample = bedpe, genome = genome, binSize = bin.size, min_length = min_length, max_length = max_length, tag_overlap = tag_overlap, use_external = ext_correction)
    reads.cor$ratio <- out$ratio
    reads.cor$est_ploidy <- out$est_ploidy
    reads.cor$breadth <- out$breadth
    reads.cor$coverage <- out$coverage
    reads.cor$prop_doublet_tags <- out$prop_doublet_tags
    reads.cor$avg_length <- out$avg_length
    reads.cor$overlap_genome_size <- out$overlap_genome_size
  }
  return(reads.cor)
}


#' #' Title
#' #'
#' #' @param ploidy estimated ploidy
#' #' @param prop_doublets proportion of doublets
#' #' @param coverage coverage of the cell
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' correct_ratio <- function(ratio, prop_doublets, use_external = NULL){
#'   # Right now its just a place holder function
#'   if(is.null(use_external)){
#'     corr_data <- Songbird::correction_data
#'   }
#'   else{
#'     corr_data <- read.table(use_external, header = T, sep = '\t')
#'   }
#'
#'   # Process the correction table with a known ploidy
#'   corr_data$true_ratio <- (corr_data$ploidy-1)/corr_data$ploidy
#'   corr_data$correction <- corr_data$est_ratio/corr_data$true_ratio
#'   cor_function <- stats::lm(correction ~ doublet_prop, data = corr_data)
#'
#'   # Predict the correction factors for each ratio using the linear function
#'   correction_factor <- stats::predict(cor_function, newdata = data.frame(doublet_prop = prop_doublets))
#'   return(ratio/correction_factor)
#' }

#' convert_long
#'
#' @param reads reads object from QDNASeq
#'
#' @return table of reads for each cell
#'
#' @examples
convert_long <- function(reads){
  uncorReads <- reads$reads
  reads <- reads$reads.cor
  metadata <- reads@featureData@data

  # Add the corrected reads
  table <- reads@assayData$copynumber
  out_table <- c()
  for(i in 1:ncol(table)){
    new_dat <- metadata
    new_dat$reads <- table[,i]
    new_dat$cell_id <- colnames(table)[i]
    out_table <- rbind(out_table, new_dat)
  }

  # Add the raw reads
  table <- uncorReads@assayData$counts
  out_table$uncorrected.reads <- table[,1]
  return(out_table)
}


#' load_cell
#'
#' @param bamPath path to the bam file
#' @param binSize size of the bins in # of nucleotides
#'
#' @return
#'
#' @examples
load_cell <- function(bamPath, binSize, genome){

  bins <- QDNAseq::getBinAnnotations(binSize = binSize/1000, genome = genome)
  bins@data$mappability <- as.numeric(bins@data$mappability)*100
  reads <- QDNAseq::binReadCounts(bins, bamfiles = bamPath, pairedEnds = T)
  #reads <- QDNAseq::applyFilters(reads, residual = F, blacklist = T)
  reads <- QDNAseq::estimateCorrection(reads)
  reads.cor <- QDNAseq::correctBins(reads)
  return(list(reads = reads, reads.cor = reads.cor))
}

#' wh_transform
#'
#' @param vals values for the walsh haddamayer transform
#' @param sv_binSize minimum confident structural variant size
#'
#' @return
#'
#' @examples
ubh_segment <- function(values, use, min_svSize){

  # Strip NAs and save their position
  out <- rep(NA, length(values))
  real_idxs <- which(!is.na(values))
  values <- values[!is.na(values)]
  sigma <- calc_madOffset(values)

  # Run UBH
  transform <- unbalhaar::best.unbal.haar(values)
  transform.filt <- unbalhaar::hard.thresh(transform, sigma)
  #transform.filt <- prune_offsets(transform.filt, values, min_svSize)
  reconstruct <- unbalhaar::reconstr(transform.filt)
  out[real_idxs] <- reconstruct
  return(out)
}

#' prune_offsets
#'
#' @param ubh_obj path to the bedpe file produced by the preprocessing pipeline
#' @param reads bin data from the QDNASeq object
#' @param min_svSize minimum read length (default 30 nt)
#'
#' @return a data frame set up as a traditional bed file
#'
#' @examples

prune_offsets <- function(ubh_obj, reads, min_svSize){
  tree <- ubh_obj$tree
  split_data <- c()
  for (i in 1:length(tree)) {
    entry <- tree[[i]]
    s <- entry[5,] - entry[3,]
    left_bias <- (entry[4,] - entry[3,])/s
    right_bias <- (entry[5,] - entry[4,])/s
    max_bias <- apply(cbind(left_bias, right_bias), 1, max)

    dat <- data.frame(split_bias = max_bias, span = s,
                      tree_level = i, column = seq_along(max_bias),
                      coeff = entry[2,])
    split_data <- rbind(split_data, dat)
  }

  split_data$threshold <- (split_data$span-0.05*split_data$span)/split_data$span
  split_data$offset_entry <- (split_data$split_bias >= split_data$threshold) & (split_data$span > min_svSize)

  offset_data <- split_data[split_data$offset_entry, ]
  offset_data$ubh_vector_peak <- NA
  for(i in 1:nrow(offset_data)){
    start <- 1
    split <- (1-offset_data$split_bias[i]) * offset_data$span[i] + 1
    end <- offset_data$span[i]
    offset_data$ubh_vector_peak[i] <- max(unbalhaar::unbal.haar.vector(c(start, split, end)))
  }

  thresh <- sd(diff(reads))
  offset_data$coeff[offset_data$coeff < -thresh] <- -thresh
  offset_data$coeff[offset_data$coeff > thresh] <- thresh
  offset_data$coeff <- 0

  for(i in 1:nrow(offset_data)){
    tree[[offset_data$tree_level[i]]][2, offset_data$column[i]] <- offset_data$coeff[i]
  }

  ubh_obj$tree <- tree
  return(ubh_obj)
}

#' load.preprocess.bed
#'
#' @param bedpe_file path to the bedpe file produced by the preprocessing pipeline
#' @param bin_data bin data from the QDNASeq object
#' @param min_length minimum read length (default 30 nt)
#' @param max_length maximum read length (default 1000 nt)
#'
#' @return a data frame set up as a traditional bed file
#'
#' @examples
load.preprocess.bed <- function(bedpe_file, bin_data, min_length = 30, max_length = 1000){
  bedpe <- utils::read.table(bedpe_file, sep = '\t')
  colnames(bedpe) <- c('Chr1', 'Start1', 'End1', 'Chr2', 'Start2', 'End2', 'Name', 'Score', 'R1_direction', 'R2_direction')
  bed <- data.frame(Chr = bedpe$Chr1, Start = bedpe$Start1, End = bedpe$End1,
                    Name = paste0(bedpe$Name, ':', bedpe$R1_direction, '/', bedpe$R2_direction))


  # Set the read start as the minimum of start1 or start2, and read end as maximum of end1 and end2
  bed$Start <- apply(cbind(bedpe$Start1, bedpe$Start2), 1, min)
  bed$End <- apply(cbind(bedpe$End1, bedpe$End2), 1, max)

  bed$Length <- bed$End - bed$Start

  bed <- bed[(bed$Length > min_length) & (bed$Length < max_length),]
  bed$Strandedness <- extractStrandedness(bed$Name, 3)

  # Grab the first three characters of the chromosome names to see if its prefixed properly. If not add chr
  chr_prefix <- substr(bed$Chr, 1, 3)
  if(sum(grepl('chr', chr_prefix)) == 0){
    bed$Chr <- paste0('chr', bed$Chr)
  }

  bed <- bed[grep('(chr[0-9]+|X|Y)$', bed$Chr),]

  bed <- filter.bed(bed, bin_data)
  return(bed)
}

#' filter.bed
#'
#' @param bed a data frame formatted as a bed file with read positions
#' @param bin_data the metadata generated by QDNAseq for the bins
#'
#' @return a data frame with duplicate reads (identical start or end sites) removed
#'
#' @examples
#'
filter.bed <- function(bed, bin_data) {
  # Remove duplicated reads that aren't caught by the standard dedup tools
  bed <- bed[with(bed, order(Chr, Start, End, Strandedness)),]
  tmp.start <- paste(bed$Chr, bed$Start, bed$Strandedness)
  tmp.end <- paste(bed$Chr, bed$End, bed$Strandedness)

  are.duplicated <- duplicated(tmp.start, fromLast = TRUE) | duplicated(tmp.end)
  bed <- bed[!are.duplicated,]

  # Remove reads that are in bins with low mappability
  hq_bins <- GenomicRanges::GRanges(paste0('chr', bin_data$chromosome),
                                    IRanges::IRanges(bin_data$start, bin_data$end))
  bed_ranges <- GenomicRanges::GRanges(bed$Chr, IRanges::IRanges(bed$Start, bed$End))
  #bed <- bed[bed_ranges %over% hq_bins,]
  bed <- bed[GenomicRanges::countOverlaps(bed_ranges, hq_bins)>0,]
  return(bed)
}

#' extractStrandedness
#'
#' @param readName name for each read pulled from the bed file
#' @param n number of characters to count back to get strandedness
#'
#' @return read orientation from the last n characters of the read name
#'
#' @examples
extractStrandedness <- function(readName, n){
  readName.count <- nchar(readName)
  return(substr(readName, readName.count - n + 1, readName.count))
}

#' count.doublets
#'
#' @param bed a data frame formatted as a bed file with read positions must be filtered by filter.bed
#' @param min.tag.overlap smallest length of the tagmentation read overlap (default 9 nt)
#' @param max.tag.overlap largest length of the tagmentation read overlap (default 10 nt)
#'
#' @return
#'
#' @examples
count.doublets <- function(bed, min.tag.overlap = 9, max.tag.overlap = 10){
  read.starts <- GenomicRanges::GRanges(seqnames = bed$Chr,
                         ranges = IRanges::IRanges(start = bed$Start+min.tag.overlap, end = bed$Start+max.tag.overlap))
  read.ends <- GenomicRanges::GRanges(seqnames = bed$Chr,
                       ranges = IRanges::IRanges(start = bed$End, end = bed$End))

  num.doublets <- sum(GenomicRanges::countOverlaps(read.starts, read.ends)>0)
  return(num.doublets/nrow(bed))
}

#' count.overlaps
#'
#' @param bed a data frame formatted as a bed file with read positions must be filtered by filter.bed
#' @param min.size minimum read length (default 50 nt)
#' @param max.size maximum read length (default 1000 nt)
#' @param tag.overlap largest length of the tagmentation read overlap (default 10 nt)
#'
#' @return the bed file with the number of reads overlapping the upstream and over regions for each read
#' @export
#'
#' @examples
count.overlaps <- function(bed, min.size = 50, max.size = 1000, tag.overlap = 10) {
  # bed <- tmp2
  upstream.ranges <- IRanges::IRanges(start = bed$Start - max.size + tag.overlap,
                                      end = bed$Start - 1)
  upstream.counting <- GenomicRanges::GRanges(seqnames = bed$Chr,
                               ranges = upstream.ranges)

  overlap.ranges <- IRanges::IRanges(start = bed$Start - min.size + tag.overlap,
                                      end = bed$Start - 1)
  over.counting <- GenomicRanges::GRanges(seqnames = bed$Chr,
                           ranges = overlap.ranges)

  bed.sized <- bed[bed$Length > min.size & bed$Length < max.size, ]
  counting.ranges <- IRanges::IRanges(start = bed.sized$Start,
                                      end = bed.sized$Start)
  for.counting <- GenomicRanges::GRanges(seqnames = bed.sized$Chr,
                          ranges = counting.ranges)

  counted.over.overlaps <- GenomicRanges::countOverlaps(over.counting, for.counting)
  counted.upstream.overlaps <- GenomicRanges::countOverlaps(upstream.counting, for.counting)

  data.out <- data.frame(Count.Upstream = counted.upstream.overlaps,
                         Count.Over = counted.over.overlaps,
                         Norm.Count.Upstream = counted.upstream.overlaps / (max.size - tag.overlap),
                         Norm.Count.Over = counted.over.overlaps / (min.size - tag.overlap))
  data.out <- cbind(bed, data.out)
  return(data.out)
}

#' calc.breadth
#'
#' @param bed imported bed data frame from bedpe object
#'
#' @return the breadth of the bed file
#'
#' @examples
calc.breadth <- function(bed){
  read_ranges <- GenomicRanges::GRanges(seqnames = bed$Chr,
                                        ranges = IRanges::IRanges(start = bed$Start, end = bed$End))
  coverage <- GenomicRanges::coverage(read_ranges)

  breadth <- 0
  for(i in 1:length(coverage)){
    breadth <- breadth + sum(coverage[[i]]@lengths[coverage[[i]]@values>0])
  }
  return(breadth)
}

#' Title
#'
#' @param sample path to bedpe file
#' @param binSize binning size in nucleotides
#' @param min_length minimum read length (default 50 nt)
#' @param max_length maximum read length (default 1000 nt)
#'
#' @return
#' @export
#'
#' @examples
estimate.ploidy <- function(sample, binSize, genome, min_length = 50, max_length = 1000, tag_overlap = 9, use_external = NULL){
  if(genome == 'hg19'){
    stop("hg19 is not supported. Please align your data to chm13v2 (best performance) or hg38 (good performance)")
  }

  # Create the bin aggregated data to identify regions where we want to measure the overlap statistics
  bins <- QDNAseq::getBinAnnotations(binSize/1000, genome = genome)
  bin_data <- bins@data

  if(genome == 'chm13v2'){
    bin_data$mappability <- as.numeric(bin_data$mappability) * 100
  }

  bin_data$use[bin_data$mappability < 90 | bin_data$bases < 90] <- FALSE
  bin_data <- bin_data[bin_data$use,]
  bin_data$binName <- paste0('chr', bin_data$chromosome, '_', bin_data$start)

  # Load bed and get QC and Ploidy Estimation Metrics
  bed <- load.preprocess.bed(sample, bin_data, min_length, max_length)
  if(nrow(bed)==0){
    out <- data.frame(ratio = NA, breadth = NA, coverage = NA, prop_doublet_tags = NA,
                      avg_length = NA, overlap_genome_size = NA, ploidy_readCount = NA)
    return(out)
  }
  prop_doublets <- count.doublets(bed, min.tag.overlap = tag_overlap, max.tag.overlap = tag_overlap+1)
  bed <- count.overlaps(bed, min.size = 100, max.size = max_length, tag.overlap = tag_overlap)

  # Merge the bed data with the bin data
  bed$binStart <- floor(bed$Start/binSize)*binSize+1
  bed$binName <- paste0(bed$Chr, '_', bed$binStart)
  bed_aggregate <- stats::aggregate(bed[,c("Count.Upstream",
                                           "Norm.Count.Upstream",
                                           "Count.Over",
                                           "Norm.Count.Over")],
                             by = list(bed$binName), function(x) sum(x, na.rm = T))
  bed_counts <- stats::aggregate(bed$binName, by = list(bed$binName), length)
  bed_coverage <- stats::aggregate(bed$Length, by = list(bed$binName), function(x) sum(x, na.rm = T))

  match_idx <- match(bin_data$binName, bed_aggregate$Group.1)
  bin_data <- cbind(bin_data, bed_aggregate[match_idx,])
  bin_data$reads <- bed_counts$x[match_idx]

  bin_data$coverage <- bed_coverage$x[match_idx]
  bin_data$prop_doublets <- prop_doublets
  bin_data <- bin_data[,!grepl('Group', colnames(bin_data))]

  out <- data.frame(ratio = mean(bin_data$Norm.Count.Over, na.rm = T)/mean(bin_data$Norm.Count.Upstream, na.rm = T),
                    breadth = calc.breadth(bed),
                    coverage = sum(bin_data$coverage, na.rm = T),
                    prop_doublet_tags = prop_doublets,
                    avg_length = mean(bed$Length, na.rm = T),
                    overlap_genome_size = nrow(bin_data)*binSize,
                    ploidy_readCount = nrow(bed))

  #if(genome == 'hg38'){
  #  #out$corrected_ratio <- correct_ratio(out$ratio, out$prop_doublet_tags, use_external = use_external)
  #  out$corrected_ratio <- out$ratio
  #}else{
  #  out$corrected_ratio <- out$ratio
  #}
  #out$est_ploidy <- 1/(1-out$corrected_ratio)
  out$est_ploidy <- 1/(1-out$ratio)
  return(out)
}
