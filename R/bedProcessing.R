#' process.batch
#'
#' @param bams path to a bam file
#' @param genome genome assembly to use (e.g. hg38, hg19)
#' @param bedpes path to the accompanying bedfile
#' @param bin.size number of nucleotides per bin
#' @param min_length minimum length for reads
#' @param max_length maximum length for reads
#' @param tag_overlap size of the read overlap created by a tagmentation event
#' @param n_cpu number of CPUs to use for parallel processing, default is NULL which uses all but one CPU
#'
#' @return songbird object
#' @export
#'
#' @examples

process.batch <- function(bams, genome = 'hg38', bedpes = NULL, bin.size = 500000, min_length = 50, max_length = NULL, tag_overlap = 9, n_cpu=NULL){
  if(is.null(n_cpu)){
    n_cpu <- parallel::detectCores() - 1
  }

  if(is.null(bedpes)){
    res <- pbmcapply::pbmclapply(1:length(bams), function(i) process.cell(bams[i], genome = genome, bin.size = bin.size, min_length = min_length, tag_overlap = tag_overlap), mc.cores = n_cpu)
  }else{
    res <- pbmcapply::pbmclapply(1:length(bams), function(i) process.cell(bams[i], genome = genome, bedpe = bedpes[i], bin.size = bin.size, min_length = min_length, max_length = max_length, tag_overlap = tag_overlap), mc.cores = n_cpu)
    #res <- lapply(1:length(bams), function(i) process.cell(bams[i], genome = genome, bedpe = bedpes[i], bin.size = bin.size, min_length = min_length, tag_overlap = tag_overlap))
  }
  return(create_sce(res))
}

#' process.cell
#'
#' @param bam path to a bam file
#' @param genome genome assembly to use (e.g. hg38, hg19)
#' @param bedpe path to the accompanying bedfile
#' @param bin.size number of nucleotides per bin
#' @param min_length minimum length for reads
#' @param max_length maximum length for reads
#' @param tag_overlap size of the read overlap created by a tagmentation event
#'
#' @return corrected reads
process.cell <- function(bam, genome, bedpe = NULL, bin.size = 500000, min_length = 50, max_length = NULL, tag_overlap = 9){
  if(!genome %in% c('hg38', 'hg19', 'chm13v2')){
    stop('Genome must be one of hg38, hg19, or chm13v2')
  }

  reads <- load_cell(bam, binSize = bin.size, genome)
  reads.cor <- convert_long(reads)
  # Lets try normalizing the read depth?
  reads.cor$reads <- (reads.cor$reads/sum(reads.cor$reads, na.rm = T))*10000
  num_reads <- c()

  # Get mode CN indices
  if(genome == 'chm13v2'){
    reads.cor$mappability <- reads.cor$mappability * 100
  }
  reads.cor$ubh_tx <- ubh_segment(reads.cor$reads)
  reads.cor$bin.depth <- reads.cor$reads/reads.cor$ubh_tx
  reads.cor$bam_file <- bam
  reads.cor$bedpe_file <- bedpe

  if(is.null(bedpe)){
    reads.cor$est_ploidy <- NA
  }else{
    reads.cor <- estimate.ploidy(sample = bedpe, bin_data = reads.cor, min_length = min_length, max_length = NULL, tag_overlap = tag_overlap)
  }
  return(reads.cor)
}

#' convert_long
#'
#' @param reads reads object from QDNASeq
#'
#' @return table of reads for each cell
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
#' @param genome genome to use for binning
#'
#' @return
#'
#' @examples
load_cell <- function(bamPath, binSize, genome){
  bins <- QDNAseq::getBinAnnotations(binSize = binSize/1000, genome = genome, verbose = F)
  bins@data$mappability <- as.numeric(bins@data$mappability)
  reads <- QDNAseq::binReadCounts(bins, bamfiles = bamPath, pairedEnds = T, verbose = F)
  reads <- QDNAseq::applyFilters(reads, verbose = F)
  reads <- QDNAseq::estimateCorrection(reads)
  reads <- QDNAseq::applyFilters(reads, chromosomes = NA, verbose = F)
  reads.cor <- QDNAseq::correctBins(reads)
  return(list(reads = reads, reads.cor = reads.cor))
}

#' ubh_segment
#'
#' @param values GC + Map corrected reads per bin
#' @param use bins to use
#'
#' @return
#'
#' @examples
ubh_segment <- function(values){

  # Strip NAs and save their position
  out <- rep(NA, length(values))
  real_idxs <- which(!is.na(values))
  filt_values <- values[!is.na(values)]
  init_sigma <- calc_madOffset(filt_values)

  # Run UBH
  transform <- unbalhaar::best.unbal.haar(filt_values)
  transform.filt <- unbalhaar::hard.thresh(transform, init_sigma)
  reconstr <- unbalhaar::reconstr(transform.filt)
  reconstr[(filt_values == 0) | (reconstr < 0)] = 0
  out[real_idxs] <- reconstr

  # Bridge the exclusion regions by setting NAs to their next value
  na_idx <- which(is.na(values))
  next_idx <- na_idx + 1
  next_idx[next_idx > length(out)] <- length(out)-1
  out[na_idx] <- out[next_idx]

  # Singlets which are far higher than the main distribution are likely true values, but make sure that ubh hasn't mis-estimated them
  # These are probs from focal amplifications where evidence is low, but CN fidelity is important
  bin_var <- var(out[(out > quantile(out, 0.05, na.rm = T)) & (out < quantile(out, 0.95, na.rm = T))])
  shift_idx <- which(c(0, diff(out)) != 0)
  singlets <- shift_idx[out[shift_idx] != out[shift_idx + 1]]

  # Make sure that high singlets actually reflect the underlying value - UBH can overestimate them
  peak_singlets <- singlets[(abs(out[singlets] - values[singlets]) > bin_var) & (out[singlets] > quantile(out, 0.9, na.rm = T))]
  peak_singlets <- peak_singlets[!is.na(peak_singlets)]
  out[peak_singlets] <- values[peak_singlets]

  # Singlets which are within the main distribution are likely errors, so smooth them out
  main_singlets <- singlets[out[singlets] < quantile(out, 0.9, na.rm = T)]
  main_singlets <- main_singlets[!is.na(main_singlets)]
  main_singlets[main_singlets >= length(out)] <- length(out)-1
  out[main_singlets] <- out[main_singlets + 1]
  return(out)
}

#' calc_madOffset
#'
#' @param values vector of values to calculate the MAD offset
#' @return MAD offset
calc_madOffset <- function(values){
  values <- values[!is.na(values)]
  offset <- utils::tail(values, -1)
  values <- utils::head(values, -1)
  return(stats::mad(sqrt(2)*abs(offset - values)))
}

#' load.preprocess.bed
#' @import data.table
#'
#' @param bedpe_file path to the bedpe file produced by the preprocessing pipeline
#' @param bin_data bin data from the QDNASeq object
#' @param min_length minimum read length (default 30 nt)
#' @param max_length maximum read length (default 1000 nt)
#'
#' @return a data frame set up as a traditional bed file
load.preprocess.bed <- function(bedpe_file, bin_data, tag_overlap = 10){
  # Load bed and get QC and Ploidy Estimation Metrics
  bedpe <- tryCatch({data.table::fread(bedpe_file, sep = '\t')},
                    error = function(e) {
                      warning(paste0("Empty Bedpe file :", bedpe_file, ". Returning empty data frame."))
                      return(data.frame(Chr1 = character(0), Start1 = integer(0), End1 = integer(0),
                                        Chr2 = character(0), Start2 = integer(0), End2 = integer(0),
                                        Name = character(0), Score = numeric(0),
                                        R1_direction = character(0), R2_direction = character(0)))})
  if (nrow(bedpe) == 0) {return(bedpe)}
  colnames(bedpe) <-  c("Chr1", "Start1", "End1", "Chr2", "Start2", "End2",
                        "Name", "Score", "R1_direction", "R2_direction")

  # Create bed format
  bed <- bedpe[, .(Chr   = Chr1, Start = pmin(Start1, Start2), End = pmax(End1, End2), Name  = paste0(Name, ":", R1_direction, "/", R2_direction))]

  # Grab the first three characters of the chromosome names to see if its prefixed properly. If not add chr
  chr_prefix <- substr(bed$Chr, 1, 3)
  if(sum(grepl('chr', chr_prefix)) == 0){
    bed$Chr <- paste0('chr', bed$Chr)
  }

  bed <- bed[grep('(chr[0-9]+|X|Y)$', bed$Chr),] # Remove all decoy and alt contigs

  # Count doublets prior to trimming the Tn5 overlap
  bed$doublet <- is.doublet(bed, min.tag.overlap = tag_overlap-1, max.tag.overlap = tag_overlap)
  bed$End <- bed$End - tag_overlap
  bed[, Length := End - Start]
  bed[, Strandedness := extractStrandedness(Name, 3)]

  # Remove the artifactual duplications
  bed <- filter.bed(bed, bin_data)
  return(bed)
}

#' filter.bed
#'
#' @param bed a data frame formatted as a bed file with read positions
#' @param bin_data the metadata generated by QDNAseq for the bins
#'
#' @return a data frame with duplicate reads (identical start or end sites) removed
filter.bed <- function(bed, bin_data) {
  # Remove duplicated reads that aren't caught by the standard dedup tools
  bed <- bed[with(bed, order(Chr, Start, End, Strandedness)),]
  tmp.start <- paste(bed$Chr, bed$Start, bed$Strandedness)
  tmp.end <- paste(bed$Chr, bed$End, bed$Strandedness)

  are.duplicated <- duplicated(tmp.start, fromLast = TRUE) | duplicated(tmp.end)
  bed <- bed[!are.duplicated,]
  return(bed)
}

#' extractStrandedness
#'
#' @param readName name for each read pulled from the bed file
#' @param n number of characters to count back to get strandedness
#'
#' @return read orientation from the last n characters of the read name
extractStrandedness <- function(readName, n){
  readName.count <- nchar(readName)
  return(substr(readName, readName.count - n + 1, readName.count))
}

#' is.doublet
#'
#' @param bed a data frame formatted as a bed file with read positions must be filtered by filter.bed
#' @param min.tag.overlap smallest length of the tagmentation read overlap (default 9 nt)
#' @param max.tag.overlap largest length of the tagmentation read overlap (default 10 nt)
#'
#' @return a boolean for each read if it is a doublet, i.e. comes from 2 successful sequential tagmentation events
is.doublet <- function(bed, min.tag.overlap = 9, max.tag.overlap = 10){
  read.starts <- GenomicRanges::GRanges(seqnames = bed$Chr,
                                        ranges = IRanges::IRanges(start = bed$Start+min.tag.overlap, end = bed$Start+max.tag.overlap))
  read.ends <- GenomicRanges::GRanges(seqnames = bed$Chr,
                                      ranges = IRanges::IRanges(start = bed$End, end = bed$End))

  return(GenomicRanges::countOverlaps(read.starts, read.ends)>0)
}

#' count.overlaps
#'
#' @param bed a data frame formatted as a bed file with read positions must be filtered by filter.bed
#' @param min.size minimum read length (default 50 nt)
#' @param max.size maximum read length (default 1000 nt)
#' @param tag.overlap largest length of the tagmentation read overlap (default 10 nt)
#'
#' @return the bed file with the number of reads overlapping the upstream and over regions for each read
count.overlaps <- function(bed, min.size = 50, max.size = 1000) {
  # bed <- tmp2
  upstream.ranges <- IRanges::IRanges(start = bed$Start - max.size,
                                      end = bed$Start - min.size - 1)
  upstream.counting <- GenomicRanges::GRanges(seqnames = bed$Chr,
                                              ranges = upstream.ranges)

  overlap.ranges <- IRanges::IRanges(start = bed$Start - min.size,
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
                         Norm.Count.Upstream = counted.upstream.overlaps / (max.size - min.size),
                         Norm.Count.Over = counted.over.overlaps / (min.size))
  data.out <- cbind(bed, data.out)
  return(data.out)
}

#' calc.breadth
#'
#' @param bed imported bed data frame from bedpe object
#'
#' @return the breadth of the bed file
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

#' estimate.ploidy
#' @import data.table
#'
#' @param sample path to bedpe file
#' @param binSize binning size in nucleotides
#' @param genome genome version to use for binning
#' @param min_length minimum read length (default 50 nt)
#' @param max_length maximum read length (default 1000 nt)
#' @param tag_overlap number of nucleotides that the Tn5 overlaps (default 9 nt)
#'
#' @return a data frame with observations derived from the bedpe file
estimate.ploidy <- function(sample, bin_data, genome, min_length = 50, max_length = NULL, tag_overlap = 10){

  bin_data$binName <- paste0('chr', bin_data$chromosome, '_', bin_data$start)

  # Load bed and get QC and Ploidy Estimation Metrics
  bed <- load.preprocess.bed(sample, bin_data, tag_overlap)
  if(nrow(bed)==0){
    out <- data.frame(ratio = NA, breadth = NA, coverage = NA, prop_doublet_tags = NA,
                      avg_length = NA, overlap_genome_size = NA, ploidy_readCount = NA)
    return(out)
  }

  if(is.null(max_length)){
    max_length <- min(max(bed$Length, na.rm = T), quantile(bed$Length, .999)) # Strip out any extremely long reads (misaligns)
  }
  bed <- bed[(bed$Length>min_length) & (bed$Length<max_length),]
  bed <- count.overlaps(bed, min.size = min_length, max.size = max_length)

  # Merge the bed data with the bin data
  binSize <- max(bin_data$end-bin_data$start) + 1
  bed$binStart <- floor(bed$Start/binSize)*binSize+1
  bed$binName <- paste0(bed$Chr, '_', bed$binStart)

  bed_summary <- bed[
    , .(
      Bin.Reads           = .N,
      Count.Upstream      = sum(Count.Upstream,      na.rm = TRUE),
      Norm.Count.Upstream = sum(Norm.Count.Upstream, na.rm = TRUE),
      Count.Over          = sum(Count.Over,          na.rm = TRUE),
      Norm.Count.Over     = sum(Norm.Count.Over,     na.rm = TRUE),
      Bin.Coverage        = sum(Length,              na.rm = TRUE),
      Avg.Length          = mean(Length,             na.rm = TRUE),
      Num.Doublets        = sum(doublet,             na.rm = TRUE)
    ),
    by = binName
  ]

  # Graveyard of attempts to calculate breadth per bin fast
  #bed_ranges <- GenomicRanges::GRanges(seqnames = bed$Chr, ranges = IRanges::IRanges(start = bed$Start, end = bed$End))
  #mcols(bed_ranges)$binName <- bed$binName
  #gr_list <- split(bed_ranges, mcols(bed_ranges)$binName)
  #breadths <- lapply(gr_list, calc.breadth.fast)
  #calc.breadth.fast <- function(gr) sum(width(reduce(gr)))
  #breadths <- vapply(gr_list, calc.breadth.fast, numeric(1))

  match_idx <- match(bin_data$binName, bed_summary$binName)
  bin_data <- cbind(bin_data, bed_summary[match_idx,-c('binName')])


  # Get Doublet Proportion per bin & calc prop available genome
  bin_data$Bin.Coverage <- bin_data$Bin.Coverage / (bin_data$end - bin_data$start)
  bin_data$Bin.Prop.Doublets <- bin_data$Num.Doublets / bin_data$Bin.Reads
  bin_data$Bin.Quality <- bin_data$Bin.Coverage/bin_data$Bin.Prop.Doublets

  # Estimate the ploidy using the top 50th percentile bins
  bin_selection <-  (bin_data$mappability >= 90) & (bin_data$bases >= 90) & bin_data$use
  selector <- (bin_data$bin.depth > median(bin_data$bin.depth, na.rm = T)) & bin_selection

  # Resorting to this hybrid method until we can figure the math for this a bit better
  if(mean(bin_data$Bin.Coverage, na.rm = T) < 0.075){
    over_counts <- bin_data$Norm.Count.Over[selector]
    upstream_counts <- bin_data$Norm.Count.Upstream[selector]
    over_counts[!is.finite(over_counts)] <- NA
    upstream_counts[!is.finite(upstream_counts)] <- NA
    ratios <- mean(over_counts, na.rm = T) / mean(upstream_counts, na.rm = T)
  }else{
    ratios <- bin_data$Norm.Count.Over / bin_data$Norm.Count.Upstream
    ratios[!selector] <- NA
    ratios[!is.finite(ratios)] <- NA
  }
  bin_data$bin_ratio <- ratios
  bin_data$est_ploidy <- 1/(1-(mean(bin_data$bin_ratio, na.rm = T)))

  return(bin_data)
}
