#' load.preprocess.bed
#'
#' @param bedpe_file path to the bedpe file produced by the preprocessing pipeline
#' @param min_length minimum read length (default 30 nt)
#' @param max_length maximum read length (default 1000 nt)
#'
#' @return a data frame set up as a traditional bed file
#' @export
#'
#' @examples
load.preprocess.bed <- function(bedpe_file, min_length = 30, max_length = 1000){
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
  bed <- bed[grep('(chr[0-9]+|X|Y)$', bed$Chr),]
  return(bed)
}

#' Title
#'
#' @param bed a data frame formatted as a bed file with read positions
#'
#' @return a data frame with duplicate reads not caught by samtools removed
#' @export
#'
#' @examples
#'
remove.duplicates <- function(bed) {
  bed <- bed[with(bed, order(Chr, Start, End, Strandedness)),]
  tmp.start <- paste(bed$Chr, bed$Start, bed$Strandedness)
  tmp.end <- paste(bed$Chr, bed$End, bed$Strandedness)

  are.duplicated <- base::duplicated(tmp.start, fromLast = TRUE) | base::duplicated(tmp.end)
  return(bed[!are.duplicated,])
}

#' Title
#'
#' @param readName name for each read pulled from the bed file
#' @param n number of characters to count back to get strandedness
#'
#' @return read orientation from the last n characters of the read name
#' @export
#'
#' @examples
extractStrandedness <- function(readName, n){
  readName.count <- base::nchar(readName)
  return(base::substr(readName, readName.count - n + 1, readName.count))
}

#' Title
#'
#' @param bed a data frame formatted as a bed file with read positions must be filtered by remove.duplicates
#' @param tag.overlap smallest length of the tagmentation read overlap (default 9 nt)
#'
#' @return
#' @export
#'
#' @examples
count.doublets <- function(bed, tag.overlap = 9){
  read.starts <- GenomicRanges::GRanges(seqnames = bed$Chr,
                         ranges = IRanges::IRanges(start = bed$Start+tag.overlap, end = bed$Start+tag.overlap))
  read.ends <- GenomicRanges::GRanges(seqnames = bed$Chr,
                       ranges = IRanges::IRanges(start = bed$End, end = bed$End))

  num.doublets <- sum(GenomicRanges::countOverlaps(read.starts, read.ends)>0)
  return(num.doublets/nrow(bed))
}

#' Title
#'
#' @param bed a data frame formatted as a bed file with read positions must be filtered by remove.duplicates
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
estimate.ploidy <- function(sample, binSize, min_length = 50, max_length = 1000){
  bed <- load.preprocess.bed(sample, min_length, max_length)
  bed <- remove.duplicates(bed)

  prop_doublets <- count.doublets(bed)
  bed <- count.overlaps(bed, min.size = 100, max.size = max_length)

  bed$binStart <- base::floor(bed$Start/binSize)*binSize+1
  bed$binName <- base::paste0(bed$Chr, '_', bed$binStart)

  bed_aggregate <- stats::aggregate(bed[,c("Count.Upstream",
                                           "Norm.Count.Upstream",
                                           "Count.Over",
                                           "Norm.Count.Over")],
                             by = list(bed$binName), function(x) mean(x, na.rm = T))
  bed_counts <- stats::aggregate(bed$binName, by = list(bed$binName), length)
  bed_breadth <- stats::aggregate(bed$Length, by = list(bed$binName), function(x) sum(x, na.rm = T))


  bins <- QDNAseq::getBinAnnotations(binSize/1000, genome = 'hg38')
  bin_data <- bins@data
  bin_data$binName <- base::paste0('chr', bin_data$chromosome, '_', bin_data$start)
  bin_data$use[bin_data$mappability < 90 | bin_data$bases < 90] <- FALSE

  match_idx <- base::match(bin_data$binName, bed_aggregate$Group.1)
  bin_data <- base::cbind(bin_data, bed_aggregate[match_idx,])
  bin_data$reads <- bed_counts$x[match_idx]
  bin_data$breadth <- bed_breadth$x[match_idx]
  bin_data$breadth <- bin_data$breadth/binSize

  coverage <- bin_data$Coverage
  coverage[!is.finite(coverage)] <- NA
  coverage <- mean(coverage, na.rm = T)

  bin_data <- bin_data[,!grepl('Group', colnames(bin_data))]
  bin_data <- bin_data[bin_data$use,]
  bin_data <- bin_data[!is.na(bin_data$reads),]
}
