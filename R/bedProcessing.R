#' process.batch
#'
#' @param bams path to bam files
#' @param genome genome assembly to use (e.g. hg38, hg19)
#' @param bedpes path to accompanying bedpe files
#' @param bin.size number of nucleotides per bin
#' @param min_length minimum length for reads
#' @param max_length maximum length for reads
#' @param tag_overlap size of the read overlap created by a tagmentation event
#' @param n_cpu number of CPUs to use for parallel processing
#' @param focal_amps whether to correct for focal amplifications
#' @param background_method background window method: "original" (narrow upstream) or
#'   "symmetric" (symmetric window around each read start). Default "original".
#' @param symmetric_x half-width of the symmetric background window in bp. Default 25000.
#' @param background_exclusion inner exclusion radius of the symmetric window in bp.
#'   If NULL, defaults to max_length. Use optimise_symmetric_x() to find empirically.
#' @param near_start_exclusion bp to exclude immediately upstream of each fragment
#'   start from the doublet overlap window. Default 0 (original behaviour).
#'
#' @return songbird object
#' @export
#'
#' @examples
process.batch <- function(bams = NULL, genome = 'hg38', bedpes = NULL, bin.size = 500000,
                          min_length = 50, max_length = NULL, tag_overlap = 10, n_cpu = NULL,
                          focal_amps = TRUE,
                          background_method    = c('original', 'symmetric'),
                          symmetric_x          = 25000,
                          background_exclusion = NULL,
                          near_start_exclusion = 0L) {

  background_method <- match.arg(background_method)

  if (is.null(n_cpu)) {
    n_cpu <- max(1L, parallel::detectCores() - 1L)
  }

  if (is.null(bams) & is.null(bedpes)) {
    stop('Either bams or bedpes must be provided')
  }

  if (!genome %in% c('hg38', 'hg19', 'hs1')) {
    stop('Genome must be one of hg38, hg19, or hs1')
  }

  if (is.null(bedpes)) {
    res <- parallel::mclapply(1:length(bams), function(i)
      process.bam(bams[i],
                  genome      = genome,
                  bin.size    = bin.size,
                  min_length  = min_length,
                  tag_overlap = tag_overlap,
                  focal_amps  = focal_amps),
      mc.cores = n_cpu)
  } else if (is.null(bams)) {
    res <- parallel::mclapply(1:length(bedpes), function(i)
      process.bedpe(genome               = genome,
                    bedpe                = bedpes[i],
                    bin.size             = bin.size,
                    min_length           = min_length,
                    tag_overlap          = tag_overlap,
                    focal_amps           = focal_amps,
                    background_method    = background_method,
                    symmetric_x          = symmetric_x,
                    background_exclusion = background_exclusion,
                    near_start_exclusion = near_start_exclusion),
      mc.cores = n_cpu)
  } else {
    res <- parallel::mclapply(1:length(bams), function(i)
      process.bam.bedpe(bam                  = bams[i],
                        genome               = genome,
                        bedpe                = bedpes[i],
                        bin.size             = bin.size,
                        min_length           = min_length,
                        tag_overlap          = tag_overlap,
                        focal_amps           = focal_amps,
                        background_method    = background_method,
                        symmetric_x          = symmetric_x,
                        background_exclusion = background_exclusion,
                        near_start_exclusion = near_start_exclusion),
      mc.cores = n_cpu)
  }

  # Detect worker errors
  is_error <- sapply(res, inherits, 'try-error')
  if (any(is_error)) {
    stop("process.batch: ", sum(is_error), " cell(s) failed.\nFirst error: ",
         conditionMessage(attr(res[[which(is_error)[1]]], 'condition')))
  }

  sce <- create_sce(res)
  S4Vectors::metadata(sce)$min_length           <- min_length
  S4Vectors::metadata(sce)$tag_overlap          <- tag_overlap
  S4Vectors::metadata(sce)$background_method    <- background_method
  S4Vectors::metadata(sce)$symmetric_x          <- symmetric_x
  S4Vectors::metadata(sce)$background_exclusion <- background_exclusion
  S4Vectors::metadata(sce)$near_start_exclusion <- near_start_exclusion
  S4Vectors::metadata(sce)$overlap_window_size  <- min_length - near_start_exclusion
  return(sce)
}

#' process.bam.bedpe
#'
#' @param bam path to a bam file
#' @param genome genome assembly to use
#' @param bedpe path to the accompanying bedpe file
#' @param bin.size number of nucleotides per bin
#' @param min_length minimum length for reads
#' @param max_length maximum length for reads
#' @param tag_overlap size of the Tn5 tagmentation overlap
#'
#' @return corrected reads data frame
process.bam.bedpe <- function(bam, genome, bedpe, bin.size = 500000, min_length = 50,
                              max_length = NULL, tag_overlap = 10, focal_amps = TRUE,
                              background_method    = 'original',
                              symmetric_x          = 25000,
                              background_exclusion = NULL,
                              near_start_exclusion = 0L) {
  reads     <- load_cell(bam, binSize = bin.size, genome)
  reads.cor <- convert_long(reads)
  reads.cor$reads <- (reads.cor$reads / sum(reads.cor$reads, na.rm = TRUE)) * 10000

  reads.cor$ubh_tx     <- ubh_segment(reads.cor$reads)
  reads.cor$bin.depth  <- reads.cor$reads / reads.cor$ubh_tx
  reads.cor$bam_file   <- bam
  reads.cor$bedpe_file <- bedpe
  reads.cor <- estimate.ploidy(bedpe                = bedpe,
                               bin_data             = reads.cor,
                               min_length           = min_length,
                               max_length           = max_length,
                               tag_overlap          = tag_overlap,
                               focal_amps           = focal_amps,
                               background_method    = background_method,
                               symmetric_x          = symmetric_x,
                               background_exclusion = background_exclusion,
                               near_start_exclusion = near_start_exclusion)
  return(reads.cor)
}

#' process.bam
#'
#' @param bam path to a bam file
#' @param genome genome assembly to use
#' @param bin.size number of nucleotides per bin
#' @param min_length minimum length for reads
#' @param max_length maximum length for reads
#' @param tag_overlap size of the Tn5 tagmentation overlap
#'
#' @return corrected reads data frame
process.bam <- function(bam, genome, bin.size = 500000, min_length = 50,
                        max_length = NULL, tag_overlap = 10, focal_amps = TRUE) {
  reads     <- load_cell(bam, binSize = bin.size, genome)
  reads.cor <- convert_long(reads)
  reads.cor$reads <- (reads.cor$reads / sum(reads.cor$reads, na.rm = TRUE)) * 10000

  reads.cor$ubh_tx     <- ubh_segment(reads.cor$reads)
  reads.cor$bin.depth  <- reads.cor$reads / reads.cor$ubh_tx
  reads.cor$bam_file   <- bam
  reads.cor$bedpe_file <- NA
  reads.cor$est_ploidy <- NA
  return(reads.cor)
}


#' process.bedpe
#'
#' @param genome genome assembly to use
#' @param bedpe path to the bedpe file
#' @param bin.size number of nucleotides per bin
#' @param min_length minimum length for reads
#' @param max_length maximum length for reads
#' @param tag_overlap size of the Tn5 tagmentation overlap
#'
#' @return corrected reads data frame
process.bedpe <- function(genome, bedpe, bin.size = 500000, min_length = 50,
                          max_length = NULL, tag_overlap = 10, focal_amps = TRUE,
                          background_method    = 'original',
                          symmetric_x          = 25000,
                          background_exclusion = NULL,
                          near_start_exclusion = 0L) {
  bins <- QDNAseq::getBinAnnotations(binSize = bin.size / 1000, genome = genome, verbose = FALSE)
  bins@data$mappability <- as.numeric(bins@data$mappability)
  reads.cor <- estimate.ploidy(bedpe                = bedpe,
                               bin_data             = bins,
                               min_length           = min_length,
                               max_length           = max_length,
                               tag_overlap          = tag_overlap,
                               focal_amps           = focal_amps,
                               background_method    = background_method,
                               symmetric_x          = symmetric_x,
                               background_exclusion = background_exclusion,
                               near_start_exclusion = near_start_exclusion)

  reads.cor$ubh_tx     <- ubh_segment(reads.cor$reads)
  reads.cor$bin.depth  <- reads.cor$reads / reads.cor$ubh_tx
  reads.cor$bedpe_file <- bedpe
  return(reads.cor)
}

#' estimate.ploidy
#' @import data.table
#'
#' @param bedpe path to bedpe file
#' @param bin_data bin data from the QDNASeq object or S4 bins object
#' @param min_length minimum read length (default 50 nt)
#' @param max_length maximum read length
#' @param tag_overlap number of nucleotides that the Tn5 overlaps (default 10)
#' @param focal_amps logical — correct for focal amplifications
#' @param background_method "original" or "symmetric"
#' @param symmetric_x half-width of symmetric background window in bp
#' @param background_exclusion inner exclusion radius of symmetric window in bp
#' @param near_start_exclusion bp to exclude near fragment start from doublet window
#'
#' @return a data frame with per-bin observations
estimate.ploidy <- function(bedpe, bin_data, min_length, max_length, tag_overlap,
                            focal_amps, genome = NULL,
                            background_method    = 'original',
                            symmetric_x          = 25000,
                            background_exclusion = NULL,
                            near_start_exclusion = 0L) {
  no_bam <- FALSE
  if (typeof(bin_data) == 'S4') {
    bins     <- bin_data
    bin_data <- bin_data@data
    no_bam   <- TRUE
  }

  bin_data$binName <- paste0('chr', bin_data$chromosome, '_', bin_data$start)

  bed <- load.preprocess.bed(bedpe, bin_data, tag_overlap)
  if (nrow(bed) == 0) {
    out <- data.frame(ratio = NA, breadth = NA, coverage = NA, prop_doublet_tags = NA,
                      avg_length = NA, overlap_genome_size = NA, ploidy_readCount = NA)
    return(out)
  }

  if (is.null(max_length)) {
    max_length <- min(max(bed$Length, na.rm = TRUE), quantile(bed$Length, .99))
  }
  bed <- bed[Length > min_length & Length < max_length]

  if (nrow(bed) == 0) {
    warning("No fragments remaining after length filtering for: ", bedpe)
    return(bin_data)
  }

  bed <- count.overlaps(bed, min.size = min_length, max.size = max_length,
                        background_method    = background_method,
                        symmetric_x          = symmetric_x,
                        background_exclusion = background_exclusion,
                        near_start_exclusion = near_start_exclusion)
  bed[, OlapCount := Count.Over + Count.Upstream]

  binSize <- max(bin_data$end - bin_data$start) + 1
  bed[, binStart := floor(Start / binSize) * binSize + 1L]
  bed[, binName  := paste0(Chr, '_', binStart)]

  bed_summary <- bed[
    , .(
      Bin.Reads           = .N,
      Count.Upstream      = sum(Count.Upstream,      na.rm = TRUE),
      Norm.Count.Upstream = sum(Norm.Count.Upstream, na.rm = TRUE),
      Count.Over          = sum(Count.Over,          na.rm = TRUE),
      Norm.Count.Over     = sum(Norm.Count.Over,     na.rm = TRUE),
      Bin.Coverage        = sum(Length,              na.rm = TRUE),
      Avg.Length          = mean(Length,             na.rm = TRUE),
      Num.Doublets        = sum(doublet,             na.rm = TRUE),
      Olap.Count          = sum(OlapCount,           na.rm = TRUE)
    ),
    by = binName
  ]

  match_idx <- match(bin_data$binName, bed_summary$binName)
  bin_data  <- cbind(bin_data, as.data.frame(bed_summary)[match_idx, -1])

  if (no_bam) {
    counts <- matrix(bin_data$Bin.Reads, ncol = 1)
    counts[is.na(counts)] <- 0
    pheno_data <- data.frame(
      name       = bedpe,
      reads      = sum(bin_data$Bin.Reads, na.rm = TRUE),
      used.reads = sum(bin_data$Bin.Reads, na.rm = TRUE),
      row.names  = bedpe,
      stringsAsFactors = FALSE
    )
    pheno_data_anno <- Biobase::AnnotatedDataFrame(pheno_data)
    reads <- new("QDNAseqReadCounts",
                 bins      = bins,
                 counts    = counts,
                 phenodata = pheno_data_anno)
    reads <- QDNAseq::applyFilters(reads, verbose = FALSE)
    reads <- QDNAseq::estimateCorrection(reads)
    reads <- QDNAseq::applyFilters(reads, chromosomes = NA, verbose = FALSE)
    reads <- QDNAseq::correctBins(reads)
    bin_data$reads             <- reads@assayData$copynumber
    bin_data$uncorrected.reads <- bin_data$Bin.Reads
  }

  bin_data$Bin.Coverage      <- bin_data$Bin.Coverage / (bin_data$end - bin_data$start)
  bin_data$Bin.Prop.Doublets <- bin_data$Num.Doublets / bin_data$Bin.Reads
  bin_data$Bin.Quality       <- bin_data$Bin.Coverage / bin_data$Bin.Prop.Doublets

  bin_selection <- (bin_data$mappability >= 90) & (bin_data$bases >= 90) & bin_data$use

  if (focal_amps) {
    depth_selection <- (bin_data$Olap.Count > quantile(bin_data$Olap.Count, 0.1, na.rm = TRUE)) &
                       (bin_data$Olap.Count < quantile(bin_data$Olap.Count, 0.9, na.rm = TRUE))
    if (sum(depth_selection, na.rm = TRUE) < 50) {
      depth_selection <- rep(TRUE, nrow(bin_data))
    }
  } else {
    depth_selection <- (bin_data$bin.depth > median(bin_data$bin.depth, na.rm = TRUE)) &
                       bin_selection
  }

  selector <- depth_selection & bin_selection
  selector[is.na(selector)] <- FALSE

  if (mean(bin_data$Bin.Coverage, na.rm = TRUE) < 0.075) {
    over_counts     <- bin_data$Norm.Count.Over[selector]
    upstream_counts <- bin_data$Norm.Count.Upstream[selector]
    over_counts[!is.finite(over_counts)]         <- NA
    upstream_counts[!is.finite(upstream_counts)] <- NA
    ratios <- mean(over_counts, na.rm = TRUE) / mean(upstream_counts, na.rm = TRUE)
  } else {
    ratios <- bin_data$Norm.Count.Over / bin_data$Norm.Count.Upstream
    ratios[!selector]          <- NA
    ratios[!is.finite(ratios)] <- NA
  }

  bin_data$bin_ratio  <- ratios
  bin_data$est_ploidy <- 1 / (1 - mean(bin_data$bin_ratio, na.rm = TRUE))
  return(bin_data)
}

#' convert_long
#'
#' @param reads reads object from QDNASeq
#'
#' @return table of reads for each cell
convert_long <- function(reads) {
  uncorReads <- reads$reads
  reads      <- reads$reads.cor
  metadata   <- reads@featureData@data

  table     <- reads@assayData$copynumber
  out_table <- c()
  for (i in 1:ncol(table)) {
    new_dat         <- metadata
    new_dat$reads   <- table[, i]
    new_dat$cell_id <- colnames(table)[i]
    out_table       <- rbind(out_table, new_dat)
  }

  table                       <- uncorReads@assayData$counts
  out_table$uncorrected.reads <- table[, 1]
  return(out_table)
}


#' load_cell
#'
#' @param bamPath path to the bam file
#' @param binSize size of the bins in nucleotides
#' @param genome genome to use for binning
#'
#' @return list with reads and corrected reads
load_cell <- function(bamPath, binSize, genome) {
  bins <- QDNAseq::getBinAnnotations(binSize = binSize / 1000, genome = genome, verbose = FALSE)
  bins@data$mappability <- as.numeric(bins@data$mappability)
  reads     <- QDNAseq::binReadCounts(bins, bamfiles = bamPath, pairedEnds = TRUE, verbose = FALSE)
  reads     <- QDNAseq::applyFilters(reads, verbose = FALSE)
  reads     <- QDNAseq::estimateCorrection(reads)
  reads     <- QDNAseq::applyFilters(reads, chromosomes = NA, verbose = FALSE)
  reads.cor <- QDNAseq::correctBins(reads)
  return(list(reads = reads, reads.cor = reads.cor))
}

#' ubh_segment
#'
#' @param values GC and mappability corrected reads per bin
#'
#' @return wavelet-segmented signal
ubh_segment <- function(values) {

  out         <- rep(NA, length(values))
  real_idxs   <- which(!is.na(values))
  filt_values <- values[!is.na(values)]
  init_sigma  <- calc_madOffset(filt_values)

  transform      <- unbalhaar::best.unbal.haar(filt_values)
  transform.filt <- unbalhaar::hard.thresh(transform, init_sigma)
  reconstr       <- unbalhaar::reconstr(transform.filt)
  reconstr[(filt_values == 0) | (reconstr < 0)] <- 0
  out[real_idxs] <- reconstr

  na_idx   <- which(is.na(values))
  next_idx <- na_idx + 1
  next_idx[next_idx > length(out)] <- length(out) - 1
  out[na_idx] <- out[next_idx]

  bin_var   <- var(out[(out > quantile(out, 0.05, na.rm = TRUE)) &
                       (out < quantile(out, 0.95, na.rm = TRUE))])
  shift_idx <- which(c(0, diff(out)) != 0)
  singlets  <- shift_idx[out[shift_idx] != out[shift_idx + 1]]

  peak_singlets <- singlets[(abs(out[singlets] - values[singlets]) > bin_var) &
                            (out[singlets] > quantile(out, 0.9, na.rm = TRUE))]
  peak_singlets <- peak_singlets[!is.na(peak_singlets)]
  out[peak_singlets] <- values[peak_singlets]

  main_singlets <- singlets[out[singlets] < quantile(out, 0.9, na.rm = TRUE)]
  main_singlets <- main_singlets[!is.na(main_singlets)]
  main_singlets[main_singlets >= length(out)] <- length(out) - 1
  out[main_singlets] <- out[main_singlets + 1]
  return(out)
}

#' calc_madOffset
#'
#' @param values vector of values
#' @return MAD offset
calc_madOffset <- function(values) {
  values <- values[!is.na(values)]
  offset <- utils::tail(values, -1)
  values <- utils::head(values, -1)
  return(stats::mad(sqrt(2) * abs(offset - values)))
}

#' load.preprocess.bed
#' @import data.table
#'
#' @param bedpe_file path to the bedpe file
#' @param bin_data bin data from QDNASeq
#' @param tag_overlap Tn5 footprint size in bp
#'
#' @return a data.table in bed format
load.preprocess.bed <- function(bedpe_file, bin_data, tag_overlap) {
  bedpe <- tryCatch(
    data.table::fread(bedpe_file, sep = '\t'),
    error = function(e) {
      warning(paste0("Empty Bedpe file: ", bedpe_file, ". Returning empty data frame."))
      return(data.frame(Chr1 = character(0), Start1 = integer(0), End1 = integer(0),
                        Chr2 = character(0), Start2 = integer(0), End2 = integer(0),
                        Name = character(0), Score = numeric(0),
                        R1_direction = character(0), R2_direction = character(0)))
    })
  if (nrow(bedpe) == 0) { return(bedpe) }
  colnames(bedpe) <- c("Chr1", "Start1", "End1", "Chr2", "Start2", "End2",
                       "Name", "Score", "R1_direction", "R2_direction")

  bed <- bedpe[, .(Chr   = Chr1,
                   Start = pmin(Start1, Start2),
                   End   = pmax(End1, End2),
                   Name  = paste0(Name, ":", R1_direction, "/", R2_direction))]

  chr_prefix <- substr(bed$Chr, 1, 3)
  if (sum(grepl('chr', chr_prefix)) == 0) {
    bed[, Chr := paste0('chr', Chr)]
  }

  bed <- bed[grepl('(chr[0-9]+|X|Y)$', Chr)]

  # Count doublets BEFORE trimming the Tn5 overlap
  bed[, doublet := is.doublet(bed, min.tag.overlap = tag_overlap - 1,
                                   max.tag.overlap = tag_overlap)]
  bed[, End    := End - tag_overlap]
  bed[, Length := End - Start]
  bed[, Strandedness := extractStrandedness(Name, 3)]

  bed <- filter.bed(bed, bin_data)
  return(bed)
}

#' filter.bed
#'
#' @param bed a data.table in bed format
#' @param bin_data the metadata generated by QDNAseq for the bins
#'
#' @return a data.table with duplicate reads removed
filter.bed <- function(bed, bin_data) {
  data.table::setorder(bed, Chr, Start, End, Strandedness)
  tmp.start <- paste(bed$Chr, bed$Start, bed$Strandedness)
  tmp.end   <- paste(bed$Chr, bed$End,   bed$Strandedness)

  are.duplicated <- duplicated(tmp.start, fromLast = TRUE) | duplicated(tmp.end)
  bed <- bed[!are.duplicated]
  return(bed)
}

#' extractStrandedness
#'
#' @param readName name for each read
#' @param n number of characters from end
#'
#' @return read orientation string
extractStrandedness <- function(readName, n) {
  readName.count <- nchar(readName)
  return(substr(readName, readName.count - n + 1, readName.count))
}

#' is.doublet
#'
#' @param bed a data.table in bed format
#' @param min.tag.overlap smallest overlap length
#' @param max.tag.overlap largest overlap length
#'
#' @return logical vector — TRUE if read is a tagmentation doublet
is.doublet <- function(bed, min.tag.overlap, max.tag.overlap) {
  read.starts <- GenomicRanges::GRanges(
    seqnames = bed$Chr,
    ranges   = IRanges::IRanges(start = bed$Start + min.tag.overlap,
                                end   = bed$Start + max.tag.overlap))
  read.ends <- GenomicRanges::GRanges(
    seqnames = bed$Chr,
    ranges   = IRanges::IRanges(start = bed$End, end = bed$End))
  return(GenomicRanges::countOverlaps(read.starts, read.ends) > 0)
}

#' count.overlaps
#'
#' @param bed a data.table in bed format
#' @param min.size minimum read length
#' @param max.size maximum read length
#' @param background_method "original" (narrow upstream) or "symmetric"
#' @param symmetric_x half-width of symmetric background window in bp. Default 25000.
#' @param background_exclusion inner exclusion radius in bp. If NULL, uses max.size.
#' @param near_start_exclusion bp to exclude near start from doublet window. Default 0.
#'
#' @return the bed data.table with Count.Upstream, Count.Over, Norm.* columns added
count.overlaps <- function(bed, min.size, max.size,
                           background_method    = c('original', 'symmetric'),
                           symmetric_x          = 25000,
                           background_exclusion = NULL,
                           near_start_exclusion = 0L) {

  background_method <- match.arg(background_method)

  # Fragment starts for counting (sized fragments only)
  bed.sized    <- bed[Length > min.size & Length < max.size]
  for.counting <- GenomicRanges::GRanges(
    seqnames = bed.sized$Chr,
    ranges   = IRanges::IRanges(start = bed.sized$Start, end = bed.sized$Start))

  # ---- Doublet overlap window ----------------------------------------------
  # [Start - min.size, Start - near_start_exclusion - 1]
  # near_start_exclusion removes reads very close to the reference (same molecule risk)
  over_end   <- bed$Start - as.integer(near_start_exclusion) - 1L
  over_start <- bed$Start - as.integer(min.size)
  over_valid <- over_end >= over_start

  over.counting <- GenomicRanges::GRanges(
    seqnames = bed$Chr[over_valid],
    ranges   = IRanges::IRanges(start = over_start[over_valid],
                                end   = over_end[over_valid]))
  counted.over.overlaps <- rep(0L, nrow(bed))
  counted.over.overlaps[over_valid] <- GenomicRanges::countOverlaps(
    over.counting, for.counting)
  overlap_window_size <- as.integer(min.size) - as.integer(near_start_exclusion)

  # ---- Background window --------------------------------------------------
  if (background_method == 'original') {

    # Classic narrow upstream window: [Start - max.size, Start - min.size - 1]
    upstream_start <- pmax(1L, bed$Start - max.size)
    upstream_end   <- bed$Start - min.size - 1L
    valid          <- upstream_end >= upstream_start
    upstream_start[!valid] <- 1L
    upstream_end[!valid]   <- 1L

    upstream.counting <- GenomicRanges::GRanges(
      seqnames = bed$Chr,
      ranges   = IRanges::IRanges(start = upstream_start, end = upstream_end))
    counted.upstream.overlaps <- GenomicRanges::countOverlaps(
      upstream.counting, for.counting)
    counted.upstream.overlaps[!valid] <- NA_integer_
    window_size <- max.size - min.size

  } else {

    # Symmetric window with central exclusion zone.
    # background_exclusion = inner boundary (where accessibility decay ends).
    # symmetric_x          = outer boundary (total window reach).
    # If background_exclusion is NULL, use max.size for backward compatibility.
    excl_zone <- if (!is.null(background_exclusion)) as.integer(background_exclusion)
                 else as.integer(max.size)

    # Left arm: [Start - symmetric_x, Start - excl_zone - 1]
    left_start <- pmax(1L, bed$Start - symmetric_x)
    left_end   <- bed$Start - excl_zone - 1L
    left_valid <- left_end >= left_start
    left_start_safe        <- left_start
    left_end_safe          <- left_end
    left_start_safe[!left_valid] <- 1L
    left_end_safe[!left_valid]   <- 1L

    left.counting <- GenomicRanges::GRanges(
      seqnames = bed$Chr,
      ranges   = IRanges::IRanges(start = left_start_safe, end = left_end_safe))
    left_counts <- GenomicRanges::countOverlaps(left.counting, for.counting)
    left_counts[!left_valid] <- 0L

    # Right arm: [Start + excl_zone, Start + symmetric_x]
    right_start <- bed$Start + excl_zone
    right_end   <- bed$Start + symmetric_x
    right_valid <- right_start <= right_end
    right.counting <- GenomicRanges::GRanges(
      seqnames = bed$Chr,
      ranges   = IRanges::IRanges(start = right_start, end = right_end))
    right_counts <- GenomicRanges::countOverlaps(right.counting, for.counting)
    right_counts[!right_valid] <- 0L

    counted.upstream.overlaps <- left_counts + right_counts

    left_width  <- pmax(0L, left_end  - left_start  + 1L)
    left_width[!left_valid]   <- 0L
    right_width <- pmax(0L, right_end - right_start + 1L)
    right_width[!right_valid] <- 0L
    window_size <- pmax(1L, left_width + right_width)
  }

  bed[, Count.Upstream      := counted.upstream.overlaps]
  bed[, Count.Over          := counted.over.overlaps]
  bed[, Norm.Count.Upstream := counted.upstream.overlaps / window_size]
  bed[, Norm.Count.Over     := counted.over.overlaps / overlap_window_size]
  return(bed)
}

#' calc.breadth
#'
#' @param bed imported bed data frame
#'
#' @return the breadth of the bed file
calc.breadth <- function(bed) {
  read_ranges <- GenomicRanges::GRanges(
    seqnames = bed$Chr,
    ranges   = IRanges::IRanges(start = bed$Start, end = bed$End))
  coverage <- GenomicRanges::coverage(read_ranges)
  breadth  <- 0
  for (i in 1:length(coverage)) {
    breadth <- breadth + sum(coverage[[i]]@lengths[coverage[[i]]@values > 0])
  }
  return(breadth)
}

#' optimise_symmetric_x
#'
#' Profiles the background read density decay as a function of distance from
#' fragment starts and returns the distance at which density plateaus —
#' i.e., the suggested \code{background_exclusion} for use in
#' \code{\link{process.batch}} with \code{background_method = "symmetric"}.
#'
#' @param bedpe_file character — path to a single .bedpe or .bedpe.gz file.
#' @param min_length integer — minimum fragment length. Default 50.
#' @param max_length integer — maximum fragment length. Default 1000.
#' @param tag_overlap integer — Tn5 footprint size in bp. Default 10.
#' @param max_dist    integer — maximum distance to profile in bp. Default 150000.
#' @param band_width  integer — width of each distance band in bp. Default 1000.
#' @param n_sample    integer — number of reference reads to sample. Default 5000.
#' @param plateau_tol numeric — fraction of far-field mean within which density is
#'   considered flat. Default 0.05 (5 percent).
#'
#' @return integer — suggested \code{background_exclusion} in bp.
#' @export
optimise_symmetric_x <- function(bedpe_file,
                                  min_length  = 50L,
                                  max_length  = 1000L,
                                  tag_overlap = 10L,
                                  max_dist    = 150000L,
                                  band_width  = 1000L,
                                  n_sample    = 5000L,
                                  plateau_tol = 0.05) {
  prof <- plot_background_decay(
    bedpe_file  = bedpe_file,
    min_length  = min_length,
    max_length  = max_length,
    tag_overlap = tag_overlap,
    max_dist    = max_dist,
    band_width  = band_width,
    n_sample    = n_sample,
    plateau_tol = plateau_tol,
    return_data = TRUE
  )
  far_field    <- prof$density[prof$mid_dist > max_dist * 0.8]
  ff_mean      <- mean(far_field, na.rm = TRUE)
  ff_upper     <- ff_mean * (1 + plateau_tol)
  plateau_band <- max(which(prof$density > ff_upper), na.rm = TRUE)
  suggested_x  <- as.integer(prof$mid_dist[plateau_band] + band_width / 2)
  message("Optimal background_exclusion: ", suggested_x, " bp")
  suggested_x
}
