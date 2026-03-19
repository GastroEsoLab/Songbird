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
#' @param focal_amps whether to filter focal amplifications during ploidy estimation
#' @param ploidy_method character — ploidy estimation method to use downstream.
#'   \code{"ratio"} (default) uses the original global ratio average via
#'   \code{\link{ploidy_correction}}.
#'   \code{"hmm"} uses the Poisson-HMM segmentation via
#'   \code{\link{ploidy_correction_hmm}}.
#'   This argument is stored in \code{metadata(sce)$ploidy_method} so that
#'   downstream functions can dispatch correctly.
#' @param background_method character — background window design for overlap
#'   counting. \code{"original"} (default) uses the narrow upstream window
#'   \code{[Start - max_length, Start - min_length - 1]}, preserving the
#'   original behaviour exactly.
#'   \code{"symmetric"} uses a window extending \code{symmetric_x} bp on both
#'   sides of each fragment start, with the doublet zone excised, reducing
#'   sensitivity to local accessibility biases.
#'   Only relevant when bedpe files are provided.
#' @param symmetric_x integer — half-width (bp) of the symmetric background
#'   window. Only used when \code{background_method = "symmetric"}.
#'   Default 25000.
#'
#' @return songbird object
#' @export
#'
#' @examples

process.batch <- function(bams = NULL, genome = 'hg38', bedpes = NULL,
                          bin.size = 500000, min_length = 50, max_length = NULL,
                          tag_overlap = 10, n_cpu = NULL, focal_amps = TRUE,
                          ploidy_method    = c('ratio', 'hmm'),
                          background_method = c('original', 'symmetric'),
                          symmetric_x      = 25000) {

  ploidy_method     <- match.arg(ploidy_method)
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

  if (ploidy_method == 'hmm' && is.null(bedpes)) {
    stop("ploidy_method = 'hmm' requires bedpe files for overlap count data.\n",
         "Re-run with bedpes provided, or use ploidy_method = 'ratio'.")
  }

  if (background_method == 'symmetric' && is.null(bedpes)) {
    warning("background_method = 'symmetric' has no effect without bedpe files.")
  }

  if (is.null(bedpes)) {
    res <- parallel::mclapply(seq_along(bams), function(i)
      process.bam(bams[i],
                  genome      = genome,
                  bin.size    = bin.size,
                  min_length  = min_length,
                  tag_overlap = tag_overlap,
                  focal_amps  = focal_amps),
      mc.cores = n_cpu)
  } else if (is.null(bams)) {
    res <- parallel::mclapply(seq_along(bedpes), function(i)
      process.bedpe(genome             = genome,
                    bedpe              = bedpes[i],
                    bin.size           = bin.size,
                    min_length         = min_length,
                    tag_overlap        = tag_overlap,
                    focal_amps         = focal_amps,
                    background_method  = background_method,
                    symmetric_x        = symmetric_x),
      mc.cores = n_cpu)
  } else {
    message('Running bam + bedpe pipeline')
    res <- parallel::mclapply(seq_along(bams), function(i)
      process.bam.bedpe(bam              = bams[i],
                        genome           = genome,
                        bedpe            = bedpes[i],
                        bin.size         = bin.size,
                        min_length       = min_length,
                        tag_overlap      = tag_overlap,
                        focal_amps       = focal_amps,
                        background_method = background_method,
                        symmetric_x      = symmetric_x),
      mc.cores = n_cpu)
  }

  # Detect any worker errors from mclapply before passing to create_sce.
  # mclapply returns try-error objects for failed workers rather than stopping;
  # these cause cryptic $ errors inside reads_to_matrix.
  is_error <- sapply(res, inherits, 'try-error')
  if(any(is_error)){
    bad_idx <- which(is_error)
    bad_files <- basename(bams[bad_idx])
    stop(
      length(bad_idx), " cell(s) failed during processing:\n",
      paste0("  [", bad_idx, "] ", bad_files, collapse = "\n"), "\n",
      "First error message:\n  ",
      as.character(res[[bad_idx[1]]])
    )
  }

  # Also detect stub returns from estimate.ploidy (empty bedpe / bad length filter)
  # which return a data frame lacking the expected columns
  is_stub <- sapply(res, function(x) !is.data.frame(x) || !'chromosome' %in% colnames(x))
  if(any(is_stub)){
    bad_idx <- which(is_stub)
    bad_files <- basename(bams[bad_idx])
    stop(
      length(bad_idx), " cell(s) returned incomplete data (empty bedpe or no fragments after filtering):\n",
      paste0("  [", bad_idx, "] ", bad_files, collapse = "\n"), "\n",
      "Check that bedpe files exist and contain fragments longer than min_length = ", min_length, " bp."
    )
  }

  sce <- create_sce(res)

  # Store choices in metadata so downstream functions have full provenance
  S4Vectors::metadata(sce)$ploidy_method     <- ploidy_method
  S4Vectors::metadata(sce)$background_method <- background_method
  S4Vectors::metadata(sce)$symmetric_x       <- symmetric_x
  S4Vectors::metadata(sce)$tag_overlap       <- tag_overlap
  S4Vectors::metadata(sce)$min_length        <- min_length

  return(sce)
}

#' process.bam.bedpe
#'
#' @param bam path to a bam file
#' @param genome genome assembly to use (e.g. hg38, hg19)
#' @param bedpe path to the accompanying bedfile
#' @param bin.size number of nucleotides per bin
#' @param min_length minimum length for reads
#' @param max_length maximum length for reads
#' @param tag_overlap size of the read overlap created by a tagmentation event
#' @param background_method background window method — \code{"original"} or \code{"symmetric"}
#' @param symmetric_x half-width (bp) of the symmetric background window
#'
#' @return corrected reads
process.bam.bedpe <- function(bam, genome, bedpe, bin.size = 500000,
                              min_length = 50, max_length = NULL,
                              tag_overlap = 10, focal_amps = TRUE,
                              background_method = 'original',
                              symmetric_x = 25000){
  reads <- load_cell(bam, binSize = bin.size, genome)
  reads.cor <- convert_long(reads)
  reads.cor$reads <- (reads.cor$reads/sum(reads.cor$reads, na.rm = T))*10000
  reads.cor$ubh_tx <- ubh_segment(reads.cor$reads)
  reads.cor$bin.depth <- reads.cor$reads/reads.cor$ubh_tx
  reads.cor$bam_file <- bam
  reads.cor$bedpe_file <- bedpe
  reads.cor <- estimate.ploidy(bedpe = bedpe, bin_data = reads.cor,
                               min_length = min_length, max_length = max_length,
                               tag_overlap = tag_overlap, focal_amps = focal_amps,
                               background_method = background_method,
                               symmetric_x = symmetric_x)
  return(reads.cor)
}

#' process.bam
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
process.bam <- function(bam, genome, bin.size = 500000, min_length = 50, max_length = NULL, tag_overlap = 10, focal_amps = TRUE){
  reads <- load_cell(bam, binSize = bin.size, genome)
  reads.cor <- convert_long(reads)
  reads.cor$reads <- (reads.cor$reads/sum(reads.cor$reads, na.rm = T))*10000
  reads.cor$ubh_tx <- ubh_segment(reads.cor$reads)
  reads.cor$bin.depth <- reads.cor$reads/reads.cor$ubh_tx
  reads.cor$bam_file <- bam
  reads.cor$bedpe_file <- NA
  reads.cor$est_ploidy <- NA
  return(reads.cor)
}


#' process.bedpe
#'
#' @param genome genome assembly to use (e.g. hg38, hg19)
#' @param bedpe path to the accompanying bedfile
#' @param bin.size number of nucleotides per bin
#' @param min_length minimum length for reads
#' @param max_length maximum length for reads
#' @param tag_overlap size of the read overlap created by a tagmentation event
#' @param background_method background window method — \code{"original"} or \code{"symmetric"}
#' @param symmetric_x half-width (bp) of the symmetric background window
#'
#' @return corrected reads
process.bedpe <- function(genome, bedpe, bin.size = 500000,
                          min_length = 50, max_length = NULL,
                          tag_overlap = 10, focal_amps = TRUE,
                          background_method = 'original',
                          symmetric_x = 25000){
  bins <- QDNAseq::getBinAnnotations(binSize = bin.size/1000, genome = genome, verbose = F)
  bins@data$mappability <- as.numeric(bins@data$mappability)
  reads.cor <- estimate.ploidy(bedpe = bedpe, bin_data = bins,
                               min_length = min_length, max_length = max_length,
                               tag_overlap = tag_overlap, focal_amps = focal_amps,
                               background_method = background_method,
                               symmetric_x = symmetric_x)
  reads.cor$ubh_tx <- ubh_segment(reads.cor$reads)
  reads.cor$bin.depth <- reads.cor$reads/reads.cor$ubh_tx
  reads.cor$bedpe_file <- bedpe
  return(reads.cor)
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
#' @param background_method character — background window method passed to
#'   \code{\link{count.overlaps}}. \code{"original"} (default) or
#'   \code{"symmetric"}. See \code{\link{count.overlaps}} for details.
#' @param symmetric_x integer — half-width (bp) of the symmetric background
#'   window. Only used when \code{background_method = "symmetric"}.
#'   Default 25000.
#'
#' @return a data frame with observations derived from the bedpe file
estimate.ploidy <- function(bedpe, bin_data, min_length, max_length,
                            tag_overlap, focal_amps, genome = NULL,
                            background_method = 'original',
                            symmetric_x = 25000){
  no_bam <- FALSE
  if(typeof(bin_data)=='S4'){
    bins <- bin_data
    bin_data <- bin_data@data
    no_bam <- TRUE
  }

  bin_data$binName <- paste0('chr', bin_data$chromosome, '_', bin_data$start)

  # Load bed and get QC and Ploidy Estimation Metrics
  bed <- load.preprocess.bed(bedpe, bin_data, tag_overlap)
  if(nrow(bed)==0){
    out <- data.frame(ratio = NA, breadth = NA, coverage = NA, prop_doublet_tags = NA,
                      avg_length = NA, overlap_genome_size = NA, ploidy_readCount = NA)
    return(out)
  }

  if(is.null(max_length)){
    max_length <- min(max(bed$Length, na.rm = T), quantile(bed$Length, .99)) # Strip out any extremely long reads (misaligns)
  }
  bed <- bed[Length > min_length & Length < max_length]

  if(nrow(bed) == 0){
    warning(paste0("No fragments remaining after length filtering [", min_length, ", ", max_length, "] for: ", bedpe))
    bin_data$est_ploidy <- NA_real_
    return(bin_data)
  }

  # Count overlaps & strip the outliers to correct for focal amplifications
  bed <- count.overlaps(bed, min.size = min_length, max.size = max_length,
                        background_method = background_method,
                        symmetric_x = symmetric_x)
  bed[, OlapCount := Count.Over + Count.Upstream]

  # Merge the bed data with the bin data
  binSize <- max(bin_data$end-bin_data$start) + 1
  bed[, binStart := floor(Start/binSize)*binSize+1L]
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
      Total.Window.Size   = sum(Window.Size,         na.rm = TRUE),
      Num.Doublets        = sum(doublet,             na.rm = TRUE),
      Olap.Count          = sum(OlapCount,           na.rm = TRUE)
    ),
    by = binName
  ]

  match_idx <- match(bin_data$binName, bed_summary$binName)
  bin_data <- cbind(bin_data, as.data.frame(bed_summary)[match_idx, setdiff(names(bed_summary), "binName"), drop = FALSE])

  if(no_bam){
    counts <- matrix(bin_data$Bin.Reads, ncol = 1)
    is.na(counts) <- 0
    pheno_data <- data.frame(
      name = bedpe,
      reads = sum(bin_data$Bin.Reads, na.rm = T),
      used.reads = sum(bin_data$Bin.Reads, na.rm = T),
      row.names = bedpe,
      stringsAsFactors = FALSE
    )
    pheno_data_anno <- Biobase::AnnotatedDataFrame(pheno_data)

    reads <- new("QDNAseqReadCounts",
                 bins = bins,
                 counts = counts,
                 phenodata = pheno_data_anno)
    reads <- QDNAseq::applyFilters(reads, verbose = F)
    reads <- QDNAseq::estimateCorrection(reads)
    reads <- QDNAseq::applyFilters(reads, chromosomes = NA, verbose = F)
    reads <- QDNAseq::correctBins(reads)
    bin_data$reads <- reads@assayData$copynumber
    bin_data$uncorrected.reads <- bin_data$Bin.Reads
  }

  # Get Doublet Proportion per bin & calc prop available genome
  bin_data$Bin.Coverage <- bin_data$Bin.Coverage / (bin_data$end - bin_data$start)
  bin_data$Bin.Prop.Doublets <- bin_data$Num.Doublets / bin_data$Bin.Reads
  bin_data$Bin.Quality <- bin_data$Bin.Coverage/bin_data$Bin.Prop.Doublets

  # Estimate the ploidy using high quality bins
  bin_selection <-  (bin_data$mappability >= 90) & (bin_data$bases >= 90) & bin_data$use

  # If we suspect that there are focal amplifications, we will try to strip those
  if(focal_amps){
    depth_selection <- (bin_data$Olap.Count > quantile(bin_data$Olap.Count, 0.1, na.rm = T)) & (bin_data$Olap.Count < quantile(bin_data$Olap.Count, 0.9, na.rm = T))
    if(sum(depth_selection, na.rm = T) < 50){
      depth_selection <- rep(TRUE, nrow(bin_data))
    }
  }else{ # Otherwise, look for bins with a higher depth than median for a more stable estimate
    depth_selection <- (bin_data$bin.depth > median(bin_data$bin.depth, na.rm = T)) & bin_selection
  }

  selector <- depth_selection & bin_selection
  selector[is.na(selector)] <- FALSE

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
load.preprocess.bed <- function(bedpe_file, bin_data, tag_overlap){
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
    bed[, Chr := paste0('chr', Chr)]
  }

  bed <- bed[grepl('(chr[0-9]+|X|Y)$', Chr)] # Remove all decoy and alt contigs

  # Count doublets prior to trimming the Tn5 overlap
  bed[, doublet := is.doublet(bed, min.tag.overlap = tag_overlap-1, max.tag.overlap = tag_overlap)]
  bed[, End := End - tag_overlap]
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
  data.table::setorder(bed, Chr, Start, End, Strandedness)
  tmp.start <- paste(bed$Chr, bed$Start, bed$Strandedness)
  tmp.end <- paste(bed$Chr, bed$End, bed$Strandedness)

  are.duplicated <- duplicated(tmp.start, fromLast = TRUE) | duplicated(tmp.end)
  bed <- bed[!are.duplicated]
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
is.doublet <- function(bed, min.tag.overlap, max.tag.overlap){
  read.starts <- GenomicRanges::GRanges(seqnames = bed$Chr,
                                        ranges = IRanges::IRanges(start = bed$Start+min.tag.overlap, end = bed$Start+max.tag.overlap))
  read.ends <- GenomicRanges::GRanges(seqnames = bed$Chr,
                                      ranges = IRanges::IRanges(start = bed$End, end = bed$End))

  return(GenomicRanges::countOverlaps(read.starts, read.ends)>0)
}

#' count.overlaps
#'
#' Counts reads falling in the background (upstream) and doublet (overlap)
#' windows for each fragment.
#'
#' Two background window designs are supported:
#'
#' \describe{
#'   \item{\code{"original"}}{The window immediately upstream of each fragment:
#'     \code{[Start - max.size, Start - min.size - 1]}.
#'     Width = \code{max.size - min.size - 1} bp.
#'     Default. Fully backward-compatible.}
#'   \item{\code{"symmetric"}}{A window centred on \code{Start} that extends
#'     \code{symmetric_x} bp in both directions, with the doublet zone
#'     \code{[Start - min.size, Start - 1]} excised.  This produces two arms:
#'     left = \code{[Start - symmetric_x, Start - min.size - 1]} and
#'     right = \code{[Start, Start + symmetric_x]}.
#'     Width = combined arm widths (varies per fragment near chromosome ends).
#'     Avoids local accessibility biases in the narrow upstream window.}
#' }
#'
#' @param bed           data.table — fragment-level data, output of
#'   \code{load.preprocess.bed}.
#' @param min.size      integer — minimum fragment length (bp).
#' @param max.size      integer — maximum fragment length (bp).
#' @param background_method character — \code{"original"} (default) or
#'   \code{"symmetric"}.
#' @param symmetric_x   integer — half-width of the symmetric window in bp.
#'   Only used when \code{background_method = "symmetric"}. Default 25000.
#'
#' @return The input \code{bed} data.table with four new columns appended:
#'   \code{Count.Upstream}, \code{Count.Over},
#'   \code{Norm.Count.Upstream} (per-bp density), \code{Norm.Count.Over}.
count.overlaps <- function(bed, min.size, max.size,
                           background_method = c('original', 'symmetric'),
                           symmetric_x = 25000) {

  background_method <- match.arg(background_method)

  # ---- query ranges: fragments to count (start-point queries) -------------
  bed.sized     <- bed[Length > min.size & Length < max.size]
  for.counting  <- GenomicRanges::GRanges(
    seqnames = bed.sized$Chr,
    ranges   = IRanges::IRanges(start = bed.sized$Start, end = bed.sized$Start)
  )

  # ---- doublet (overlap) window — identical in both methods ---------------
  over.counting <- GenomicRanges::GRanges(
    seqnames = bed$Chr,
    ranges   = IRanges::IRanges(
      start = bed$Start - min.size,
      end   = bed$Start - 1L
    )
  )
  counted.over.overlaps <- GenomicRanges::countOverlaps(over.counting, for.counting)

  # ---- background window --------------------------------------------------
  if (background_method == 'original') {

    # Classic narrow upstream window: [Start - max.size, Start - min.size - 1]
    upstream_start <- pmax(1L, bed$Start - max.size)
    upstream_end   <- bed$Start - min.size - 1L

    # Drop any windows pushed off the chromosome start
    valid <- upstream_end >= upstream_start
    upstream_start[!valid] <- 1L
    upstream_end[!valid]   <- 1L   # will get 0 counts, handled below

    upstream.counting <- GenomicRanges::GRanges(
      seqnames = bed$Chr,
      ranges   = IRanges::IRanges(start = upstream_start, end = upstream_end)
    )
    counted.upstream.overlaps <- GenomicRanges::countOverlaps(
      upstream.counting, for.counting
    )
    counted.upstream.overlaps[!valid] <- NA_integer_

    window_size <- max.size - min.size   # constant, same as original

  } else {
    # Symmetric window centred on Start, doublet zone excised.
    # Left arm:  [Start - symmetric_x, Start - min.size - 1]
    # Right arm: [Start,               Start + symmetric_x ]
    # The doublet zone [Start - min.size, Start - 1] is intentionally excluded.

    # Exclusion zone: [Start - max.size, Start + max.size]
    # max.size is the 99th percentile of fragment lengths, so excluding +/-max.size
    # around the reference start ensures no read from the same molecule or the
    # same local accessibility peak contaminates the background estimate.
    left_start <- pmax(1L, bed$Start - symmetric_x)
    left_end   <- bed$Start - max.size - 1L
    left_valid <- left_end >= left_start

    right_start <- bed$Start + max.size
    right_end   <- bed$Start + symmetric_x

    # Left arm counts
    left_start_safe        <- left_start
    left_end_safe          <- left_end
    left_start_safe[!left_valid] <- 1L
    left_end_safe[!left_valid]   <- 1L

    left.counting <- GenomicRanges::GRanges(
      seqnames = bed$Chr,
      ranges   = IRanges::IRanges(start = left_start_safe, end = left_end_safe)
    )
    left_counts <- GenomicRanges::countOverlaps(left.counting, for.counting)
    left_counts[!left_valid] <- 0L

    # Right arm counts
    right.counting <- GenomicRanges::GRanges(
      seqnames = bed$Chr,
      ranges   = IRanges::IRanges(start = right_start, end = right_end)
    )
    right_valid <- right_start <= right_end
    right_counts <- GenomicRanges::countOverlaps(right.counting, for.counting)
    right_counts[!right_valid] <- 0L

    counted.upstream.overlaps <- left_counts + right_counts

    # Per-fragment window size (varies near chrom starts due to left arm clipping,
    # and if symmetric_x <= max.size the right arm collapses to zero)
    left_width  <- pmax(0L, left_end - left_start + 1L)
    left_width[!left_valid] <- 0L
    right_width <- pmax(0L, right_end - right_start + 1L)
    right_width[!right_valid] <- 0L
    window_size <- left_width + right_width
    window_size <- pmax(1L, window_size)   # avoid /0
  }

  # Assign directly onto bed with data.table := to preserve the data.table
  # class. cbind(bed, data.frame(...)) silently drops the data.table class,
  # which breaks the := and by= operations in estimate.ploidy that follow.
  bed[, Count.Upstream      := counted.upstream.overlaps]
  bed[, Count.Over          := counted.over.overlaps]
  bed[, Norm.Count.Upstream := counted.upstream.overlaps / window_size]
  bed[, Norm.Count.Over     := counted.over.overlaps     / min.size]
  bed[, Window.Size          := window_size]
  return(bed)
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
