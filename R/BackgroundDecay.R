#' plot_background_decay
#'
#' Diagnoses the local chromatin accessibility decay around Tn5 insertion sites
#' by computing the read density at successive distance bands from a sample of
#' reference fragment start positions.
#'
#' The density profile reveals three zones:
#' \enumerate{
#'   \item \strong{Doublet zone} (0 to ~\code{min_length} bp): extremely high
#'     density from reads on the same or adjacent molecules. This is what the
#'     doublet signal exploits.
#'   \item \strong{Decay zone} (\code{min_length} to ~\code{plateau_start} bp):
#'     elevated density due to local chromatin accessibility. Using this zone as
#'     background overestimates background density, causing the doublet ratio to
#'     be underestimated.
#'   \item \strong{Plateau zone} (> \code{plateau_start} bp): flat background
#'     density representing true local read depth, independent of the reference
#'     read's chromatin context. The background window should sit here.
#' }
#'
#' The vertical dashed line marks the suggested \code{symmetric_x} — the
#' distance at which density has decayed to within \code{plateau_tol} of the
#' far-field mean.
#'
#' @param bedpe_file character — path to a single \code{.bedpe.gz} file.
#' @param bin_data   data.frame — QDNAseq bin annotations (from
#'   \code{process.batch} internals), used to assign fragments to bins.
#' @param min_length integer — minimum fragment length (same as in
#'   \code{process.batch}). Default 50.
#' @param max_length integer — maximum fragment length. Default 1000.
#' @param tag_overlap integer — Tn5 footprint size in bp. Default 10.
#' @param max_dist    integer — maximum distance from reference start to
#'   profile. Default 200000 bp.
#' @param band_width  integer — width of each distance band in bp. Default 500.
#' @param n_sample    integer — number of reference reads to sample. Default
#'   5000. More gives a smoother profile but takes longer.
#' @param plateau_tol numeric — fraction of far-field mean within which density
#'   is considered flat. Default 0.05 (5%).
#' @param return_data logical — if TRUE, return the profile data.frame instead
#'   of plotting. Default FALSE.
#'
#' @return A ggplot object (printed as side-effect) or a data.frame if
#'   \code{return_data = TRUE}, invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' # Run on one cell's bedpe file to find optimal symmetric_x
#' plot_background_decay('path/to/cell.bedpe.gz',
#'                       bin_data   = QDNAseq::getBinAnnotations(500, 'hg38')@data,
#'                       max_dist   = 100000,
#'                       band_width = 1000)
#' }
plot_background_decay <- function(bedpe_file,
                                  bin_data,
                                  min_length   = 50L,
                                  max_length   = 1000L,
                                  tag_overlap  = 10L,
                                  max_dist     = 200000L,
                                  band_width   = 500L,
                                  n_sample     = 5000L,
                                  plateau_tol  = 0.05,
                                  return_data  = FALSE) {

  # ---- load and preprocess bed ----------------------------------------------
  bedpe <- tryCatch(
    data.table::fread(bedpe_file, sep = '\t'),
    error = function(e) stop("Cannot read bedpe file: ", conditionMessage(e))
  )
  colnames(bedpe) <- c("Chr1","Start1","End1","Chr2","Start2","End2",
                        "Name","Score","R1_direction","R2_direction")

  bed <- bedpe[, .(Chr   = Chr1,
                   Start = pmin(Start1, Start2),
                   End   = pmax(End1,   End2))]

  if (!any(grepl('^chr', bed$Chr))) bed[, Chr := paste0('chr', Chr)]
  bed <- bed[grepl('(chr[0-9]+|chrX|chrY)$', Chr)]
  bed[, End    := End - as.integer(tag_overlap)]
  bed[, Length := End - Start]
  bed <- bed[Length > min_length & Length < max_length]

  if (nrow(bed) == 0) stop("No fragments remaining after length filtering.")

  # ---- build GRanges of all fragment starts (for counting) -----------------
  all_starts <- GenomicRanges::GRanges(
    seqnames = bed$Chr,
    ranges   = IRanges::IRanges(start = bed$Start, end = bed$Start)
  )

  # ---- sample reference reads ----------------------------------------------
  n_sample  <- min(n_sample, nrow(bed))
  ref_idx   <- sample(nrow(bed), n_sample)
  ref_bed   <- bed[ref_idx]

  # ---- compute density in successive distance bands ------------------------
  # For band i covering distance [d_near, d_far] from the reference start:
  #   Left arm:  [Start - d_far, Start - d_near - 1]   (upstream)
  #   Right arm: [Start + d_near + 1, Start + d_far]   (downstream)
  # Both arms are 1 bp outside the reference read at d_near = 0.
  bands    <- seq(0L, as.integer(max_dist), by = as.integer(band_width))
  n_bands  <- length(bands) - 1L
  mid_dist <- (bands[-length(bands)] + bands[-1L]) / 2

  message("Profiling background density across ", n_bands,
          " distance bands using ", n_sample, " reference reads...")

  density <- vapply(seq_len(n_bands), function(i) {
    d_near <- as.integer(bands[i])
    d_far  <- as.integer(bands[i + 1L]) - 1L   # inclusive far edge

    # Left arm: [Start - d_far, Start - d_near - 1]
    l_end   <- ref_bed$Start - d_near - 1L
    l_start <- pmax(1L, ref_bed$Start - d_far)
    valid_left <- l_end >= l_start & l_end >= 1L

    if (any(valid_left)) {
      l_start[!valid_left] <- 1L
      l_end[!valid_left]   <- 1L
      left.gr <- GenomicRanges::GRanges(
        seqnames = ref_bed$Chr,
        ranges   = IRanges::IRanges(start = l_start, end = l_end)
      )
      left_counts <- GenomicRanges::countOverlaps(left.gr, all_starts)
      left_counts[!valid_left] <- 0L
    } else {
      left_counts <- rep(0L, nrow(ref_bed))
    }

    # Right arm: [Start + d_near + 1, Start + d_far]
    r_start <- ref_bed$Start + d_near + 1L
    r_end   <- ref_bed$Start + d_far
    right.gr <- GenomicRanges::GRanges(
      seqnames = ref_bed$Chr,
      ranges   = IRanges::IRanges(start = r_start, end = r_end)
    )
    right_counts <- GenomicRanges::countOverlaps(right.gr, all_starts)

    mean(left_counts + right_counts, na.rm = TRUE) / (2 * band_width)
  }, numeric(1))

  # ---- find plateau --------------------------------------------------------
  # Far-field mean: average of the outer 20% of the profile
  far_field  <- density[mid_dist > max_dist * 0.8]
  ff_mean    <- mean(far_field, na.rm = TRUE)
  ff_upper   <- ff_mean * (1 + plateau_tol)

  # Plateau start: first band (from far field inward) where density exceeds
  # ff_upper — the suggested symmetric_x sits just beyond this
  plateau_band <- max(which(density > ff_upper), na.rm = TRUE)
  suggested_x  <- bands[plateau_band + 1L]

  message(sprintf(
    "Far-field mean density: %.4f reads/bp | Suggested symmetric_x: %d bp (within %.0f%% of far-field)",
    ff_mean, suggested_x, plateau_tol * 100
  ))

  # ---- data frame ----------------------------------------------------------
  prof <- data.frame(
    mid_dist = mid_dist,
    density  = density,
    zone     = ifelse(mid_dist <= min_length,    'doublet',
               ifelse(mid_dist <= suggested_x,   'decay',
                                                  'plateau'))
  )

  if (return_data) return(invisible(prof))

  # ---- plot ----------------------------------------------------------------
  zone_colors <- c('doublet' = '#b3402e', 'decay' = '#e9c47e', 'plateau' = '#9fbdd7')

  p <- ggplot2::ggplot(prof, ggplot2::aes(x = mid_dist / 1000, y = density,
                                          colour = zone)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = ff_mean,    linetype = 'dashed',
                        colour = 'grey40', linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = ff_upper,   linetype = 'dotted',
                        colour = 'grey60', linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = suggested_x / 1000, linetype = 'dashed',
                        colour = '#b3402e', linewidth = 0.6) +
    ggplot2::annotate('text',
                      x     = suggested_x / 1000,
                      y     = max(density) * 0.95,
                      label = paste0('suggested\nsymmetric_x\n', suggested_x, ' bp'),
                      hjust = -0.1, size = 3, colour = '#b3402e') +
    ggplot2::scale_colour_manual(values = zone_colors, name = 'Zone') +
    ggplot2::scale_x_continuous(name = 'Distance from reference read start (kb)') +
    ggplot2::scale_y_continuous(name = 'Read density (reads / bp)') +
    ggplot2::labs(
      title    = 'Background read density vs distance from reference read',
      subtitle = paste0('Dashed grey = far-field mean  |  ',
                        'Red dashed = suggested symmetric_x  |  ',
                        'n = ', n_sample, ' reference reads')
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 8, colour = 'grey40'))

  plot(p)
  invisible(prof)
}


#' optimise_symmetric_x
#'
#' Finds the optimal \code{symmetric_x} background window size by profiling
#' the background read density as a function of distance from reference read
#' starts (via \code{\link{plot_background_decay}}) and returning the distance
#' at which density has decayed to within \code{plateau_tol} of the far-field
#' mean.
#'
#' This is the same logic as \code{\link{plot_background_decay}} but returns
#' the numeric value directly for use in \code{\link{process.batch}}.
#'
#' @inheritParams plot_background_decay
#'
#' @return integer — suggested \code{symmetric_x} in bp.
#' @export
#'
#' @examples
#' \dontrun{
#' opt_x <- optimise_symmetric_x('path/to/cell.bedpe.gz',
#'                                bin_data   = bins@data,
#'                                max_dist   = 100000)
#' sbird_sce <- process.batch(bams, bedpes, genome = 'hg38',
#'                            background_method = 'symmetric',
#'                            symmetric_x       = opt_x)
#' }
optimise_symmetric_x <- function(bedpe_file,
                                  bin_data,
                                  min_length  = 50L,
                                  max_length  = 1000L,
                                  tag_overlap = 10L,
                                  max_dist    = 200000L,
                                  band_width  = 500L,
                                  n_sample    = 5000L,
                                  plateau_tol = 0.05) {

  prof <- plot_background_decay(
    bedpe_file  = bedpe_file,
    bin_data    = bin_data,
    min_length  = min_length,
    max_length  = max_length,
    tag_overlap = tag_overlap,
    max_dist    = max_dist,
    band_width  = band_width,
    n_sample    = n_sample,
    plateau_tol = plateau_tol,
    return_data = TRUE
  )

  far_field   <- prof$density[prof$mid_dist > max_dist * 0.8]
  ff_mean     <- mean(far_field, na.rm = TRUE)
  ff_upper    <- ff_mean * (1 + plateau_tol)
  plateau_band <- max(which(prof$density > ff_upper), na.rm = TRUE)
  suggested_x  <- as.integer(prof$mid_dist[plateau_band] + (band_width / 2))

  message("Optimal symmetric_x: ", suggested_x, " bp")
  suggested_x
}
