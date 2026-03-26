#' plot_overlap_fit
#'
#' Plots observed overlap read counts against expected counts derived from the
#' background density, for a single cell. This is the key diagnostic for the
#' doublet-based ploidy model: under a diploid genome the points should fall
#' on a line with slope 0.5. Aneuploid regions scatter above or below
#' according to their copy number state.
#'
#' Expected count per bin:
#' \deqn{\hat{y}_b = \frac{\text{Norm.Count.Upstream}_b}{\text{min\_length}} \times
#'   \text{Bin.Reads}_b \times \frac{k-1}{k}}
#'
#' @param sbird_sce  SingleCellExperiment — songbird object after
#'   \code{\link{process.batch}} with bedpe files.
#' @param cell       character — cell name (column name in assay matrices).
#' @param min_length integer — overlap window size in bp. If NULL, read from
#'   \code{metadata(sbird_sce)$overlap_window_size}, fallback 50.
#' @param cn_states  integer vector — CN states for reference lines. Default 1:6.
#' @param alpha      numeric — point transparency. Default 0.4.
#' @param return_plot logical — return ggplot object instead of printing.
#'
#' @return A ggplot object (if \code{return_plot = TRUE}) or NULL.
#' @export
plot_overlap_fit <- function(sbird_sce, cell,
                             min_length  = NULL,
                             cn_states   = 1:6,
                             alpha       = 0.4,
                             return_plot = FALSE) {

  # Resolve min_length
  if (is.null(min_length)) {
    min_length <- S4Vectors::metadata(sbird_sce)$overlap_window_size
    if (is.null(min_length)) min_length <- S4Vectors::metadata(sbird_sce)$min_length
    if (is.null(min_length)) { message("min_length not in metadata, defaulting to 50"); min_length <- 50L }
  }

  assay_names  <- SummarizedExperiment::assayNames(sbird_sce)
  required     <- c('Norm.Count.Upstream', 'Count.Over', 'Bin.Reads')
  missing_assays <- setdiff(required, assay_names)
  if (length(missing_assays) > 0) {
    stop("plot_overlap_fit requires assays: ", paste(missing_assays, collapse = ', '),
         "\nThese are present only when process.batch was run with bedpe files.")
  }

  all_cells <- colnames(SummarizedExperiment::assay(sbird_sce, 'reads'))
  if (!cell %in% all_cells) {
    stop("Cell '", cell, "' not found.")
  }
  cell_idx <- which(all_cells == cell)

  norm_up   <- SummarizedExperiment::assay(sbird_sce, 'Norm.Count.Upstream')[, cell_idx]
  count_ov  <- SummarizedExperiment::assay(sbird_sce, 'Count.Over')[, cell_idx]
  bin_reads <- SummarizedExperiment::assay(sbird_sce, 'Bin.Reads')[, cell_idx]
  rd        <- as.data.frame(SummarizedExperiment::rowData(sbird_sce))

  dat <- data.frame(
    expected   = norm_up * min_length * bin_reads,
    observed   = count_ov,
    use        = if ('overlap_use' %in% names(rd)) rd$overlap_use else TRUE,
    stringsAsFactors = FALSE
  )
  dat <- dat[!is.na(dat$expected) & dat$expected > 0 & !is.na(dat$observed) & dat$use, ]
  if (nrow(dat) == 0) stop("No usable bins found for cell '", cell, "'.")

  ref_lines <- data.frame(
    cn    = cn_states,
    slope = (cn_states - 1) / cn_states
  )

  x_max <- stats::quantile(dat$expected, 0.99, na.rm = TRUE)
  y_max <- stats::quantile(dat$observed, 0.99, na.rm = TRUE)

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = expected, y = observed)) +
    ggplot2::geom_point(size = 0.8, alpha = alpha, colour = 'grey40') +
    ggplot2::geom_abline(data = ref_lines,
                         mapping = ggplot2::aes(slope = slope, intercept = 0,
                                                colour = factor(cn)),
                         linewidth = 0.7, linetype = 'dashed') +
    ggplot2::scale_colour_manual(
      name   = 'CN state',
      values = setNames(c('#9fbdd7','#c1c1c1','#e9c47e','#d6804f','#b3402e','#821010',
                          '#6a0936','#ab1964','#b6519f','#ad80b9','#c2a9d1')[
                            seq_along(cn_states)], as.character(cn_states)),
      drop = TRUE) +
    ggplot2::coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max)) +
    ggplot2::labs(x     = 'Expected overlap count',
                  y     = 'Observed overlap count',
                  title = paste0('Overlap fit — ', cell),
                  subtitle = 'Dashed lines show theoretical slope (k-1)/k per CN state') +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 9, colour = 'grey40'))

  if (return_plot) return(p) else plot(p)
}


#' plot_hmm_doublet
#'
#' Plots the observed doublet fraction per bin (Count.Over / expected) alongside
#' the HMM-decoded copy number segmentation for a single cell. Each point is
#' the empirical doublet fraction for one bin; a running mean is overlaid; the
#' HMM CN segmentation is shown as coloured horizontal line segments at the
#' theoretical (k-1)/k level for each segment.
#'
#' The doublet fraction lives in [0,1): diploid ~0.5, triploid ~0.667,
#' tetraploid ~0.75, etc.
#'
#' @param sbird_sce  SingleCellExperiment — songbird object after
#'   \code{\link{ploidy_correction}} or copy number calling.
#'   Requires assays \code{"hmm_cn"} (if available), \code{"Norm.Count.Upstream"},
#'   \code{"Count.Over"}, and \code{"Bin.Reads"}.
#' @param cell       character — cell name.
#' @param window     integer — bins for running mean smoother. Default 11.
#'   Must be odd; if even, incremented by 1.
#' @param chr        character — chromosome to plot (e.g. "1" or "chr1").
#'   NULL plots all chromosomes.
#' @param return_plot logical — return ggplot object instead of printing.
#'
#' @return A ggplot object (if \code{return_plot = TRUE}) or NULL.
#' @export
plot_hmm_doublet <- function(sbird_sce, cell,
                             window      = 11L,
                             chr         = NULL,
                             return_plot = FALSE) {

  if (window %% 2 == 0) { window <- window + 1L; message('window incremented to ', window) }

  assay_names <- SummarizedExperiment::assayNames(sbird_sce)
  required    <- c('Norm.Count.Upstream', 'Count.Over', 'Bin.Reads')
  missing_a   <- setdiff(required, assay_names)
  if (length(missing_a) > 0)
    stop("plot_hmm_doublet requires assays: ", paste(missing_a, collapse = ', '))

  all_cells <- colnames(SummarizedExperiment::assay(sbird_sce, 'reads'))
  if (!cell %in% all_cells) stop("Cell '", cell, "' not found.")
  cell_idx <- which(all_cells == cell)

  ows <- S4Vectors::metadata(sbird_sce)$overlap_window_size
  if (is.null(ows)) ows <- S4Vectors::metadata(sbird_sce)$min_length
  if (is.null(ows)) ows <- 50L

  norm_up   <- SummarizedExperiment::assay(sbird_sce, 'Norm.Count.Upstream')[, cell_idx]
  count_ov  <- SummarizedExperiment::assay(sbird_sce, 'Count.Over')[, cell_idx]
  bin_reads <- SummarizedExperiment::assay(sbird_sce, 'Bin.Reads')[, cell_idx]
  rd        <- as.data.frame(SummarizedExperiment::rowData(sbird_sce))

  expected  <- norm_up * ows * bin_reads
  obs_ratio <- ifelse(is.finite(expected) & expected > 0, count_ov / expected, NA_real_)

  # Running mean
  n        <- length(obs_ratio)
  half     <- window %/% 2L
  smoothed <- vapply(seq_len(n), function(i)
    mean(obs_ratio[max(1L, i - half):min(n, i + half)], na.rm = TRUE), numeric(1))

  # HMM CN if available
  has_hmm <- 'hmm_cn' %in% assay_names
  if (has_hmm) {
    hmm_cn <- SummarizedExperiment::assay(sbird_sce, 'hmm_cn')[, cell_idx]
  } else {
    hmm_cn <- rep(NA_integer_, n)
  }

  dat <- data.frame(
    bin_number = seq_len(n),
    chr        = as.character(rd$chr),
    obs_ratio  = obs_ratio,
    smoothed   = smoothed,
    hmm_cn     = hmm_cn,
    stringsAsFactors = FALSE
  )

  if (!is.null(chr)) {
    chr_clean <- sub('^chr', '', chr)
    dat <- dat[sub('^chr', '', dat$chr) == chr_clean, ]
  }
  dat$bin_number <- seq_len(nrow(dat))

  # Chromosome axis labels
  chr_labels  <- unique(dat$chr)
  chr_locs    <- sapply(chr_labels, function(x) mean(which(dat$chr == x)))
  border_locs <- sapply(chr_labels, function(x) max(which(dat$chr == x)))
  border_locs <- border_locs[seq_len(length(border_locs) - 1)]
  all_locs    <- c(chr_locs, border_locs)
  all_labels  <- c(chr_labels, rep('', length(border_locs)))
  all_tics    <- c(rep(NA, length(chr_locs)), rep('black', length(border_locs)))

  cn_palette        <- c('#496bab','#9fbdd7','#c1c1c1','#e9c47e','#d6804f','#b3402e',
                         '#821010','#6a0936','#ab1964','#b6519f','#ad80b9','#c2a9d1')
  names(cn_palette) <- as.character(0:11)
  dat$cn_factor     <- factor(as.character(pmin(dat$hmm_cn, 11L)), levels = as.character(0:11))

  ymax <- max(1, stats::quantile(dat$obs_ratio, 0.99, na.rm = TRUE) * 1.05)

  # HMM segments via RLE on non-NA bins
  valid_bins <- which(!is.na(dat$hmm_cn))
  if (length(valid_bins) > 0) {
    rle_cn      <- rle(dat$hmm_cn[valid_bins])
    seg_end_v   <- cumsum(rle_cn$lengths)
    seg_start_v <- c(1L, seg_end_v[-length(seg_end_v)] + 1L)
    seg_df <- data.frame(
      x_start  = dat$bin_number[valid_bins[seg_start_v]],
      x_end    = dat$bin_number[valid_bins[seg_end_v]],
      y        = (rle_cn$values - 1) / rle_cn$values,
      cn_factor = factor(as.character(pmin(rle_cn$values, 11L)), levels = as.character(0:11)),
      stringsAsFactors = FALSE
    )
  } else {
    seg_df <- data.frame(x_start = integer(), x_end = integer(), y = numeric(),
                         cn_factor = factor(character(), levels = as.character(0:11)))
  }

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = bin_number)) +
    ggplot2::geom_point(ggplot2::aes(y = obs_ratio, colour = cn_factor),
                        size = 0.4, alpha = 0.5, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = smoothed),
                       colour = 'black', linewidth = 0.4, na.rm = TRUE) +
    ggplot2::geom_segment(data = seg_df,
                          ggplot2::aes(x = x_start, xend = x_end,
                                       y = y, yend = y, colour = cn_factor),
                          linewidth = 1.2, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 'dashed',
                        colour = 'grey60', linewidth = 0.4) +
    ggplot2::scale_colour_manual(name = 'HMM CN', values = cn_palette, drop = TRUE) +
    ggplot2::scale_y_continuous(
      name   = 'Doublet fraction (Count.Over / expected)',
      limits = c(0, ymax),
      breaks = seq(0, 1, 0.1)) +
    ggplot2::scale_x_continuous(name = 'Chromosome',
                                 breaks = all_locs, labels = all_labels) +
    ggplot2::labs(
      title    = paste0('Doublet signal \u2014 ', cell),
      subtitle = paste0('Points: obs/exp ratio | Line: running mean (window=', window,
                        ') | Segments: HMM (k-1)/k | Dashed: diploid baseline (0.5)')) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.ticks.x  = ggplot2::element_line(color = all_tics),
                   plot.subtitle = ggplot2::element_text(size = 8, colour = 'grey40'))

  if (return_plot) return(p) else plot(p)
}


#' inspect_cell
#'
#' Plots three panels for a single cell side by side for diagnostic inspection:
#' (1) read depth with copy number overlay (\code{\link{plot_cell}}),
#' (2) the doublet signal (\code{\link{plot_hmm_doublet}}),
#' (3) the overlap fit (\code{\link{plot_overlap_fit}}).
#'
#' Requires the \pkg{cowplot} package. If \pkg{cowplot} is not installed, the
#' three plots are printed individually.
#'
#' @param sbird_sce SingleCellExperiment — songbird object.
#' @param cell      character — cell name.
#' @param chr       character — chromosome to plot, or NULL for all.
#' @param assay     character — assay for \code{plot_cell}. Default 'copy'.
#' @param window    integer — running mean window for doublet plot. Default 11.
#'
#' @return invisibly returns a list of the three ggplot objects.
#' @export
inspect_cell <- function(sbird_sce, cell, chr = NULL, assay = 'copy', window = 11L) {

  p1 <- plot_cell(sbird_sce, cell, assay = assay, chr = chr, return_plot = TRUE)

  has_doublet_assays <- all(c('Norm.Count.Upstream', 'Count.Over', 'Bin.Reads') %in%
                              SummarizedExperiment::assayNames(sbird_sce))
  if (has_doublet_assays) {
    p2 <- plot_hmm_doublet(sbird_sce, cell, chr = chr, window = window, return_plot = TRUE)
    p3 <- plot_overlap_fit(sbird_sce, cell, return_plot = TRUE)
    plots <- list(copy = p1, doublet = p2, overlap_fit = p3)
  } else {
    message("inspect_cell: doublet assays not found — plotting read depth only.")
    plots <- list(copy = p1)
  }

  if (requireNamespace('cowplot', quietly = TRUE) && length(plots) == 3) {
    print(cowplot::plot_grid(plotlist = plots, ncol = 1, align = 'v'))
  } else {
    for (p in plots) print(p)
  }
  invisible(plots)
}


#' plot_background_decay
#'
#' Diagnoses the local chromatin accessibility decay around Tn5 insertion sites
#' by computing the read density at successive distance bands from a sample of
#' reference fragment starts. Use this to choose \code{background_exclusion}
#' for \code{\link{process.batch}} with \code{background_method = "symmetric"}.
#'
#' Three zones are identified:
#' \enumerate{
#'   \item \strong{Doublet zone} (0 to \code{min_length} bp): very high density
#'     from same-molecule reads — this is what the doublet signal exploits.
#'   \item \strong{Decay zone}: elevated density from local chromatin accessibility.
#'     Using this region as background overestimates background density.
#'   \item \strong{Plateau zone}: flat far-field density — the correct background.
#'     Set \code{background_exclusion} to the start of this zone.
#' }
#'
#' @param bedpe_file character — path to a .bedpe or .bedpe.gz file.
#' @param min_length integer — minimum fragment length. Default 50.
#' @param max_length integer — maximum fragment length. Default 1000.
#' @param tag_overlap integer — Tn5 footprint size in bp. Default 10.
#' @param max_dist    integer — maximum distance to profile in bp. Default 150000.
#' @param band_width  integer — width of each distance band in bp. Default 1000.
#' @param n_sample    integer — reference reads to sample. Default 5000.
#' @param plateau_tol numeric — fraction of far-field mean within which density is
#'   considered flat. Default 0.05.
#' @param return_data logical — if TRUE return the data.frame instead of plotting.
#'
#' @return A ggplot object (printed) or a data.frame if \code{return_data = TRUE}.
#' @export
plot_background_decay <- function(bedpe_file,
                                   min_length   = 50L,
                                   max_length   = 1000L,
                                   tag_overlap  = 10L,
                                   max_dist     = 150000L,
                                   band_width   = 1000L,
                                   n_sample     = 5000L,
                                   plateau_tol  = 0.05,
                                   return_data  = FALSE) {

  # Load and preprocess fragments (same as estimate.ploidy)
  bedpe <- tryCatch(
    data.table::fread(bedpe_file, sep = '\t'),
    error = function(e) stop("Cannot read bedpe file: ", conditionMessage(e)))
  colnames(bedpe) <- c("Chr1","Start1","End1","Chr2","Start2","End2",
                        "Name","Score","R1_direction","R2_direction")

  bed <- bedpe[, .(Chr   = Chr1,
                   Start = pmin(Start1, Start2),
                   End   = pmax(End1, End2))]

  chr_prefix <- substr(bed$Chr, 1, 3)
  if (sum(grepl('chr', chr_prefix)) == 0) bed[, Chr := paste0('chr', Chr)]
  bed <- bed[grepl('(chr[0-9]+|chrX|chrY)$', Chr)]
  bed[, End    := End - as.integer(tag_overlap)]
  bed[, Length := End - Start]
  bed <- bed[Length > min_length & Length < max_length]

  if (nrow(bed) == 0) stop("No fragments remaining after length filtering.")

  all_starts <- GenomicRanges::GRanges(
    seqnames = bed$Chr,
    ranges   = IRanges::IRanges(start = bed$Start, end = bed$Start))

  n_sample  <- min(n_sample, nrow(bed))
  ref_idx   <- sample(nrow(bed), n_sample)
  ref_bed   <- bed[ref_idx]

  bands    <- seq(0L, as.integer(max_dist), by = as.integer(band_width))
  n_bands  <- length(bands) - 1L
  mid_dist <- (bands[-length(bands)] + bands[-1L]) / 2

  message("Profiling background density across ", n_bands,
          " distance bands using ", n_sample, " reference reads...")

  density <- vapply(seq_len(n_bands), function(i) {
    d_near <- as.integer(bands[i])
    d_far  <- as.integer(bands[i + 1L]) - 1L

    # Left arm: [Start - d_far, Start - d_near - 1]
    l_end   <- ref_bed$Start - d_near - 1L
    l_start <- pmax(1L, ref_bed$Start - d_far)
    valid_l <- l_end >= l_start & l_end >= 1L
    if (any(valid_l)) {
      ls <- l_start; le <- l_end
      ls[!valid_l] <- 1L; le[!valid_l] <- 1L
      left.gr     <- GenomicRanges::GRanges(bed$Chr[ref_idx],
                       IRanges::IRanges(start = ls, end = le))
      left_counts <- GenomicRanges::countOverlaps(left.gr, all_starts)
      left_counts[!valid_l] <- 0L
    } else {
      left_counts <- rep(0L, nrow(ref_bed))
    }

    # Right arm: [Start + d_near + 1, Start + d_far]
    r_start <- ref_bed$Start + d_near + 1L
    r_end   <- ref_bed$Start + d_far
    right.gr     <- GenomicRanges::GRanges(bed$Chr[ref_idx],
                      IRanges::IRanges(start = r_start, end = r_end))
    right_counts <- GenomicRanges::countOverlaps(right.gr, all_starts)

    mean(left_counts + right_counts, na.rm = TRUE) / (2 * band_width)
  }, numeric(1))

  far_field    <- density[mid_dist > max_dist * 0.8]
  ff_mean      <- mean(far_field, na.rm = TRUE)
  ff_upper     <- ff_mean * (1 + plateau_tol)
  above        <- which(density > ff_upper)
  plateau_band <- if (length(above) > 0) max(above) else 1L
  suggested_x  <- bands[plateau_band + 1L]

  message(sprintf(
    "Far-field mean density: %.4f reads/bp | Suggested background_exclusion: %d bp (within %.0f%% of far-field)",
    ff_mean, suggested_x, plateau_tol * 100))

  prof <- data.frame(
    mid_dist = mid_dist,
    density  = density,
    zone     = ifelse(mid_dist <= min_length,   'doublet',
               ifelse(mid_dist <= suggested_x,  'decay', 'plateau'))
  )

  if (return_data) return(invisible(prof))

  zone_colors <- c('doublet' = '#b3402e', 'decay' = '#e9c47e', 'plateau' = '#9fbdd7')

  p <- ggplot2::ggplot(prof, ggplot2::aes(x = mid_dist / 1000, y = density,
                                           colour = zone)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = ff_mean,  linetype = 'dashed',
                        colour = 'grey40', linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = ff_upper, linetype = 'dotted',
                        colour = 'grey60', linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = suggested_x / 1000, linetype = 'dashed',
                        colour = '#b3402e', linewidth = 0.6) +
    ggplot2::annotate('text', x = suggested_x / 1000, y = max(density) * 0.95,
                      label = paste0('suggested\nbackground_exclusion\n', suggested_x, ' bp'),
                      hjust = -0.1, size = 3, colour = '#b3402e') +
    ggplot2::scale_colour_manual(values = zone_colors, name = 'Zone') +
    ggplot2::scale_x_continuous(name = 'Distance from reference read start (kb)') +
    ggplot2::scale_y_continuous(name = 'Read density (reads / bp)') +
    ggplot2::labs(
      title    = 'Background read density vs distance from reference read',
      subtitle = paste0('Grey dashed = far-field mean | Red dashed = suggested background_exclusion | ',
                        'n = ', n_sample, ' reference reads')) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 8, colour = 'grey40'))

  plot(p)
  invisible(prof)
}
