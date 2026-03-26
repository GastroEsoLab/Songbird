
#' Data Heatmap Plotter
#'
#' @param sce songbird object
#' @param assay_name the name of the assay to plot
#' @param cell_attribs column name of the metadata you want plotted next to the cells
#' @param row_split vector of identities to group the cells - typically subclonal id
#' @param rowClust boolean indicating whether to cluster the rows of the heatmap
#' @param use_raster boolean indicating whether to use raster graphics for the heatmap
#' @param bin_attribs column name of the column attributes to plot above the heatmap
#' @param return boolean indicating whether to return the heatmap object or just plot it out
#' @param row_title title for the row labels of the heatmap
#' @param plot_title title for the heatmap
#' @param legend_title title for the heatmap legend
#' @param show_cellNames boolean indicating whether to show the cell names on the heatmap
#'
#' @return a heatmap object if return = TRUE, otherwise plots the heatmap
#' @export
#'
#' @examples
plot_heatmap <- function(sce, assay_name, cell_attribs = NULL, row_split = NULL,
                         rowClust = T, use_raster = TRUE, bin_attribs = NULL, return = FALSE,
                         row_title = NULL, plot_title = NULL, legend_title = NULL,
                         show_cellNames = FALSE){

  # Extract the matrix and split by chromosomes
  matrix <- t(SummarizedExperiment::assay(sce, assay_name))
  matrix[is.na(matrix)] <- 0
  chr_names <- SingleCellExperiment::rowData(sce)$chr

  #Put a new line after odd chromosome names and before even chromosome names
  chr_names <- ifelse(grepl('(1$|3$|5$|7$|9$|11$|13$|15$|17$|19$|21$|X)', chr_names), paste0(chr_names, '\n'), paste0('\n', chr_names))
  chr_names <- factor(chr_names, levels = c("1\n", "\n2", "3\n", "\n4", "5\n",
                                            "\n6", "7\n", "\n8", "9\n", "\n10",
                                            "11\n", "\n12", "13\n", "\n14", "15\n",
                                            "\n16", "17\n", "\n18", "19\n", "\n20",
                                            "21\n", "\n22", "X\n", "\nY" ))

  # Set colormap
  if(grepl('copy', assay_name)){
    colScale <- circlize::colorRamp2(breaks = 0:11,
                                     colors = c('#496bab', '#9fbdd7', '#c1c1c1', '#e9c47e',
                                                '#d6804f', '#b3402e', '#821010', '#6a0936',
                                                '#ab1964', '#b6519f', '#ad80b9', '#c2a9d1'))
    matrix[matrix > 12] = 12
  }else{
    minMax <- stats::quantile(matrix, c(0.01, .99), na.rm = T)
    colScale <- circlize::colorRamp2(seq(from = minMax[1], to = minMax[2], length.out = 1000),
                                     scico::scico(1000, palette = 'batlow'))
  }

  #Create cell annotations
  if(!is.null(cell_attribs)){
    # Check that cell_attribs is a vector
    metadata <- as.data.frame(sce@colData)

    if(!is.vector(cell_attribs)){
      stop('cell_attribs must be a named vector of metadata columns')
    }

    if(!all(cell_attribs %in% names(metadata))){
      stop('cell_attribs must be a named vector of metadata columns')
    }

    if(is.null(names(cell_attribs))){
      names(cell_attribs) <- cell_attribs
    }

    lt <- lapply(metadata[,cell_attribs,drop=F], function(x) gen_annobar(x, 'row'))
    names(lt) <- names(cell_attribs)
    cell_anno <- do.call(ComplexHeatmap::rowAnnotation, lt)

    for(i in 1:length(cell_anno@anno_list)){
      cell_anno@anno_list[[i]]@name_param$rot <- 45
    }
  }else{
    cell_anno <- cell_attribs
  }

  # Repeat for the bin attribs
  if(!is.null(bin_attribs)){
    metadata <- get_binMetadata(sce)
    metadata[is.na(metadata)] <- 0
    if(!is.vector(bin_attribs)){
      stop('bin_anno must be a named vector of metadata columns')
    }

    if(!all(bin_attribs %in% names(metadata))){
      stop('bin_anno must be a named vector of metadata columns')
    }

    if(is.null(names(bin_attribs))){
      names(bin_attribs) <- bin_attribs
    }

    lt <- lapply(metadata[,bin_attribs,drop=F], function(x) gen_annobar(x, 'column'))
    names(lt) <- names(bin_attribs)
    bin_anno <- do.call(ComplexHeatmap::HeatmapAnnotation, lt)
  }else{
    bin_anno <- bin_attribs
  }


  # Plot the heatmap
  p <- ComplexHeatmap::Heatmap(matrix,
               cluster_columns = F,
               cluster_rows = T,
               show_row_names = show_cellNames,
               show_column_names = F,
               column_split = chr_names,
               show_row_dend = F,
               row_split = row_split,
               left_annotation = cell_anno,
               top_annotation = bin_anno,
               col = colScale,
               use_raster = use_raster,
               row_title = row_title,
               column_title = plot_title,
               name = legend_title)

  if(return){
    return(p)
  }else{
    plot(p)
  }
}

#' gen_annobar
#'
#' @param values a vector of values to plot as an annotation bar
#' @param orientation a string indicating the orientation of the annotation bar, either 'row' or 'column'
#' @param ylim a numeric vector of length 2 indicating the limits of the y-axis for the barplot. If NULL, the limits will be set to the 1st and 99th percentiles of the values.
#'
#' @return a ComplexHeatmap annotation object
gen_annobar <- function(values, orientation, ylim = NULL){
  if(!is.numeric(values)){
    return(ComplexHeatmap::anno_simple(values, which = orientation))
  }else{
    if(is.null(ylim)){
      ylim <- minMax <- stats::quantile(values, c(0.01, .99), na.rm = T)
    }
    values[values > ylim[2]] <- ylim[2]
    values[values < ylim[1]] <- ylim[1]
    return(ComplexHeatmap::anno_barplot(values, ylim = ylim, which = orientation))
  }
}


#' get_binMetadata
#'
#' @param sce songbird object
#'
#' @return column metadata for the songbird object
get_binMetadata <- function(sce){
  return(as.data.frame(sce@rowRanges@elementMetadata))
}

#' plot_cell
#'
#' @param sbird_sce songbird object
#' @param cell cell name to plot
#' @param assay assay to plot. Will be overlaid over reads
#' @param chr chromosome to plot, if NULL all chromosomes will be plotted
#' @param return_plot if TRUE, returns the plot object instead of plotting it
#'
#' @return a ggplot object or plots the data
#' @export
#'
#' @examples
plot_cell <- function(sbird_sce, cell, assay = 'copy', chr = NULL, return_plot = FALSE){
  copy_mtx <- SummarizedExperiment::assay(sbird_sce, assay)
  reads_mtx <- SummarizedExperiment::assay(sbird_sce, 'reads')

  dat <- data.frame(copy = copy_mtx[,colnames(copy_mtx)==cell,drop=T],
                    reads = reads_mtx[,colnames(reads_mtx)==cell,drop=T],
                    bin_number = seq(1, nrow(reads_mtx)),
                    bin_name = rownames(reads_mtx))
  dat$outlier <- dat$copy > stats::quantile(dat$copy, 0.95, na.rm = T)
  # Scale the reads to fit on the plot grid
  deviation <- mean(dat$copy, na.rm = TRUE) / mean(dat$reads, na.rm = TRUE)
  dat$adj_reads <- dat$reads*deviation

  # Set plotting max based on outliers
  if(assay == 'copy'){
    ymax <- max(10,max(dat$copy[!dat$outlier], na.rm = T))

    colors <- c('#496bab', '#9fbdd7', '#c1c1c1', '#e9c47e',
                '#d6804f', '#b3402e', '#821010', '#6a0936',
                '#ab1964', '#b6519f', '#ad80b9', '#c2a9d1')
    names(colors) <- c(0:11)
    dat$copy[dat$copy>11] <- 11
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = bin_number, y = adj_reads, color = as.factor(copy))) + ggplot2::geom_point(size = 0.3) +
      ggplot2::geom_point(ggplot2::aes(y = copy))
  }else{
    ymax <- stats::quantile(dat$copy, 0.99, na.rm = T)
    colors <- NULL
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = bin_number, y = adj_reads)) + ggplot2::geom_point(size = 0.3, color = 'grey') +
      ggplot2::geom_point(ggplot2::aes(y = copy))
  }

  # Get the plotting position for the chromsomes
  dat$chr <- gsub('^([0-9]+|X|Y)_.*', '\\1', dat$bin_name)
  chr_labels <- unique(dat$chr)
  chr_locs <- sapply(chr_labels, function(x) mean(which(dat$chr == x)))
  border_locs <- sapply(chr_labels, function(x) max(which(dat$chr == x)))
  border_locs <- border_locs[1:(length(border_locs)-1)]
  all_locs <- c(chr_locs, border_locs)
  all_labels <- c(chr_labels, rep('', length(border_locs)))
  all_tics <- c(rep(NA, length(chr_locs)), rep('black', length(border_locs)))

  if(!is.null(chr)){dat <- dat[dat$chr==chr,]}
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_color_manual(name = 'Copy\nNumber', values = colors)
  }
  p <- p +
    ggplot2::scale_y_continuous(name = 'Copy Number', breaks = seq(0, ymax, 1), labels = seq(0, ymax, 1), limits = c(0, ymax)) +
    ggplot2::scale_x_continuous(name = 'Chromosome', breaks = all_locs, labels = all_labels) + ggplot2::theme_classic() +
    ggplot2::theme(axis.ticks.x = ggplot2::element_line(color = all_tics))

  if(return_plot){return(p)}
  else{plot(p)}
}


# =============================================================================
# NEW DIAGNOSTIC FUNCTIONS
# =============================================================================

#' inspect_cell
#'
#' Side-by-side comparison of copy number calls and wavelet segmentation
#' for a single cell.
#'
#' @param sbird_sce songbird object (after \code{\link{copyCall}})
#' @param cell      character — cell name
#' @param chr       character — chromosome to plot, or NULL for all
#' @param return_plot logical — return a named list of ggplot objects
#'
#' @return invisibly NULL or named list of ggplot objects
#' @export
#'
#' @examples
#' \dontrun{
#' inspect_cell(sbird_sce, cell = colnames(sbird_sce)[1])
#' }
inspect_cell <- function(sbird_sce, cell, chr = NULL, return_plot = FALSE) {
  p_copy <- plot_cell(sbird_sce, cell, assay = 'copy',      chr = chr, return_plot = TRUE)
  p_seg  <- plot_cell(sbird_sce, cell, assay = 'segmented', chr = chr, return_plot = TRUE)

  p_copy <- p_copy + ggplot2::ggtitle(paste0(cell, '  \u2014  copy number'))
  p_seg  <- p_seg  + ggplot2::ggtitle(paste0(cell, '  \u2014  wavelet segmentation'))

  if (return_plot) return(list(copy = p_copy, segmented = p_seg))
  gridExtra::grid.arrange(p_copy, p_seg, ncol = 1)
  invisible(NULL)
}


#' plot_overlap_fit
#'
#' Diagnostic scatter plot of observed vs expected doublet overlap counts.
#' Expected = Norm.Count.Upstream * min_length * Bin.Reads.  Reference lines
#' show the theoretical slope \code{(k-1)/k} for each CN state.  The solid
#' black line is the observed mean slope across all bins.
#'
#' Requires bedpe-derived columns in \code{rowData}
#' (Norm.Count.Upstream, Norm.Count.Over, Bin.Reads).
#'
#' @param sbird_sce songbird object
#' @param cell      character — cell name
#' @param min_length integer — overlap window size in bp. Read from
#'   \code{metadata(sbird_sce)$overlap_window_size} if NULL.
#' @param cn_states integer vector — CN states for reference lines. Default 1:6.
#' @param alpha     numeric — point transparency. Default 0.4.
#' @param return_plot logical — return ggplot object. Default FALSE.
#'
#' @return ggplot object or NULL
#' @export
#'
#' @examples
#' \dontrun{
#' plot_overlap_fit(sbird_sce, cell = colnames(sbird_sce)[1])
#' }
plot_overlap_fit <- function(sbird_sce, cell,
                              min_length  = NULL,
                              cn_states   = 1:6,
                              alpha       = 0.4,
                              return_plot = FALSE) {

  if (is.null(min_length)) {
    min_length <- S4Vectors::metadata(sbird_sce)$overlap_window_size
    if (is.null(min_length)) min_length <- S4Vectors::metadata(sbird_sce)$min_length
    if (is.null(min_length)) { message("min_length not found in metadata, defaulting to 50"); min_length <- 50L }
  }

  rd <- as.data.frame(SummarizedExperiment::rowData(sbird_sce))
  required_rd <- c('Norm.Count.Upstream', 'Norm.Count.Over', 'Bin.Reads', 'use', 'mappability', 'bases')
  missing_rd  <- setdiff(required_rd, colnames(rd))
  if (length(missing_rd) > 0)
    stop("plot_overlap_fit requires rowData columns: ", paste(missing_rd, collapse = ', '),
         "\nThese are present only when process.batch was run with bedpe files.")

  all_cells <- colnames(SummarizedExperiment::assay(sbird_sce, 'reads'))
  if (!cell %in% all_cells) stop("Cell '", cell, "' not found.")

  dat <- data.frame(
    Norm.Count.Upstream = rd$Norm.Count.Upstream,
    Norm.Count.Over     = rd$Norm.Count.Over,
    Bin.Reads           = rd$Bin.Reads,
    use                 = rd$use,
    mappability         = rd$mappability,
    bases               = rd$bases,
    stringsAsFactors    = FALSE
  )
  dat$expected <- dat$Norm.Count.Upstream * min_length * dat$Bin.Reads

  dat <- dat[!is.na(dat$expected) & dat$expected > 0 &
               !is.na(dat$Norm.Count.Over) & dat$use &
               !is.na(dat$mappability) & dat$mappability >= 90 &
               !is.na(dat$bases)       & dat$bases       >= 90, ]
  if (nrow(dat) == 0) stop("No usable bins found for cell '", cell, "'.")

  ref_lines <- data.frame(cn = cn_states, slope = (cn_states - 1) / cn_states)
  obs_slope <- sum(dat$Norm.Count.Over, na.rm = TRUE) / sum(dat$expected, na.rm = TRUE)

  x_max <- stats::quantile(dat$expected,       0.99, na.rm = TRUE)
  y_max <- stats::quantile(dat$Norm.Count.Over, 0.99, na.rm = TRUE)

  cn_palette <- c('0'='#496bab','1'='#9fbdd7','2'='#c1c1c1','3'='#e9c47e',
                  '4'='#d6804f','5'='#b3402e','6'='#821010','7'='#6a0936',
                  '8'='#ab1964','9'='#b6519f','10'='#ad80b9','11'='#c2a9d1')

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = expected, y = Norm.Count.Over)) +
    ggplot2::geom_point(size = 0.8, alpha = alpha, colour = 'grey50') +
    ggplot2::geom_abline(data = ref_lines,
      mapping = ggplot2::aes(slope = slope, intercept = 0, colour = as.character(cn)),
      linewidth = 0.7, linetype = 'dashed', show.legend = TRUE) +
    ggplot2::geom_abline(slope = obs_slope, intercept = 0,
      colour = 'black', linewidth = 0.9, linetype = 'solid') +
    ggplot2::annotate('text', x = x_max * 0.02, y = y_max * 0.97,
      label = sprintf('obs slope = %.3f', obs_slope), hjust = 0, size = 3.5) +
    ggplot2::scale_colour_manual(name = 'CN state',
      values = cn_palette[as.character(cn_states)],
      labels = paste0('CN', cn_states, '  (', round((cn_states-1)/cn_states, 3), ')')) +
    ggplot2::coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max)) +
    ggplot2::labs(
      x = paste0('Expected overlap count  (Norm.Count.Upstream \u00d7 ', min_length, ' bp)'),
      y = 'Observed Norm.Count.Over',
      title = paste0('Overlap fit  \u2014  ', cell),
      subtitle = 'Dashed = theoretical slope (k-1)/k  |  Solid = observed mean slope') +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 9, colour = 'grey40'))

  if (return_plot) return(p) else plot(p)
}


#' plot_hmm_doublet
#'
#' Plots the observed doublet ratio (Norm.Count.Over / Norm.Count.Upstream)
#' per genomic bin with a running mean overlay.  The y-axis is the empirical
#' doublet fraction in [0,1): diploid ~0.5, triploid ~0.667, tetraploid ~0.75.
#'
#' Requires bedpe-derived columns in \code{rowData}.
#'
#' @param sbird_sce songbird object
#' @param cell      character — cell name
#' @param window    integer — running mean window width in bins (odd). Default 5.
#' @param chr       character — chromosome to plot, or NULL for all
#' @param return_plot logical — return ggplot object. Default FALSE.
#'
#' @return ggplot object or NULL
#' @export
#'
#' @examples
#' \dontrun{
#' plot_hmm_doublet(sbird_sce, cell = colnames(sbird_sce)[1], window = 41)
#' }
plot_hmm_doublet <- function(sbird_sce, cell,
                              window      = 5L,
                              chr         = NULL,
                              return_plot = FALSE) {

  rd        <- as.data.frame(SummarizedExperiment::rowData(sbird_sce))
  all_cells <- colnames(SummarizedExperiment::assay(sbird_sce, 'reads'))

  required_rd <- c('Norm.Count.Upstream', 'Norm.Count.Over')
  missing_rd  <- setdiff(required_rd, colnames(rd))
  if (length(missing_rd) > 0)
    stop("plot_hmm_doublet requires rowData columns: ", paste(missing_rd, collapse = ', '),
         "\nThese are present only when process.batch was run with bedpe files.")

  if (!cell %in% all_cells) stop("Cell '", cell, "' not found.")

  if (window %% 2 == 0) { window <- window + 1L; message('window incremented to ', window) }

  dat <- data.frame(
    bin_number = seq_len(nrow(rd)),
    bin_name   = rd$bin_name,
    chr        = rd$chr,
    obs_ratio  = ifelse(rd$Norm.Count.Upstream > 0,
                        rd$Norm.Count.Over / rd$Norm.Count.Upstream, NA_real_),
    stringsAsFactors = FALSE
  )

  if (!is.null(chr)) dat <- dat[dat$chr == chr, ]
  dat$bin_number <- seq_len(nrow(dat))

  half <- window %/% 2L; n <- nrow(dat)
  dat$smoothed <- vapply(seq_len(n), function(i) {
    lo <- max(1L, i - half); hi <- min(n, i + half)
    mean(dat$obs_ratio[lo:hi], na.rm = TRUE)
  }, numeric(1))

  chr_labels  <- unique(dat$chr)
  chr_locs    <- sapply(chr_labels, function(x) mean(which(dat$chr == x)))
  border_locs <- sapply(chr_labels, function(x) max(which(dat$chr == x)))
  border_locs <- border_locs[seq_len(length(border_locs) - 1)]
  all_locs    <- c(chr_locs, border_locs)
  all_labels  <- c(chr_labels, rep('', length(border_locs)))
  all_tics    <- c(rep(NA, length(chr_locs)), rep('black', length(border_locs)))

  ymax <- max(1, stats::quantile(dat$obs_ratio, 0.99, na.rm = TRUE) * 1.05)

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = bin_number)) +
    ggplot2::geom_point(ggplot2::aes(y = obs_ratio),
      size = 0.4, alpha = 0.4, colour = 'grey60', na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = smoothed),
      colour = 'black', linewidth = 0.5, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 'dashed',
      colour = 'grey40', linewidth = 0.4) +
    ggplot2::scale_y_continuous(
      name = 'Doublet fraction  (Norm.Count.Over / Norm.Count.Upstream)',
      limits = c(0, ymax), breaks = seq(0, 1, 0.1)) +
    ggplot2::scale_x_continuous(name = 'Chromosome',
      breaks = all_locs, labels = all_labels) +
    ggplot2::labs(
      title = paste0('Doublet signal  \u2014  ', cell),
      subtitle = paste0('Points: obs/exp ratio per bin  |  ',
        'Line: running mean (window = ', window, ')  |  ',
        'Dashed: diploid baseline (0.5)')) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.ticks.x = ggplot2::element_line(color = all_tics),
                   plot.subtitle = ggplot2::element_text(size = 8, colour = 'grey40'))

  if (return_plot) return(p) else plot(p)
}


#' plot_background_decay
#'
#' Profiles the read density at successive distance bands from sampled
#' fragment starts to identify the local chromatin accessibility decay and
#' suggest an optimal \code{symmetric_x} / \code{background_exclusion}.
#'
#' @param bedpe_file character — path to a bedpe or bedpe.gz file
#' @param min_length integer — minimum fragment length. Default 50.
#' @param max_length integer — maximum fragment length. Default 1000.
#' @param tag_overlap integer — Tn5 footprint in bp. Default 10.
#' @param max_dist    integer — maximum distance to profile in bp. Default 150000.
#' @param band_width  integer — width of each distance band in bp. Default 1000.
#' @param n_sample    integer — number of reference reads to sample. Default 5000.
#' @param plateau_tol numeric — fraction of far-field mean for plateau detection. Default 0.05.
#' @param return_data logical — return profile data.frame instead of plotting. Default FALSE.
#'
#' @return ggplot (side-effect) or data.frame (invisibly) if \code{return_data = TRUE}
#' @export
#'
#' @examples
#' \dontrun{
#' plot_background_decay('path/to/cell.bedpe.gz')
#' }
plot_background_decay <- function(bedpe_file,
                                   min_length   = 50L,
                                   max_length   = 1000L,
                                   tag_overlap  = 10L,
                                   max_dist     = 150000L,
                                   band_width   = 1000L,
                                   n_sample     = 5000L,
                                   plateau_tol  = 0.05,
                                   return_data  = FALSE) {

  bedpe <- tryCatch(data.table::fread(bedpe_file, sep = '\t'),
    error = function(e) stop("Cannot read bedpe file: ", conditionMessage(e)))

  colnames(bedpe) <- c("Chr1","Start1","End1","Chr2","Start2","End2",
                       "Name","Score","R1_direction","R2_direction")
  bed <- bedpe[, .(Chr = Chr1, Start = pmin(Start1, Start2), End = pmax(End1, End2))]
  if (!any(grepl('^chr', bed$Chr))) bed[, Chr := paste0('chr', Chr)]
  bed <- bed[grepl('(chr[0-9]+|chrX|chrY)$', Chr)]
  bed[, End    := End - as.integer(tag_overlap)]
  bed[, Length := End - Start]
  bed <- bed[Length > min_length & Length < max_length]
  if (nrow(bed) == 0) stop("No fragments remaining after length filtering.")

  all_starts <- GenomicRanges::GRanges(seqnames = bed$Chr,
    ranges = IRanges::IRanges(start = bed$Start, end = bed$Start))

  n_sample <- min(n_sample, nrow(bed))
  ref_bed  <- bed[sample(nrow(bed), n_sample)]

  bands    <- seq(0L, as.integer(max_dist), by = as.integer(band_width))
  n_bands  <- length(bands) - 1L
  mid_dist <- (bands[-length(bands)] + bands[-1L]) / 2

  message("Profiling ", n_bands, " distance bands using ", n_sample, " reference reads...")

  density <- vapply(seq_len(n_bands), function(i) {
    d_near <- as.integer(bands[i]); d_far <- as.integer(bands[i + 1L]) - 1L

    l_end <- ref_bed$Start - d_near - 1L; l_start <- pmax(1L, ref_bed$Start - d_far)
    l_valid <- l_end >= l_start & l_end >= 1L
    if (any(l_valid)) {
      ls <- l_start; le <- l_end; ls[!l_valid] <- 1L; le[!l_valid] <- 1L
      lc <- GenomicRanges::countOverlaps(
        GenomicRanges::GRanges(seqnames = ref_bed$Chr,
                               ranges = IRanges::IRanges(start = ls, end = le)),
        all_starts)
      lc[!l_valid] <- 0L
    } else { lc <- rep(0L, nrow(ref_bed)) }

    rc <- GenomicRanges::countOverlaps(
      GenomicRanges::GRanges(seqnames = ref_bed$Chr,
        ranges = IRanges::IRanges(start = ref_bed$Start + d_near + 1L,
                                  end   = ref_bed$Start + d_far)),
      all_starts)
    mean(lc + rc, na.rm = TRUE) / (2 * band_width)
  }, numeric(1))

  far_field  <- density[mid_dist > max_dist * 0.8]
  ff_mean    <- mean(far_field, na.rm = TRUE)
  ff_upper   <- ff_mean * (1 + plateau_tol)
  exceed_idx <- which(density > ff_upper)
  plat_band  <- if (length(exceed_idx) > 0) max(exceed_idx) else 1L
  suggested_x <- as.integer(mid_dist[plat_band] + band_width / 2)

  message(sprintf("Far-field mean: %.5f reads/bp | Suggested symmetric_x: %d bp",
                  ff_mean, suggested_x))

  prof <- data.frame(mid_dist = mid_dist, density = density,
    zone = ifelse(mid_dist <= min_length, 'doublet',
           ifelse(mid_dist <= suggested_x, 'decay', 'plateau')))

  if (return_data) return(invisible(prof))

  zone_colors <- c('doublet'='#b3402e', 'decay'='#e9c47e', 'plateau'='#9fbdd7')

  p <- ggplot2::ggplot(prof, ggplot2::aes(x = mid_dist/1000, y = density, colour = zone)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = ff_mean,  linetype='dashed',  colour='grey40', linewidth=0.4) +
    ggplot2::geom_hline(yintercept = ff_upper, linetype='dotted',  colour='grey60', linewidth=0.3) +
    ggplot2::geom_vline(xintercept = suggested_x/1000, linetype='dashed', colour='#b3402e', linewidth=0.6) +
    ggplot2::annotate('text', x=suggested_x/1000, y=max(density)*0.95,
      label=paste0('suggested\n', suggested_x, ' bp'), hjust=-0.1, size=3, colour='#b3402e') +
    ggplot2::scale_colour_manual(values = zone_colors, name = 'Zone') +
    ggplot2::scale_x_continuous(name = 'Distance from reference read start (kb)') +
    ggplot2::scale_y_continuous(name = 'Read density (reads / bp)') +
    ggplot2::labs(title = 'Background read density vs distance from reference read',
      subtitle = paste0('Grey dashed = far-field mean  |  Red dashed = suggested symmetric_x  |  ',
                        'n = ', n_sample, ' reference reads')) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 8, colour = 'grey40'))

  plot(p)
  invisible(prof)
}
