
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
  if(grepl('copy', assay_name) || assay_name == 'hmm_cn'){
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
      names(bin_anno) <- bin_attribs
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
  p <- p + ggplot2::scale_color_manual(name = 'Copy\nNumber', values = colors) +
    ggplot2::scale_y_continuous(name = 'Copy Number', breaks = seq(0, ymax, 1), labels = seq(0, ymax, 1), limits = c(0, ymax)) +
    ggplot2::scale_x_continuous(name = 'Chromosome', breaks = all_locs, labels = all_labels) + ggplot2::theme_classic() +
    ggplot2::theme(axis.ticks.x = ggplot2::element_line(color = all_tics))

  if(return_plot){return(p)}
  else{plot(p)}
}


#' plot_overlap_fit
#'
#' Plots observed overlap read counts against expected counts derived from the
#' background (upstream) density, for a single cell.  This is the key
#' diagnostic for the Poisson doublet-based ploidy model: under a pure
#' diploid genome the points should fall on the line
#' \code{observed = ratio_diploid x expected}, where
#' \code{ratio_diploid = 0.5}.  Aneuploid regions scatter above or below
#' according to their copy number state.
#'
#' When the SCE contains the \code{"hmm_cn"} assay (i.e. after
#' \code{\link{ploidy_correction_hmm}} has been run), each bin is coloured by
#' its decoded CN state and reference lines are drawn for every CN state
#' present in the cell.  Without \code{"hmm_cn"}, a single diploid reference
#' line is shown.
#'
#' Expected count per bin is computed as:
#' \deqn{\hat{y}_b = \text{ratio}_k \times
#'   \frac{\text{Count.Upstream}_b}{\text{Total.Window.Size}_b} \times
#'   \text{min\_length} \times \text{Bin.Reads}_b}
#' where \eqn{\text{ratio}_k = (k-1)/k} for CN state \eqn{k}.
#'
#' @param sbird_sce  SingleCellExperiment — songbird object after
#'   \code{\link{process.batch}} (with bedpe files).  Requires
#'   assays \code{Count.Over}, \code{Count.Upstream}, \code{Total.Window.Size},
#'   and \code{Bin.Reads}.
#' @param cell       character — cell name (column name in the assay matrices).
#' @param min_length integer — overlap window size in bp.  Should match the
#'   \code{min_length} used in \code{\link{process.batch}}.  If
#'   \code{NULL} (default), read from \code{metadata(sbird_sce)$min_length};
#'   falls back to 50 if not found.
#' @param cn_states  integer vector — CN states for which reference lines are
#'   drawn.  Default \code{1:6}.
#' @param alpha      numeric — point transparency.  Default 0.4.
#' @param return_plot logical — if \code{TRUE} returns the ggplot object
#'   instead of printing.  Default \code{FALSE}.
#'
#' @return A ggplot object (if \code{return_plot = TRUE}) or \code{NULL}
#'   (printed as side-effect).
#' @export
#'
#' @examples
#' \dontrun{
#' # Before HMM — single diploid reference line
#' plot_overlap_fit(sbird_sce, cell = 'cell_001')
#'
#' # After HMM — bins coloured by CN state, one line per state
#' sbird_sce <- ploidy_correction_hmm(sbird_sce)
#' plot_overlap_fit(sbird_sce, cell = 'cell_001')
#' }
plot_overlap_fit <- function(sbird_sce, cell,
                             min_length  = NULL,
                             cn_states   = 1:6,
                             alpha       = 0.4,
                             return_plot = FALSE) {

  # ---- resolve min_length --------------------------------------------------
  if (is.null(min_length)) {
    min_length <- S4Vectors::metadata(sbird_sce)$min_length
    if (is.null(min_length)) {
      message("min_length not found in metadata, defaulting to 50")
      min_length <- 50L
    }
  }

  # ---- check required assays -----------------------------------------------
  assay_names  <- SummarizedExperiment::assayNames(sbird_sce)
  overlap_cols <- c('Count.Over', 'Count.Upstream', 'Total.Window.Size', 'Bin.Reads')
  missing      <- setdiff(overlap_cols, assay_names)
  if (length(missing) > 0) {
    stop(
      "plot_overlap_fit requires assays: ",
      paste(missing, collapse = ', '), "\n",
      "These are present only when process.batch was run with bedpe files."
    )
  }

  rd        <- as.data.frame(SummarizedExperiment::rowData(sbird_sce))
  all_cells <- colnames(SummarizedExperiment::assay(sbird_sce, 'reads'))

  # ---- check cell exists ---------------------------------------------------
  if (!cell %in% all_cells) {
    stop("Cell '", cell, "' not found. Available cells:\n",
         paste(head(all_cells, 10), collapse = ', '),
         if (length(all_cells) > 10) paste0(' ... (', length(all_cells), ' total)') else '')
  }

  # ---- extract per-cell overlap values from assay matrices -----------------
  cell_idx   <- which(all_cells == cell)
  count_over  <- SummarizedExperiment::assay(sbird_sce, 'Count.Over')[, cell_idx]
  count_up    <- SummarizedExperiment::assay(sbird_sce, 'Count.Upstream')[, cell_idx]
  win_size    <- SummarizedExperiment::assay(sbird_sce, 'Total.Window.Size')[, cell_idx]
  bin_reads   <- SummarizedExperiment::assay(sbird_sce, 'Bin.Reads')[, cell_idx]

  # ---- build per-bin data frame --------------------------------------------
  dat <- data.frame(
    bin_name      = rd$bin_name,
    chr           = rd$chr,
    Count.Over    = count_over,
    Count.Up      = count_up,
    Win.Size      = win_size,
    Bin.Reads     = bin_reads,
    use           = rd$overlap_use,
    stringsAsFactors = FALSE
  )

  # Expected count: background density x overlap window x reads in bin
  # Expected = (Count.Upstream / Total.Window.Size) x min_length
  # Uses raw background counts and mean window size to correctly scale across
  # different background window designs (original narrow vs wide symmetric).
  dat$expected <- (dat$Count.Up / dat$Win.Size) * min_length * dat$Bin.Reads

  # ---- CN state colouring --------------------------------------------------
  has_hmm <- 'hmm_cn' %in% SummarizedExperiment::assayNames(sbird_sce)

  if (has_hmm) {
    hmm_mat    <- SummarizedExperiment::assay(sbird_sce, 'hmm_cn')
    cell_cn    <- hmm_mat[, cell]
    dat$cn     <- factor(cell_cn, levels = sort(unique(stats::na.omit(cell_cn))))
    cn_present <- sort(unique(stats::na.omit(cell_cn)))
  } else {
    dat$cn     <- factor(2L)
    cn_present <- 2L
  }

  # ---- remove bins with no usable data ------------------------------------
  dat <- dat[!is.na(dat$expected) & dat$expected > 0 &
             !is.na(dat$Count.Over) & !is.na(dat$Win.Size) &
             dat$Win.Size > 0 & dat$use, ]

  if (nrow(dat) == 0) {
    stop("No usable bins found for cell '", cell,
         "' after applying quality filters.")
  }

  # ---- reference lines: one per CN state -----------------------------------
  # Slope = ratio_k = (k-1)/k.  Line passes through origin.
  ref_lines <- data.frame(
    cn    = factor(cn_present, levels = levels(dat$cn)),
    slope = (cn_present - 1) / cn_present
  )

  # ---- observed mean slope per CN state ------------------------------------
  # Fit Count.Over ~ expected through the origin (no intercept) per CN state.
  # slope = sum(observed) / sum(expected) — the MLE for a Poisson rate ratio,
  # equivalent to lm(y ~ 0 + x) but without the normality assumption.
  # This gives the actual inferred doublet ratio for each state, which should
  # be close to ratio_k = (k-1)/k if the model fits well.
  obs_slopes <- do.call(rbind, lapply(cn_present, function(k) {
    idx <- !is.na(dat$cn) & as.integer(as.character(dat$cn)) == k
    if (sum(idx) < 3) return(NULL)
    obs_slope <- sum(dat$Count.Over[idx], na.rm = TRUE) /
                 sum(dat$expected[idx],   na.rm = TRUE)
    data.frame(
      cn         = factor(k, levels = levels(dat$cn)),
      obs_slope  = obs_slope,
      label      = sprintf('CN%d: obs=%.3f  th=%.3f', k, obs_slope, (k-1)/k)
    )
  }))

  # CN colour palette — matches the copy number palette used in plot_cell
  cn_palette <- c(
    '0'  = '#496bab', '1' = '#9fbdd7', '2' = '#c1c1c1',
    '3'  = '#e9c47e', '4' = '#d6804f', '5' = '#b3402e',
    '6'  = '#821010', '7' = '#6a0936', '8' = '#ab1964',
    '9'  = '#b6519f', '10'= '#ad80b9', '11'= '#c2a9d1'
  )
  used_palette <- cn_palette[as.character(sort(unique(as.integer(
    as.character(dat$cn)
  ))))]

  # ---- plot ----------------------------------------------------------------
  # Clip extreme outliers for display (keep 99th percentile)
  x_max <- stats::quantile(dat$expected,    0.99, na.rm = TRUE)
  y_max <- stats::quantile(dat$Count.Over,  0.99, na.rm = TRUE)

  p <- ggplot2::ggplot(dat,
         ggplot2::aes(x = expected, y = Count.Over, colour = cn)) +
    ggplot2::geom_point(size = 0.8, alpha = alpha) +
    # Theoretical slope (dashed)
    ggplot2::geom_abline(
      data    = ref_lines,
      mapping = ggplot2::aes(slope = slope, intercept = 0, colour = cn),
      linewidth = 0.7, linetype = 'dashed', show.legend = FALSE
    ) +
    # Observed mean slope (solid) — passes through origin, fitted per CN state
    ggplot2::geom_abline(
      data    = obs_slopes,
      mapping = ggplot2::aes(slope = obs_slope, intercept = 0, colour = cn),
      linewidth = 0.9, linetype = 'solid', show.legend = FALSE
    ) +
    # Annotation: observed vs theoretical slope per CN state
    ggplot2::annotate(
      'text',
      x     = x_max * 0.02,
      y     = seq(y_max * 0.98, by = -y_max * 0.07, length.out = nrow(obs_slopes)),
      label = obs_slopes$label,
      hjust = 0, vjust = 1, size = 3, colour = 'grey30'
    ) +
    ggplot2::scale_colour_manual(
      name   = 'CN state',
      values = used_palette,
      drop   = TRUE
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, x_max),
      ylim = c(0, y_max)
    ) +
    ggplot2::labs(
      x     = 'Expected overlap count  (Count.Upstream / window size \u00d7 min_length)',
      y     = 'Observed overlap count',
      title = paste0('Overlap fit  \u2014  ', cell),
      subtitle = if (has_hmm)
        'Dashed = theoretical slope (k-1)/k; solid = observed mean slope; annotations show obs vs theoretical'
      else
        'Dashed = diploid reference (0.5); solid = observed slope; run ploidy_correction_hmm() for per-state view'
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = 9, colour = 'grey40')
    )

  if (return_plot) return(p) else plot(p)
}
