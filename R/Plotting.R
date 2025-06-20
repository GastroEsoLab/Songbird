
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
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = ggplot2::.data$bin_number, y = ggplot2::.data$adj_reads, color = as.factor(ggplot2::.data$copy))) + ggplot2::geom_point(size = 0.3) +
      ggplot2::geom_point(ggplot2::aes(y = ggplot2::.data$copy))
  }else{
    ymax <- stats::quantile(dat$copy, 0.99, na.rm = T)
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = ggplot2::.data$bin_number, y = ggplot2::.data$adj_reads)) + ggplot2::geom_point(size = 0.3, color = 'grey') +
      ggplot2::geom_point(ggplot2::aes(y = ggplot2::.data$copy))
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

