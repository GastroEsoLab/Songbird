
#' Data Heatmap Plotter
#'
#' @param sce
#' @param assay_name
#' @param cell_attribs
#' @param row_split
#' @param rowClust
#' @param use_raster
#' @param bin_attribs
#' @param return
#' @param row_title
#' @param plot_title
#' @param legend_title
#'
#' @return
#' @export
#'
#' @examples
plot_heatmap <- function(sce, assay_name, cell_attribs = NULL, row_split = NULL,
                         rowClust = T, use_raster = TRUE, bin_attribs = NULL, return = FALSE,
                         row_title = NULL, plot_title = NULL, legend_title = NULL){

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
  if(assay_name == 'copy'){
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
               show_row_names = F,
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

#' Title
#'
#' @param values
#' @param orientation
#' @param ylim
#'
#' @return
#'
#' @examples
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


#' Title
#'
#' @param sce
#'
#' @return
#'
#' @examples
get_binMetadata <- function(sce){
  return(as.data.frame(sce@rowRanges@elementMetadata))
}



#' Title
#'
#' @param sce
#' @param cell_id
#' @param assay
#' @param return
#'
#' @return
#' @export
#'
#' @examples
plot_cell <- function(sce, cell_id, assay, return = F){

  data <- assay(sce, entry)
  data <- data[, cell_id]

  bin_data <- get_binMetadata(sce)
  data <- data.frame(expr = data, chrom = bin_data$chr)

  spacer <- data.frame(expr = rep(NA, 50), chrom = 'spacer')
  while(i < nrow(data)-1){
    if(data$chrom[i] != data$chrom[i+1]){
      data <- rbind(data[1:i,], spacer, data[(i+1):nrow(data),])
      print(data$chrom[i:(i+12)])
      i <- i + nrow(spacer)
    }
    i <- i + 1
  }
  data$pos <- seq(1, nrow(data))

  chr_pos <- c()
  for(chr in unique(data$chrom)){
    avg_position <- mean(data$pos[data$chrom == chr])
    chr_pos <- rbind(chr_pos, data.frame(pos = avg_position, chrom = chr))
  }

  p <- ggplot2::ggplot(data, ggplot2::aes(x = pos, y = expr)) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = 'Chromosome', y = names(assay), title = cell) +
    ggplot2::scale_x_continuous(breaks = chr_pos$pos, labels = chr_pos$chrom)

  if(return){
    return(p)
  } else {
    plot(p)
  }
}

#' Title
#'
#' @param sbird_sce
#' @param cell
#' @param chr
#' @param return_plot
#'
#' @return
#' @export
#'
#' @examples
plot_cell <- function(sbird_sce, cell, chr = NULL, return_plot = FALSE){
  copy_mtx <- SummarizedExperiment::assay(sbird_sce, 'copy')
  reads_mtx <- SummarizedExperiment::assay(sbird_sce, 'reads')

  dat <- data.frame(copy = copy_mtx[,colnames(copy_mtx)==cell,drop=T],
                    reads = reads_mtx[,colnames(reads_mtx)==cell,drop=T],
                    bin_number = seq(1, nrow(reads_mtx)),
                    bin_name = rownames(reads_mtx))
  dat$outlier <- dat$copy > quantile(dat$copy, 0.95, na.rm = T)

  # Set plotting max based on outliers
  ymax <- max(10,max(dat$copy[!dat$outlier], na.rm = T))

  # Scale the reads to fit on the plot grid
  dat$deviation <- dat$copy/dat$reads
  deviations <- dat$deviation
  deviations <- deviations[is.finite(deviations)]
  dat$adj_reads <- dat$reads*mean(deviations)

  colors <- c('#496bab', '#9fbdd7', '#c1c1c1', '#e9c47e',
              '#d6804f', '#b3402e', '#821010', '#6a0936',
              '#ab1964', '#b6519f', '#ad80b9', '#c2a9d1')
  names(colors) <- c(0:11)

  # Get the plotting position for the chromsomes
  dat$chr <- gsub('^([0-9]+|X|Y)_.*', '\\1', dat$bin_name)
  chr_labels <- unique(dat$chr)
  chr_locs <- sapply(chrs, function(x) mean(which(dat$chr == x)))
  border_locs <- sapply(chrs, function(x) max(which(dat$chr == x)))
  border_locs <- border_locs[1:(length(border_locs)-1)]
  all_locs <- c(chr_locs, border_locs)
  all_labels <- c(chr_labels, rep('', length(border_locs)))
  all_tics <- c(rep(NA, length(chr_locs)), rep('black', length(border_locs)))

  if(!is.null(chr)){dat <- dat[dat$chr==chr,]}
  p <- ggplot2::ggplot(dat, aes(x = bin_number, y = adj_reads, color = as.factor(copy))) + ggplot2::geom_point(size = 0.3) +
    ggplot2::geom_point(aes(y = copy)) +
    ggplot2::scale_color_manual(name = 'Copy\nNumber', values = colors) +
    ggplot2::scale_y_continuous(name = 'Copy Number', breaks = seq(0, ymax, 1), labels = seq(0, ymax, 1), limits = c(0, ymax)) +
    ggplot2::scale_x_continuous(name = 'Chromosome', breaks = all_locs, labels = all_labels) + theme_classic() +
    ggplot2::theme(axis.ticks.x = element_line(color = all_tics))
  #ggplot2::scale_x_continuous(name = 'Chromosome', breaks = chr_locs, labels = chrs) + theme_classic() +
  #ggplot2::theme(axis.ticks.x = element_line(color = c(rep(NA, length(border_ticks)))))


  #scale_x_continuous(breaks = c(sort(unique(data$x)), x_tick),
  #                   labels = c(sort(unique(data$name)), rep(c(""), len))) +
  #  theme(axis.ticks.x = element_line(color = c(rep(NA, len - 1), rep("black", len))))

  if(return_plot){return(p)}
  else{plot(p)}
}

