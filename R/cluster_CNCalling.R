#' Title
#'
#' @param data list of cell processing results
#' @param column column to select from
#' @param return_counts boolean to return total reads
#'
#' @return
#' @export
#'
#' @examples
#'
reads_to_matrix <- function(data, column, return_counts = FALSE) {
  all_bins <- c()
  total_reads <- c()
  for(cell in data){
    cell$binName <- paste0(cell$chromosome, "_", cell$start, "_", cell$end)
    all_bins <- c(all_bins, cell$binName)
    total_reads <- c(total_reads, sum(cell$uncorrected.reads, na.rm = T))
  }
  all_bins <- unique(all_bins)
  cell_names <- sapply(data, function(x) unique(x$cell_id))
  dat_matrix <- matrix(nrow = length(data), ncol = length(all_bins),
                       dimnames = list(cell_names, all_bins))
  for(i in 1:length(data)){
    dat <- unlist(data[[i]][column])
    col_index <- match(colnames(dat_matrix), paste0(data[[i]]$chromosome, "_", data[[i]]$start, "_", data[[i]]$end))
    dat_matrix[i,col_index] <- dat
  }
  if(return_counts){
    return(list(dat_matrix, total_reads))
  }else{
    return(dat_matrix)
  }
}

#' Title
#'
#' @param means output from the unbalanced haar transform
#' @param use boolean array marking high quality bins
#' @param expected_ploidy noisy estimator of true cell ploidy
#' @param sigma standard deviation of the normal distribution
#'
#' @return
#' @export
#'
#' @examples
fitMeans <- function(means, use, expected_ploidy, sigma = 0.5){
  # Abnormal ploidy states are considered diploid
  if(expected_ploidy < 0.5){
    expected_ploidy <- 2
  }else if(expected_ploidy > 15){
    expected_ploidy <- 15
  }

  final_cn <- rep(NA, length(means))
  all_means <- means
  use_idx <- which(use)
  means <- means[use]
  nBins <- length(means)

  # Identify largest peak and propose CN state
  dens <- stats::density(means)
  peak <- abs(dens$x[which.max(dens$y)])

  # main detected peak could be any one of 1,2,3,4,5 state - modeled by a poisson mean 2
  modeState <- seq(round(expected_ploidy*0.75), round(expected_ploidy*1.5))
  modeState <- modeState[modeState>0]

  fitScores <- c() # How well does the mu fit the fixed states applied to them (dnormal around each state)
  stateScores <- c() # How well does the inferred states fit our strong prior that cells should be diploid?
  fitStates <- matrix(nrow = length(modeState), ncol = nBins)

  for(i in 1:length(modeState)){
    # Given mode state - find uniploid state
    uniploid <- peak/modeState[i]

    # Get the total number of states (+1 to account for 0 copy regions) in the dataset
    # given uniploid state & score the gaussians for each
    numStates <- max(round(means/uniploid))+1
    if(numStates>1000){
      numStates <- 1000
    }
    scores <- matrix(nrow = numStates, ncol = nBins)
    for(j in 1:numStates){
      scores[j,] <- stats::dnorm(means,
                                 mean = uniploid*(j-1), # Account for 0 copy regions
                                 sd = sigma,
                                 log = T)
    }

    # Best fitting state is the normal that best fits each state
    states <- apply(scores, 2, which.max)
    topScores <- rep(0, length(states))
    for(j in 1:length(states)){
      topScores[j] <- scores[states[j],j]
    }
    fitScores <- c(fitScores, sum(topScores))
    stateScores <- c(stateScores, sum(stats::dnorm(mean = states-1, sd = sigma, expected_ploidy, log = T)))
    fitStates[i,] <- states - 1
  }

  # Get the highest scoring fit and apply it to the good bins
  stateScores <- sapply(stateScores, function(x) x-matrixStats::logSumExp(stateScores))
  fitScores <- sapply(fitScores, function(x) x-matrixStats::logSumExp(fitScores))
  bestFit <- which.max(stateScores + fitScores)
  final_cn[use_idx] <- fitStates[bestFit,]

  # get true uniploid & apply to the remaining bins
  best_uniploid <- peak/modeState[bestFit]
  final_cn[!use_idx] <- round(all_means[!use_idx]/best_uniploid)
  return(final_cn)
}

#' Title
#'
#' @param sbird_sce
#'
#' @return
#' @export
#'
#' @examples
identify_subclones <- function(sbird_sce){
  # Cluster the cells using the changepoint matrix and sce
  change_mtx <- generate_changepoint_matrix(SummarizedExperiment::assay(sbird_sce, 'segmented'), use_mask = SummarizedExperiment::rowData(sbird_sce)$overlap_use)
  SingleCellExperiment::reducedDim(sbird_sce, 'changepoint_mtx') <- change_mtx

  clust_res <- Rphenograph::Rphenograph(change_mtx, k = 30)
  sbird_sce$subclone <- clust_res[[2]]$membership
  return(sbird_sce)
}

#' Title
#'
#' @param sbird_sce
#' @param min_reads
#'
#' @return
#' @export
#'
#' @examples
ploidy_correction <- function(sbird_sce, min_reads = 100000){
  # For each subclone average the ploidy_estimate
  subclones <- unique(sbird_sce$subclone)
  sbird_sce$wgd <- FALSE
  sbird_sce$corr.ploidy <- NA
  for(subclone in subclones){

    selector <- sbird_sce$subclone == subclone &
      sbird_sce$total_reads>min_reads &
      sbird_sce$est_ploidy > 0 &
      sbird_sce$est_ploidy < 10 &
      !is.na(sbird_sce$est_ploidy)
    est_ploidies <- sbird_sce$est_ploidy[selector]
    avg_ploidy <- mean(est_ploidies, na.rm = T)
    if(length(est_ploidies) > 10){

      # Fit a kmeans k=2 to the ploidy estimates to identify wgd cells
      wgd_clustering <- stats::kmeans(est_ploidies, centers = 2)
      wgd_index <- which.max(wgd_clustering$centers)
      if(wgd_clustering$betweenss/wgd_clustering$totss > 0.7){
        sbird_sce$wgd[selector][wgd_clustering$cluster == wgd_index] <- TRUE
      }
    }
    sbird_sce$corr.ploidy[sbird_sce$subclone == subclone] <- avg_ploidy
  }
  sbird_sce$corr.ploidy[sbird_sce$wgd] <- sbird_sce$corr.ploidy[sbird_sce$wgd]*2
  return(sbird_sce)
}

#' Title
#'
#' @param sbird_sce the songbird single cell experiment object
#'
#' @return
#' @export
#'
#' @examples
copyCall <- function(sbird_sce){
  # Fit means to produce the final copy matrix
  cn_matrix <- c()
  segmented_matrix <- SummarizedExperiment::assay(sbird_sce, 'segmented')
  for(i in 1:ncol(segmented_matrix)){
    cn_matrix <- cbind(cn_matrix, fitMeans(segmented_matrix[,i], SummarizedExperiment::rowData(sbird_sce)$overlap_use, sbird_sce$corr.ploidy[i]))
  }
  SummarizedExperiment::assay(sbird_sce, 'copy', withDimnames = F) <- cn_matrix
  return(sbird_sce)
}

#' Title
#'
#' @param res
#'
#' @return
#' @export
#'
#' @examples
create_sce <- function(res){
  assay_cols <- c('uncorrected.reads', 'reads', 'ubh_tx')
  counts <- reads_to_matrix(res, assay_cols[1])
  reads <- reads_to_matrix(res, assay_cols[2])
  segmented <- reads_to_matrix(res, assay_cols[3])

  sbird_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = t(counts), reads = t(reads), segmented = t(segmented)))

  # Per cell file information
  sbird_sce$bam_file <- sapply(res, function(x) x$bam_file[1])
  sbird_sce$bedpe_file <- sapply(res, function(x) x$bedpe_file[1])
  sbird_sce$cell_id <- sapply(res, function(x) x$cell_id[1])

  # Global Cell Quantified Information
  sbird_sce$total_reads <- sapply(res, function(x) sum(x$uncorrected.reads))

  # Cell information derived from just the high quality regions
  sbird_sce$est_ploidy <- sapply(res, function(x) mean(x$est_ploidy))
  sbird_sce$observed_coverage <- sapply(res, function(x) mean(x$coverage, na.rm = T)/mean(x$overlap_genome_size))
  sbird_sce$overlap_genome_size <- sapply(res, function(x) mean(x$overlap_genome_size))
  sbird_sce$doublet_tag_rate <- sapply(res, function(x) mean(x$prop_doublet_tags, na.rm = T))
  sbird_sce$est_genome_size <- sbird_sce$doublet_tag_rate/sbird_sce$observed_coverage

  # Bin information
  SummarizedExperiment::rowData(sbird_sce)$chr <- res[[1]]$chromosome
  SummarizedExperiment::rowData(sbird_sce)$start <- res[[1]]$start
  SummarizedExperiment::rowData(sbird_sce)$end <- res[[1]]$end
  SummarizedExperiment::rowData(sbird_sce)$bin_name <- paste0(res[[1]]$chromosome, ':', res[[1]]$start, '-', res[[1]]$end)
  SummarizedExperiment::rowData(sbird_sce)$overlap_use <- res[[1]]$use
  SummarizedExperiment::rowData(sbird_sce)$gc <- res[[1]]$gc
  SummarizedExperiment::rowData(sbird_sce)$mappability <- res[[1]]$mappability
  return(sbird_sce)
}

# Apply a gaussian kernel to the change matrix
#' Title
#'
#' @param matrix the inputted matrix to smooth over
#' @param n_neighbors number of adjacent bins to smooth with
#'
#' @return
#' @export
#'
#' @examples
gauss_kernel <- function(matrix, n_neighbors){
  # Check that n_neighbors is odd
  if(n_neighbors %% 2 == 0){
    stop('n_neighbors must be odd')
  }

  # Define a linear gaussian kernel
  kernel <- c()
  for(i in 1:n_neighbors){
    kernel <- c(kernel, stats::dnorm(i, mean = (n_neighbors + 1)/2, sd = n_neighbors/6))
  }
  kernel <- kernel/max(kernel)
  span <- n_neighbors%/%2

  # Apply the kernel to the matrix
  kernel_mtx <- matrix(0, nrow = nrow(matrix), ncol = ncol(matrix)-span+1)
  max_idx <- 0
  for(i in 1:nrow(matrix)){
    for(j in (span+1):(ncol(matrix)-span)){
      kernel_mtx[i,j] <- sum(matrix[i, (j-span):(j+span)] * kernel)
    }
  }
  return(kernel_mtx)
}

#' Title
#'
#' @param matrix the inputted matrix to identify changepoints from
#' @param use_mask bins to use
#'
#' @return
#' @export
#'
#' @examples
generate_changepoint_matrix <- function(matrix, use_mask){
  # Find the change points in the mixed matrix
  matrix <- t(matrix)[,use_mask]
  change_mtx <- t(apply(matrix, 1, diff))

  # Smooth the change matrix by applying a gaussian kernel
  kernel_mtx <- round(gauss_kernel(change_mtx, 5), 2)

  # Select change points observed in more than 10% of cells
  bins_to_use <- which(colSums(abs(kernel_mtx)>0.1)>(nrow(kernel_mtx)*0.1))
  kernel_mtx <- kernel_mtx[,bins_to_use]
  return(kernel_mtx)
}
