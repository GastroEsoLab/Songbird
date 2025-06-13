#' Reads_to_Matrix
#' This function converts a long form data into a matrix
#'
#' @param data list of cell processing results
#' @param column column to select from
#' @param return_counts boolean to return total reads
#'
#' @return a matrix of reads or counts
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

#' Sanitize_Ploidy
#'
#' @param expected_ploidy noisy estimator of true cell ploidy
#'
#' @return a vector of possible modal ploidy for the cell
#'
sanitize_ploidy <- function(expected_ploidy){
  # If using the overlap calling pipeline
  if(!is.na(expected_ploidy)){
    if(expected_ploidy < 0.5){
      expected_ploidy <- 2
    }else if(expected_ploidy > 15){
      expected_ploidy <- 15
    }
    # main detected peak could be any one of 1,2,3,4,5 state - modeled by a poisson mean 2
    modeState <- seq(round(expected_ploidy*0.75), round(expected_ploidy*1.5))
    modeState <- modeState[modeState>0]

  }else{ # If skipping the overlap calling pipeline, assume ploidies from 2 to 8
    modeState <- seq(2, 8)
  }
  return(modeState)
}

#' est_cn
#'
#' @param values numeric vector of segmented read counts to fit
#' @param uniploid estimated ratio of reads per copy number
#' @param sigma estimate of the noise within the bins
#' @param return_cn boolean to return the final CN + LogProb, or to return the BIC
#'
#' @return either a list containing the CN and logProb of fit, or the BIC for tuning
#'
#' @examples
est_cn <- function(values, uniploid, sigma, return_cn = F){

  # Get the total number of states (+1 to account for 0 copy regions) in the dataset
  # given uniploid state & score the gaussians for each
  numStates <- max(round(values/uniploid))+1
  if(numStates>1000){
    numStates <- 1000
  }
  scores <- matrix(nrow = numStates, ncol = length(values))
  for(j in 1:numStates){
    scores[j,] <- stats::dnorm(values,
                               mean = uniploid*(j-1), # Account for 0 copy regions
                               sd = sigma,
                               log = T)
  }

  # Best fitting state is the normal that best fits each state
  states <- apply(scores, 2, which.max)
  topScores <- sapply(1:length(states), function(j) scores[states[j],j])
  if(return_cn){
    out <- list(fit = sum(topScores), states = states - 1)
    return(out)
  }else{
    # If not outputting the just return the BIC
    topScores <- sum(topScores)
    n_cna <- sum(diff(states)!=0)
    BIC <- -2 * sum(topScores) + log(length(states)) * n_cna
    return(BIC)
  }
}

#' fitMeans
#'
#' @param means numeric vector of segmented read counts to fit
#' @param use logical vector indicating which bins are good to use
#' @param sigma estimate of the noise within the bins - calculated by the residual of the reads - segmented state
#' @param expected_ploidy noisy estimator of true cell ploidy - NA results in testing all Cell Ploidies 2-8
#' @param tune_uniploid boolean to tune the reads per CN value to minimize the breakpoints & maximize the fit
#'
#' @return The final CN state for each bin
#'
#' @examples
fitMeans <- function(means, use, sigma, expected_ploidy = NA, tune_uniploid = FALSE){

  modeState <- sanitize_ploidy(expected_ploidy)

  final_cn <- rep(NA, length(means))
  all_means <- means
  use_idx <- which(use & !is.na(means))
  means <- means[use_idx]
  nBins <- length(means)

  # Identify largest peak and propose CN state
  dens <- stats::density(means)
  peak <- abs(dens$x[which.max(dens$y)])

  fitScores <- c() # How well does the mu fit the fixed states applied to them (dnormal around each state)
  stateScores <- c() # How well does the inferred states fit our strong prior that cells should be diploid?
  fitStates <- matrix(nrow = length(modeState), ncol = nBins)

  for(i in 1:length(modeState)){
    # Given mode state - find uniploid state
    uniploid <- peak/modeState[i]
    if(tune_uniploid){
      lb <- (uniploid * modeState[i])/(modeState[i] + 0.5)
      ub <- (uniploid * modeState[i])/(modeState[i] - 0.5)
      best_uniploid <- stats::optim(par = uniploid, fn = est_cn, values = means, sigma = sigma,
                             return_cn = F, method = 'Brent', lower = lb, upper = ub)
      uniploid <- best_uniploid$par
    }
    res <- est_cn(means, uniploid, sigma, return_cn = T)
    fitScores <- c(fitScores, res$fit)
    fitStates[i,] <- res$states
    stateScores <- c(stateScores, sum(stats::dnorm(mean = res$states, sd = 1, expected_ploidy, log = T)))

  }

  # Get the highest scoring fit and apply it to the good bins
  stateScores <- sapply(stateScores, function(x) x-matrixStats::logSumExp(stateScores))
  fitScores <- sapply(fitScores, function(x) x-matrixStats::logSumExp(fitScores))

  # If there is no prior expected ploidy, place no weight on the state scores
  bestFit <- which.max(sum(stateScores, fitScores, na.rm = T))
  print(bestFit)
  final_cn[use_idx] <- fitStates[bestFit,]

  # get true uniploid & apply to the remaining bins
  best_uniploid <- peak/modeState[bestFit]
  final_cn[!use_idx] <- round(all_means[!use_idx]/best_uniploid)
  return(final_cn)
}


#' ploidy_correction
#'
#' @param sbird_sce the songbird object after subclone calling
#' @param min_reads minimum number of reads to use a cell's ploidy estimation for correction
#' @param column the column in colData to use for subclone membership, default is 'subclone'
#'
#' @return songbird object with corrected ploidy estimates
#' @export
#'
#' @examples
ploidy_correction <- function(sbird_sce, min_reads = 100000, column = 'subclone'){
  # For each subclone average the ploidy_estimate
  clonal_membership <- SingleCellExperiment::colData(sbird_sce)[[column]]
  subclones <- unique(clonal_membership)
  sbird_sce$wgd <- FALSE
  sbird_sce$corr.ploidy <- NA
  for(subclone in subclones){
    clone <- clonal_membership == subclone
    est_ploidies <- sbird_sce$est_ploidy[clone]
    readCounts <- sbird_sce$total_reads[clone]

    high_qPloidies <- est_ploidies[est_ploidies > 0 & est_ploidies < 10 & is.finite(est_ploidies) & readCounts > min_reads]
    if(length(high_qPloidies) > 10){
      sbird_sce$corr.ploidy[clone] <- mean(high_qPloidies, na.rm = T)
      #sbird_sce$wgd[clone] <- detect_wgd(high_qPloidies, est_ploidies)
    }else{
      warning(paste0("Not enough high quality cells to estimate ploidy for subclone: ", subclone, '\nDefaulting to a range from 2-8'))
      sbird_sce$corr.ploidy[clone] <- NA
    }
  }
  sbird_sce$corr.ploidy[sbird_sce$wgd] <- sbird_sce$corr.ploidy[sbird_sce$wgd]*2
  return(sbird_sce)
}

#' detect_wgd
#'
#' @param high_qPloidies vector of ploidies that are used for correction
#' @param all_ploidies all ploidies
#'
#' @return the WGD identity of cells - only detects 1 whole genome duplication
#'
#' @examples
detect_wgd <- function(high_qPloidies, all_ploidies){
  wgd_clustering <- stats::kmeans(high_qPloidies, centers = 2)
  wgd_index <- which.max(wgd_clustering$centers)

  # Find which center each ploidy is closest to and assign
  distances <- sapply(wgd_clustering$centers, function(x) abs(x - all_ploidies))
  closest <- apply(distances, 1, which.min)
  closest[is.na(closest)] <- which.min(wgd_clustering$centers)
  if(wgd_clustering$betweenss/wgd_clustering$totss > 0.7){
    wgd_cells <- closest == wgd_index
  }else{
    wgd_cells <- rep(FALSE, length(all_ploidies))
  }
  return(wgd_cells)
}

#' copyCall
#'
#' @param sbird_sce the songbird single cell experiment object - requires that 'segmented' and 'reads' assays are present
#' @param num_cores number of cores to use for parallel processing, default is all but one core
#' @param tune_uniploid whether to tune the uniploid value for each bin, default is FALSE
#'
#' @return songbird object with the copy number called for each bin
#' @export
#'
#' @examples
copyCall <- function(sbird_sce, num_cores = NULL, tune_uniploid = FALSE){
  if(is.null(num_cores)){
    num_cores <- parallel::detectCores() - 1
  }
  # Fit means to produce the final copy matrix
  cn_matrix <- c()
  segmented_matrix <- SummarizedExperiment::assay(sbird_sce, 'segmented')
  reads_matrix <- SummarizedExperiment::assay(sbird_sce, 'reads')

  use <- SummarizedExperiment::rowData(sbird_sce)$overlap_use
  var_matrix <- segmented_matrix[use,] - reads_matrix[use,]
  sigmas <- apply(var_matrix, 2, function(x) stats::sd(x, na.rm = T))
  cn_matrix <- parallel::mclapply(1:ncol(segmented_matrix), function(i){
    fitMeans(segmented_matrix[,i], use, sigmas[i], sbird_sce$corr.ploidy[i], tune_uniploid)
  }, mc.cores = num_cores)

  cn_matrix <- do.call(cbind, cn_matrix)

  SummarizedExperiment::assay(sbird_sce, 'copy', withDimnames = F) <- cn_matrix
  sbird_sce$bin_noise <- sigmas
  return(sbird_sce)
}

#' create_sce
#'
#' @param res a list of results from cell processing
#'
#' @return the initial songbird object
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
  sbird_sce$state_variance <- sapply(res, function(x) sqrt(sum((x$reads - x$ubh_tx)^2, na.rm = T)/nrow(x)))

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

#' create_sce_from_res
#'
#' @param res a list of results similar to the cell processing results
#' @param n_cpu number of cores to use for parallel processing, default is all but one core
#'
#' @return songbird object with the reads, counts, and segmented data
#'
#' @examples
create_sce_from_res <- function(res, n_cpu=NULL){
  if(is.null(n_cpu)){
    n_cpu <- parallel::detectCores() - 1
  }

  reads <- reads_to_matrix(res, 'reads')
  counts <- reads_to_matrix(res, 'uncorrected.reads')
  segmented <- pbmcapply::pbmclapply(res, function(x) Songbird::ubh_segment(x$reads, use = TRUE, min_svSize = 4), mc.cores = n_cpu)
  segmented <- do.call(rbind, segmented)

  sbird_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = t(counts), reads = t(reads), segmented = t(segmented)))

  # Per cell file information
  sbird_sce$bam_file <- NA
  sbird_sce$bedpe_file <- NA
  sbird_sce$cell_id <- sapply(res, function(x) unique(x$cell_id))
  sbird_sce$total_reads <- apply(counts, 1, function(x) sum(x, na.rm = T))

  # Cell information derived from just the high quality regions
  sbird_sce$est_ploidy <- NULL
  sbird_sce$observed_coverage <- NULL
  sbird_sce$overlap_genome_size <- NULL
  sbird_sce$doublet_tag_rate <- NULL
  sbird_sce$est_genome_size <- NULL

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

