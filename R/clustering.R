#' Identify Subclones using breakpoint identification and leiden clustering. Clusters with fewer than min_size cells are marked as noisy and put in cluster 0.
#'
#' @param sbird_sce songbird object
#' @param assay the assay to use for clustering
#' @param k number of nearest neighbors to use for clustering
#' @param method the clustering method to use, default is inverse manhattan distance
#' @param min_size minimum size of a subclone to be considered, default is 5 cells
#' @param seed random seed for reproducibility
#' @param min_readCount minimum read count for a cell to be considered, default is 1e5
#' @param column_name the column name to store the subclone information in, default is 'subclone'
#'
#' @return A SingleCellExperiment object with clonal membership
#' @export
#'
#' @examples
identify_subclones <- function(sbird_sce, assay = 'copy', k = 30, res = 'auto', method = inv_manhattan, min_size = 5, seed = 1234, min_readCount = 1e5, column_name = 'subclone'){
  set.seed(seed)
  # Cluster the cells using the changepoint matrix and sce
  change_mtx <- generate_changepoint_matrix(SummarizedExperiment::assay(sbird_sce, assay),
                                            bin_mask = SummarizedExperiment::rowData(sbird_sce)$overlap_use,
                                            cell_mask = sbird_sce$total_reads > min_readCount)
  SingleCellExperiment::reducedDim(sbird_sce, paste0(assay, '_changepoint')) <- change_mtx

  # Make graph
  knn_graph <- RANN::nn2(change_mtx, k = k)
  knn_matrix <- knn_graph$nn.idx

  # perform clustering
  membership <- clustering(knn_matrix = knn_matrix, feat_matrix = change_mtx, method = method, res = res)

  # Clusters under min_size are marked as 0
  subclones <- unique(membership)
  subclone_counts <- sapply(subclones, function(x) sum(membership == x))
  membership[membership %in% subclones[subclone_counts < min_size]] <- 0

  # Add the subclone information to the SingleCellExperiment object
  if(column_name %in% colnames(SummarizedExperiment::colData(sbird_sce))){
    warning(paste0('Column ', column_name, ' already exists in the SingleCellExperiment object. Overwriting.'))
  }
  SummarizedExperiment::colData(sbird_sce)[[column_name]] <- membership
  return(sbird_sce)
}

#' gauss_kernel
#'
#' @param matrix the inputted matrix to smooth over
#' @param n_neighbors number of adjacent bins to smooth with
#'
#' @return a smoothed matrix
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

#' generate_changepoint_matrix
#'
#' @param matrix the inputted matrix to identify changepoints from
#' @param bin_mask bins to use
#' @param cell_mask cells to use
#'
#' @return a matrix of changepoints
generate_changepoint_matrix <- function(matrix, bin_mask, cell_mask){
  # Find the change points in the mixed matrix
  matrix <- t(matrix)[,bin_mask]
  change_mtx <- t(apply(matrix, 1, diff))

  # Smooth the change matrix by applying a gaussian kernel
  kernel_mtx <- round(gauss_kernel(change_mtx, 5), 2)

  # Select change points observed in more than 10% of cells
  highq_mtx <- kernel_mtx[cell_mask,,drop=F]
  bins_to_use <- which(colSums(abs(highq_mtx)>0.1)>(nrow(highq_mtx)*0.1))
  kernel_mtx <- kernel_mtx[,bins_to_use,drop=F]
  return(kernel_mtx)
}

#' inv_manhattan
#'
#' @param x a vector of changepoints
#' @param y a vector of changepoints
#'
#' @return a matrix of changepoints
#' @export
inv_manhattan <- function(x, y){
  # x and y are vectors of the same length
  l1_dist <- sum(abs(x-y))
  return(1/(1+l1_dist))
}

#' convert_mtxlong
#'
#' @param mtx a matrix to convert to long format
#' @return a long format data frame
convert_mtxlong <- function(mtx){
  mtx <- data.frame(mtx)
  mtx$from <- mtx[,1]
  mtx_long <- tidyr::pivot_longer(mtx, cols = -from, names_to = 'col', values_to = 'value')
  return(mtx_long)
}

#' leiden
#'
#' @param graph
#' @param res
#' @return_modularity (default = TRUE)
#'
#' @return either modularity or the membership
leiden <- function(graph, res, return_modularity = TRUE){
  # Run the Leiden algorithm
  clustering <- igraph::cluster_leiden(graph, resolution = res,
                                       objective_function = 'modularity',
                                       n_iterations = 5)
  membership <- clustering$membership

  if(return_modularity){
    return(-igraph::modularity(graph, membership = membership))
  } else {
    return(membership)
  }
}


#' clustering
#'
#' @param knn_matrix a matrix of k-nearest neighbors
#' @param feat_matrix a matrix of features for each cell
#' @param method a function to calculate the affinity between cells
#'
#' @return a vector of cluster membership for each cell
clustering <- function(knn_matrix, feat_matrix, method, res = 'auto'){
  affy_matrix <- matrix(0, nrow = nrow(knn_matrix), ncol = ncol(knn_matrix))
  for(i in 1:nrow(knn_matrix)){
    for(j in 1:ncol(knn_matrix)){
      neigh_idx <- knn_matrix[i, j]
      affy_matrix[i,j] <- method(feat_matrix[i,], feat_matrix[neigh_idx,])
    }
  }

  affy_long <- convert_mtxlong(affy_matrix)
  knn_long <- convert_mtxlong(knn_matrix)
  relations <- data.frame(from = knn_long$from, to = knn_long$value, weight = affy_long$value)
  relations <- relations[relations$from!=relations$to, ]  # Remove self-loops
  graph <- igraph::graph_from_data_frame(relations, directed = F)

  if(res == 'auto'){
    # Automatically determine the resolution parameter based on modularity
    optim_res <- stats::optim(par = 1, fn = leiden, graph = graph, return_modularity = TRUE, method = 'L-BFGS-B', lower = 0.0001, upper = 10)
    optim_res <- optim_res$par
    message(paste0('Optimized resolution parameter: ', round(optim_res, 3)))
  }
  else {
    optim_res <- res
  }
  membership <- leiden(graph, optim_res, return_modularity = FALSE)
  return(membership)
}
