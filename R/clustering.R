#' Title
#'
#' @param sbird_sce
#'
#' @return
#' @export
#'
#' @examples
identify_subclones <- function(sbird_sce, assay = 'segmented', k = 40, method = inv_manhattan, min_size = 5, seed = 1234){
  set.seed(seed)
  # Cluster the cells using the changepoint matrix and sce
  change_mtx <- generate_changepoint_matrix(SummarizedExperiment::assay(sbird_sce, assay), use_mask = SummarizedExperiment::rowData(sbird_sce)$overlap_use)
  SingleCellExperiment::reducedDim(sbird_sce, paste0(assay, '_changepoint')) <- change_mtx

  # Make graph
  knn_graph <- RANN::nn2(change_mtx, k = k)
  knn_matrix <- knn_graph$nn.idx

  # perform clustering
  membership <- clustering(knn_matrix = knn_matrix, feat_matrix = change_mtx, method = method)

  # Clusters under min_size are marked as 0
  subclones <- unique(membership)
  subclone_counts <- sapply(subclones, function(x) sum(membership == x))
  membership[membership %in% subclones[subclone_counts < min_size]] <- 0

  sbird_sce$subclone <- membership
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

# inverse manhattan
inv_manhattan <- function(x, y){
  # x and y are vectors of the same length
  l1_dist <- sum(abs(x-y))
  return(1/(1+l1_dist))
}

convert_mtxlong <- function(mtx){
  mtx <- data.frame(mtx)
  mtx$from <- mtx[,1]
  mtx_long <- tidyr::pivot_longer(mtx, cols = -from, names_to = 'col', values_to = 'value')
  return(mtx_long)
}

clustering <- function(knn_matrix, feat_matrix, method){
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
  graph <- igraph::graph_from_data_frame(relations, directed = F)

  # Perform leiden clustering
  target <- 150
  avg_cellsClust <- 0
  res <- 1
  while(avg_cellsClust < target & res > 0.1){
    # Perform clustering
    clustering <- igraph::cluster_louvain(graph, resolution = res)
    avg_cellsClust <- nrow(knn_matrix)/length(unique(clustering$membership))
    res <- res - 0.1
  }
  return(clustering$membership)
}
