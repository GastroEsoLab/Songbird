#' ploidy_correction_hmm_joint
#'
#' Drop-in replacement for \code{\link{ploidy_correction_hmm}} that uses a
#' joint Poisson-emission HMM combining the tagmentation doublet signal with
#' raw read depth (\code{\link{hmm_ploidy_joint}}).
#'
#' The two emissions are independent Poisson processes:
#' \itemize{
#'   \item \strong{Doublet}: \eqn{y^{\text{over}} \sim \text{Pois}(\frac{k-1}{k} \lambda^{\text{bg}})}
#'   \item \strong{Read depth}: \eqn{y^{\text{counts}} \sim \text{Pois}(k \cdot \mu^{\text{read}})}
#' }
#' The joint log emission is a weighted sum controlled by \code{alpha}.
#' Setting \code{alpha = 0} reproduces \code{\link{ploidy_correction_hmm}}.
#'
#' @param sbird_sce  SingleCellExperiment — output of \code{\link{process.batch}}.
#'   Must have assays \code{"reads"}, \code{"segmented"}, \code{"counts"}
#'   (raw integer reads), \code{"Count.Over"}, \code{"Count.Upstream"},
#'   \code{"Total.Window.Size"}, and \code{"Bin.Reads"}.
#' @param min_reads  integer — minimum reads for WGD detection. Default 10000.
#' @param k          integer — neighbours for pre-clustering. Default 45.
#' @param cn_states  integer vector — CN states modelled. Default 1:6.
#' @param p_stay     numeric — initial self-transition probability. Default 0.999.
#' @param min_length integer — effective overlap window size in bp. Read from
#'   \code{metadata(sbird_sce)$overlap_window_size} if not supplied.
#' @param alpha      numeric in [0, 1] — weight of read depth emission relative
#'   to doublet emission. Default 0.5 (equal weight). Set to 0 for doublet-only
#'   (reproduces \code{\link{ploidy_correction_hmm}}).
#' @param max_iter   integer — Baum-Welch iterations. Default 50. Set to 0
#'   to skip Baum-Welch and run Viterbi on the fixed prior directly.
#' @param n_cpu      integer — cores for parallel processing.
#'
#' @return Same structure as \code{\link{ploidy_correction_hmm}} — input SCE
#'   with \code{"hmm_cn"} assay, \code{hmm_ploidy}, \code{corr.ploidy},
#'   and \code{wgd} columns added.
#' @export
#'
#' @examples
#' \dontrun{
#' sbird_sce <- process.batch(bams, bedpes, genome = 'hg38',
#'                            background_method = 'symmetric',
#'                            ploidy_method = 'hmm')
#' # Joint model with equal weight
#' sbird_sce <- ploidy_correction_hmm_joint(sbird_sce, alpha = 0.5)
#' # Read-depth dominated
#' sbird_sce <- ploidy_correction_hmm_joint(sbird_sce, alpha = 0.8)
#' sbird_sce <- copyCall(sbird_sce)
#' sbird_sce <- identify_subclones(sbird_sce)
#' }
ploidy_correction_hmm_joint <- function(sbird_sce,
                                        min_reads  = 10000,
                                        k          = 45,
                                        cn_states  = 1:6,
                                        p_stay     = 0.999,
                                        min_length = 50,
                                        alpha      = 0.5,
                                        max_iter   = 50,
                                        n_cpu      = NULL) {

  if (is.null(n_cpu)) n_cpu <- max(1L, parallel::detectCores() - 1L)

  # ---- resolve min_length from metadata ------------------------------------
  stored_ows <- S4Vectors::metadata(sbird_sce)$overlap_window_size
  stored_ml  <- S4Vectors::metadata(sbird_sce)$min_length
  if (min_length == 50L) {
    if (!is.null(stored_ows)) {
      min_length <- stored_ows
    } else if (!is.null(stored_ml)) {
      min_length <- stored_ml
    }
  }

  # ---- check required assays -----------------------------------------------
  assay_names  <- SummarizedExperiment::assayNames(sbird_sce)
  overlap_cols <- c("Count.Over", "Count.Upstream", "Total.Window.Size", "Bin.Reads")
  required     <- c("reads", "segmented", "counts", overlap_cols)
  missing      <- setdiff(required, assay_names)
  if (length(missing) > 0) {
    stop(
      "ploidy_correction_hmm_joint requires assays: ",
      paste(missing, collapse = ", "), "\n",
      "Ensure process.batch was run with bedpe files and that the SCE ",
      "contains the 'counts' (raw reads) assay."
    )
  }

  rd           <- as.data.frame(SummarizedExperiment::rowData(sbird_sce))
  missing_rd   <- setdiff(c("mappability", "bases", "use"), colnames(rd))
  if (length(missing_rd) > 0) {
    stop("Missing rowData columns: ", paste(missing_rd, collapse = ", "))
  }

  # ---- extract assay matrices ----------------------------------------------
  reads_mat     <- SummarizedExperiment::assay(sbird_sce, "reads")
  segmented_mat <- SummarizedExperiment::assay(sbird_sce, "segmented")
  counts_mat    <- SummarizedExperiment::assay(sbird_sce, "counts")
  count_over_mat  <- SummarizedExperiment::assay(sbird_sce, "Count.Over")
  count_up_mat    <- SummarizedExperiment::assay(sbird_sce, "Count.Upstream")
  win_size_mat    <- SummarizedExperiment::assay(sbird_sce, "Total.Window.Size")
  bin_reads_mat   <- SummarizedExperiment::assay(sbird_sce, "Bin.Reads")

  n_bins     <- nrow(reads_mat)
  n_cells    <- ncol(reads_mat)
  cell_names <- colnames(reads_mat)

  # ---- bin template: bin-level constants from rowData ----------------------
  bin_template <- rd[, c("mappability", "bases", "use",
                          setdiff(colnames(rd),
                                  c("mappability", "bases", "use",
                                    overlap_cols)))]

  # ---- run joint HMM per cell ----------------------------------------------
  message("Running joint HMM ploidy estimation on ", n_cells, " cells",
          "  (alpha = ", alpha, ")...")

  hmm_results <- pbmcapply::pbmclapply(seq_len(n_cells), function(i) {

    bd <- bin_template
    bd$reads               <- reads_mat[,      i]
    bd$ubh_tx              <- segmented_mat[,   i]
    bd$raw_counts          <- counts_mat[,      i]
    bd$Count.Over          <- count_over_mat[,  i]
    bd$Count.Upstream      <- count_up_mat[,    i]
    bd$Total.Window.Size   <- win_size_mat[,    i]
    bd$Bin.Reads           <- bin_reads_mat[,   i]

    tryCatch(
      hmm_ploidy_joint(bd,
                       min_length = min_length,
                       cn_states  = cn_states,
                       p_stay     = p_stay,
                       max_iter   = max_iter,
                       alpha      = alpha),
      error = function(e) {
        warning("Joint HMM failed for cell ", cell_names[i],
                ": ", conditionMessage(e))
        list(cn_sequence = rep(NA_integer_, n_bins),
             mean_ploidy = NA_real_,
             gamma       = matrix(NA_real_, n_bins, length(cn_states)),
             transition  = NULL,
             bin_data    = NULL)
      }
    )
  }, mc.cores = n_cpu)

  # ---- store results -------------------------------------------------------
  cn_matrix <- do.call(cbind, lapply(hmm_results, `[[`, "cn_sequence"))
  colnames(cn_matrix) <- cell_names
  SummarizedExperiment::assay(sbird_sce, "hmm_cn", withDimnames = FALSE) <- cn_matrix

  sbird_sce$hmm_ploidy <- vapply(hmm_results, `[[`, numeric(1), "mean_ploidy")

  posteriors <- lapply(hmm_results, `[[`, "gamma")
  names(posteriors) <- cell_names
  S4Vectors::metadata(sbird_sce)$hmm_posteriors <- posteriors
  S4Vectors::metadata(sbird_sce)$hmm_alpha      <- alpha

  # ---- cluster-level WGD detection (same as ploidy_correction_hmm) ---------
  sbird_sce <- identify_subclones(sbird_sce,
                                  assay       = "segmented",
                                  k           = k,
                                  column_name = "pc_groups")

  clonal_membership <- SummarizedExperiment::colData(sbird_sce)[["pc_groups"]]
  subclones         <- unique(clonal_membership)

  sbird_sce$wgd         <- FALSE
  sbird_sce$corr.ploidy <- NA_real_

  for (subclone in subclones) {
    clone        <- clonal_membership == subclone
    hmm_ploidies <- sbird_sce$hmm_ploidy[clone]
    read_counts  <- sbird_sce$total_reads[clone]

    high_q <- hmm_ploidies[
      hmm_ploidies > 0 & hmm_ploidies < 10 &
      is.finite(hmm_ploidies) &
      read_counts > min_reads
    ]

    min_cluster <- max(3L, min(10L, floor(sum(clone) / 3L)))

    if (length(high_q) >= min_cluster) {
      sbird_sce$corr.ploidy[clone] <- mean(high_q, na.rm = TRUE)
      sbird_sce$wgd[clone]         <- detect_wgd(high_q, hmm_ploidies)
    } else {
      fallback <- hmm_ploidies[is.finite(hmm_ploidies) & hmm_ploidies > 0]
      if (length(fallback) > 0) {
        sbird_sce$corr.ploidy[clone] <- mean(fallback, na.rm = TRUE)
        message("Subclone ", subclone, ": using per-cell mean ploidy as fallback (",
                length(fallback), " cells).")
      } else {
        warning("Not enough high-quality cells for subclone: ", subclone,
                "\nDefaulting to range 2-8 in copyCall.")
        sbird_sce$corr.ploidy[clone] <- NA_real_
      }
    }
  }

  wgd_idx <- which(sbird_sce$wgd & !is.na(sbird_sce$wgd))
  sbird_sce$corr.ploidy[wgd_idx] <- sbird_sce$corr.ploidy[wgd_idx] * 2

  sbird_sce
}
