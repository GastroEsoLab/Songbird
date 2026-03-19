#' ploidy_correction_hmm
#'
#' Drop-in replacement for \code{\link{ploidy_correction}} that estimates
#' per-bin copy number using a Poisson-emission HMM rather than a global
#' ratio average.
#'
#' Instead of producing a single scalar \code{est_ploidy} per cell, this
#' function runs \code{\link{hmm_ploidy}} on each cell and stores the decoded
#' CN sequence as a new assay (\code{"hmm_cn"}) in the SCE.  The
#' genome-wide mean ploidy from the HMM is written to \code{colData} as
#' \code{hmm_ploidy}, and WGD detection proceeds the same way as before
#' (k-means on the per-cell mean ploidy estimates within each cluster).
#'
#' The function also stores the per-bin posterior state probabilities as a
#' list in \code{metadata(sbird_sce)$hmm_posteriors} — one matrix per cell,
#' rows = bins, columns = CN states.
#'
#' @param sbird_sce  SingleCellExperiment — output of \code{process.batch}.
#'                   Must have assays \code{"reads"} and \code{"segmented"},
#'                   and \code{colData} columns \code{total_reads} and
#'                   \code{est_ploidy}.
#'                   Additionally requires per-bin overlap count assays:
#'                   \code{Count.Over}, \code{Count.Upstream}, \code{Total.Window.Size},
#'                   \code{Bin.Reads}; and \code{rowData} columns \code{mappability},
#'                   \code{bases}, \code{use}.
#'                   These are present when \code{process.batch} was called with
#'                   bedpe files (i.e. the overlap-based ploidy pipeline ran).
#' @param min_reads  integer — minimum total reads for a cell to contribute its
#'                   HMM ploidy to the cluster-level WGD detection. Default 50000.
#' @param k          integer — neighbours for the pre-clustering used by WGD
#'                   detection (same role as in \code{ploidy_correction}). Default 45.
#' @param cn_states  integer vector — CN states modelled by the HMM. Default 1:6.
#' @param p_stay     numeric — initial HMM self-transition probability. Default 0.999.
#' @param min_length integer — overlap window size in bp, should match the
#'                   \code{tag_overlap} used in \code{process.batch}. Default 50.
#' @param n_cpu      integer — cores for parallel processing. Default all minus one.
#'
#' @return The input \code{sbird_sce} with:
#' \itemize{
#'   \item new assay \code{"hmm_cn"} — per-bin integer CN state matrix (bins x cells)
#'   \item \code{colData} column \code{hmm_ploidy} — genome-wide mean CN per cell
#'   \item \code{colData} column \code{corr.ploidy} — cluster-corrected ploidy
#'         (same semantics as original \code{ploidy_correction})
#'   \item \code{colData} column \code{wgd} — logical WGD flag per cell
#'   \item \code{metadata(sbird_sce)$hmm_posteriors} — list of posterior matrices
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' sbird_sce <- process.batch(bams = bams, bedpes = bedpes, genome = 'hg38')
#' sbird_sce <- ploidy_correction_hmm(sbird_sce)
#' sbird_sce <- copyCall(sbird_sce)
#' sbird_sce <- identify_subclones(sbird_sce)
#' }
ploidy_correction_hmm <- function(sbird_sce,
                                  min_reads  = 50000,
                                  k          = 45,
                                  cn_states  = 1:6,
                                  p_stay     = 0.999,
                                  min_length = 50,
                                  n_cpu      = NULL) {

  if (is.null(n_cpu)) n_cpu <- max(1L, parallel::detectCores() - 1L)

  # Read min_length from SCE metadata if not explicitly supplied.
  # min_length must match the value used in process.batch — it defines
  # the overlap window size in the emission model.
  if (min_length == 50L) {
    stored <- S4Vectors::metadata(sbird_sce)$min_length
    if (!is.null(stored)) min_length <- stored
  }

  # ---- check required data is present ------------------------------------
  rd           <- as.data.frame(SummarizedExperiment::rowData(sbird_sce))
  assay_names  <- SummarizedExperiment::assayNames(sbird_sce)
  overlap_cols <- c("Count.Over", "Count.Upstream", "Total.Window.Size", "Bin.Reads")

  # Overlap counts are stored as assay matrices (bins x cells) by create_sce
  # when bedpe files were provided. bases and use are scalar rowData columns.
  missing_assays <- setdiff(overlap_cols, assay_names)
  missing_rd     <- setdiff(c("mappability", "bases", "use"), colnames(rd))

  if (length(missing_assays) > 0 || length(missing_rd) > 0) {
    stop(
      "ploidy_correction_hmm requires per-bin overlap count assays and rowData columns.\n",
      if (length(missing_assays) > 0)
        paste0("Missing assays: ", paste(missing_assays, collapse = ", "), "\n"),
      if (length(missing_rd) > 0)
        paste0("Missing rowData: ", paste(missing_rd, collapse = ", "), "\n"),
      "These are only present when process.batch was run with bedpe files.\n",
      "Fall back to ploidy_correction() if bedpe data is unavailable."
    )
  }

  reads_mat     <- SummarizedExperiment::assay(sbird_sce, "reads")
  segmented_mat <- SummarizedExperiment::assay(sbird_sce, "segmented")
  count_over_mat    <- SummarizedExperiment::assay(sbird_sce, "Count.Over")
  count_up_mat      <- SummarizedExperiment::assay(sbird_sce, "Count.Upstream")
  win_size_mat      <- SummarizedExperiment::assay(sbird_sce, "Total.Window.Size")
  bin_reads_mat     <- SummarizedExperiment::assay(sbird_sce, "Bin.Reads")
  n_bins        <- nrow(reads_mat)
  n_cells       <- ncol(reads_mat)
  cell_names    <- colnames(reads_mat)

  # ---- build per-bin data template (bin-level constants only) -------------
  bin_template <- rd[, c("mappability", "bases", "use",
                          setdiff(colnames(rd),
                                  c("mappability", "bases", "use",
                                    overlap_cols)))]

  # ---- run HMM per cell in parallel ---------------------------------------
  message("Running HMM ploidy estimation on ", n_cells, " cells...")

  hmm_results <- pbmcapply::pbmclapply(seq_len(n_cells), function(i) {

    # Assemble per-cell bin_data: read counts and overlap counts all come
    # from assay matrices; rowData supplies the bin-level constants.
    bd <- bin_template

    bd$reads               <- reads_mat[, i]
    bd$ubh_tx              <- segmented_mat[, i]
    bd$Count.Over      <- count_over_mat[, i]
    bd$Count.Upstream  <- count_up_mat[, i]
    bd$Total.Window.Size <- win_size_mat[, i]
    bd$Bin.Reads       <- bin_reads_mat[, i]

    tryCatch(
      hmm_ploidy(bd,
                 min_length = min_length,
                 cn_states  = cn_states,
                 p_stay     = p_stay),
      error = function(e) {
        warning("HMM failed for cell ", cell_names[i], ": ", conditionMessage(e))
        list(cn_sequence  = rep(NA_integer_, n_bins),
             mean_ploidy  = NA_real_,
             gamma        = matrix(NA_real_, n_bins, length(cn_states)),
             transition   = NULL,
             bin_data     = NULL)
      }
    )
  }, mc.cores = n_cpu)

  # ---- store HMM CN sequences as assay ------------------------------------
  cn_matrix <- do.call(cbind, lapply(hmm_results, `[[`, "cn_sequence"))
  colnames(cn_matrix) <- cell_names

  SummarizedExperiment::assay(sbird_sce, "hmm_cn", withDimnames = FALSE) <- cn_matrix

  # ---- per-cell mean ploidy -----------------------------------------------
  sbird_sce$hmm_ploidy <- vapply(hmm_results, `[[`, numeric(1), "mean_ploidy")

  # ---- store posteriors in metadata ---------------------------------------
  posteriors <- lapply(hmm_results, `[[`, "gamma")
  names(posteriors) <- cell_names
  S4Vectors::metadata(sbird_sce)$hmm_posteriors <- posteriors

  # ---- cluster-level WGD detection (same logic as ploidy_correction) ------
  # Pre-cluster on segmented data so we can share ploidy info across noisy cells
  sbird_sce <- identify_subclones(sbird_sce,
                                  assay       = "segmented",
                                  k           = k,
                                  column_name = "pc_groups")

  clonal_membership <- SummarizedExperiment::colData(sbird_sce)[["pc_groups"]]
  subclones         <- unique(clonal_membership)

  sbird_sce$wgd        <- FALSE
  sbird_sce$corr.ploidy <- NA_real_

  for (subclone in subclones) {
    clone       <- clonal_membership == subclone
    hmm_ploidies <- sbird_sce$hmm_ploidy[clone]
    read_counts  <- sbird_sce$total_reads[clone]

    high_q <- hmm_ploidies[
      hmm_ploidies > 0 & hmm_ploidies < 10 &
      is.finite(hmm_ploidies) &
      read_counts > min_reads
    ]

    if (length(high_q) >= 10) {
      sbird_sce$corr.ploidy[clone] <- mean(high_q, na.rm = TRUE)
      sbird_sce$wgd[clone]         <- detect_wgd(high_q, hmm_ploidies)
    } else {
      warning("Not enough high-quality cells to estimate ploidy for subclone: ",
              subclone, "\nDefaulting to range 2-8 in copyCall.")
      sbird_sce$corr.ploidy[clone] <- NA_real_
    }
  }

  # Adjust for WGD
  wgd_idx <- which(sbird_sce$wgd & !is.na(sbird_sce$wgd))
  sbird_sce$corr.ploidy[wgd_idx] <- sbird_sce$corr.ploidy[wgd_idx] * 2

  sbird_sce
}
