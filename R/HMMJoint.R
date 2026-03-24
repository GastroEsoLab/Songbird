#' hmm_ploidy_joint
#'
#' Estimate per-bin copy number for a single cell using a joint Poisson-emission
#' HMM that combines two independent signals:
#'
#' \enumerate{
#'   \item \strong{Doublet ratio} — same as \code{\link{hmm_ploidy}}.
#'         \eqn{y_b^{\text{over}} \sim \text{Poisson}(\frac{k-1}{k} \cdot \lambda_b^{\text{bg}})}
#'   \item \strong{Raw read depth} — uncorrected integer read counts per bin,
#'         which are linear in CN after accounting for GC/mappability via the
#'         estimated \emph{uniploid} parameter.
#'         \eqn{y_b^{\text{counts}} \sim \text{Poisson}(k \cdot \mu_b^{\text{read}})}
#' }
#'
#' The joint log emission is the sum of the two independent log likelihoods.
#' The two signals are complementary: the doublet ratio is ploidy-scale-free
#' and immune to GC/mappability bias but saturates at high CN; the read depth
#' is linear in CN and discriminates high states well but requires normalisation.
#' Together they substantially improve segmentation at low coverage.
#'
#' The \emph{uniploid} parameter (expected counts at CN 1) is estimated from
#' the corrected \code{reads} assay as
#' \code{median(reads) / modal_cn}, where \code{modal_cn} is the middle of
#' \code{cn_states}.  This borrows the GC/mappability normalisation already
#' applied by QDNAseq without requiring re-correction.
#'
#' @param bin_data   data.frame — one row per bin.  Must contain all columns
#'   required by \code{\link{hmm_ploidy}} plus \code{raw_counts} (integer
#'   uncorrected read counts per bin).
#' @param min_length integer — effective overlap window size in bp.
#' @param cn_states  integer vector — CN states to model. Default 1:6.
#' @param p_stay     numeric — initial self-transition probability. Default 0.999.
#' @param max_iter   integer — maximum Baum-Welch iterations. Default 50.
#'   Set to 0 to skip Baum-Welch and run Viterbi on the fixed prior directly.
#' @param tol        numeric — convergence tolerance. Default 1e-4.
#' @param alpha      numeric in [0,1] — weight of the read depth emission
#'   relative to the doublet emission.  Default 0.5 (equal weight).
#'   Set to 0 to reproduce the doublet-only HMM; set to 1 for read-depth only.
#'
#' @return Same structure as \code{\link{hmm_ploidy}}.
#' @export
hmm_ploidy_joint <- function(bin_data,
                             min_length = 50,
                             cn_states  = 1:6,
                             p_stay     = 0.999,
                             max_iter   = 50,
                             tol        = 1e-4,
                             alpha      = 0.5) {

  K <- length(cn_states)

  # ---- quality mask — same as hmm_ploidy -----------------------------------
  good <- !is.na(bin_data$Count.Over) &
          !is.na(bin_data$Count.Upstream) &
          !is.na(bin_data$Total.Window.Size) &
          bin_data$Count.Upstream > 0 &
          bin_data$Total.Window.Size > 0 &
          is.finite(bin_data$Count.Upstream) &
          is.finite(bin_data$Total.Window.Size) &
          !is.na(bin_data$mappability) &
          # bin_data$mappability >= 90 &
          !is.na(bin_data$bases) &
          # bin_data$bases >= 90 &
          bin_data$use &
          !is.na(bin_data$raw_counts)

  n_good <- sum(good)
  if (n_good < 10) {
    warning("Fewer than 10 usable bins — returning NA ploidy")
    bin_data$cn_state <- NA_integer_
    return(list(cn_sequence = rep(NA_integer_, nrow(bin_data)),
                mean_ploidy = NA_real_,
                gamma       = NULL,
                transition  = NULL,
                bin_data    = bin_data))
  }

  bd <- bin_data[good, ]

  # ---- sort to genomic order -----------------------------------------------
  # Handles both prefixed (chr1) and unprefixed (1) chromosome names.
  has_prefix    <- any(grepl('^chr', bd$chr))
  chr_nums      <- if (has_prefix) paste0('chr', c(1:22, 'X', 'Y'))
                   else             c(as.character(1:22), 'X', 'Y')
  chr_factor    <- factor(bd$chr,
                          levels = chr_nums[chr_nums %in% unique(bd$chr)])
  genomic_order <- order(chr_factor, bd$start)
  already_sorted <- identical(genomic_order, seq_len(nrow(bd)))
  if (!already_sorted) {
    bd <- bd[genomic_order, ]
  }
  restore_order <- if (!already_sorted) order(genomic_order) else seq_len(nrow(bd))

  n  <- nrow(bd)

  # ---- uniploid estimation (genome-wide, used inside per-chr loop) ---------
  modal_cn     <- cn_states[round(length(cn_states) / 2)]
  uniploid     <- stats::median(bd$reads, na.rm = TRUE) / modal_cn
  if (!is.finite(uniploid) || uniploid <= 0) uniploid <- 1
  scale_factor <- mean(bd$raw_counts, na.rm = TRUE) /
                  pmax(mean(bd$reads, na.rm = TRUE), 1e-6)
  uniploid_raw <- uniploid * scale_factor

  # ---- transition matrix (shared prior across chromosomes) -----------------
  A <- .build_trans(p_stay, K)

  # ---- run HMM independently per chromosome --------------------------------
  chrs    <- unique(bd$chr)
  cn_good <- integer(nrow(bd))
  gamma   <- matrix(1/K, nrow(bd), K)

  for (this_chr in chrs) {
    idx <- which(bd$chr == this_chr)
    if (length(idx) < 3) {
      cn_good[idx] <- cn_states[round(K / 2)]
      next
    }
    bdc <- bd[idx, ]
    nc  <- nrow(bdc)

    # Joint log emission for this chromosome
    ratios_c   <- (cn_states - 1) / cn_states
    base_exp_c <- (bdc$Count.Upstream / bdc$Total.Window.Size) * min_length * bdc$Bin.Reads
    mu_over_c  <- outer(base_exp_c, ratios_c)
    mu_over_c[mu_over_c <= 0 | !is.finite(mu_over_c)] <- 1e-10

    mu_reads_c <- outer(rep(uniploid_raw, nc), cn_states)
    mu_reads_c[mu_reads_c <= 0 | !is.finite(mu_reads_c)] <- 1e-10

    log_emit_c <- matrix(NA_real_, nc, K)
    for (k in seq_len(K)) {
      log_emit_c[, k] <-
        (1 - alpha) * dpois(bdc$Count.Over, lambda = mu_over_c[, k], log = TRUE) +
        alpha       * dpois(as.integer(round(bdc$raw_counts)),
                            lambda = mu_reads_c[, k], log = TRUE)
    }

    # Chromosome warm start
    init_cn_c <- .wavelet_to_cn(bdc$reads, bdc$ubh_tx, cn_states, chr = bdc$chr)
    pi_init_c <- tabulate(match(init_cn_c, cn_states), nbins = K) / nc
    pi_init_c <- pmax(pi_init_c, 1e-6); pi_init_c <- pi_init_c / sum(pi_init_c)

    res <- .hmm_bw_viterbi(log_emit_c, A, pi_init_c, K, nc, max_iter, tol)
    cn_good[idx] <- cn_states[res$path]
    gamma[idx, ] <- res$gamma
  }

  # ---- restore original bin order ------------------------------------------
  cn_good <- cn_good[restore_order]
  gamma   <- gamma[restore_order, , drop = FALSE]

  # ---- broadcast back ------------------------------------------------------
  cn_full <- rep(NA_integer_, nrow(bin_data))
  cn_full[good] <- cn_good

  gamma_full <- matrix(NA_real_, nrow(bin_data), K)
  gamma_full[good, ] <- gamma
  colnames(gamma_full) <- paste0("P.CN", cn_states)

  bin_data$cn_state <- cn_full
  bin_data          <- cbind(bin_data, gamma_full)

  list(
    cn_sequence = cn_full,
    mean_ploidy = mean(cn_good, na.rm = TRUE),
    gamma       = gamma_full,
    transition  = A,
    bin_data    = bin_data
  )
}
