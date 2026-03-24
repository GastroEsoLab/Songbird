#' hmm_ploidy
#'
#' Estimate per-bin copy number for a single cell using a Poisson-emission HMM.
#' Wavelet breakpoints from ubh_segment are used to initialise the Viterbi path,
#' giving the HMM a warm start that dramatically reduces the risk of getting
#' stuck in poor local optima at very low coverage.
#'
#' The emission model is:
#'   Count.Over_b ~ Poisson( ratio_k * Expected.Over_b )
#' where ratio_k = (k-1)/k for CN state k, and
#'   Expected.Over_b = Count.Upstream_b * min_length / Total.Window.Size_b
#'
#' where Count.Upstream_b is the total number of background reads observed in
#' the upstream window across all reads in bin b, Total.Window.Size_b is the mean
#' background window width in bp across those reads, and min_length is the
#' overlap window width in bp. This correctly scales the background density to
#' the overlap window without the asymmetry introduced by dividing at fragment
#' level when background windows are much larger than the overlap window.
#'
#' The transition matrix has a high self-transition probability (p_stay) that is
#' refined by Baum-Welch. This encodes the prior that CN changes are rare, so
#' the HMM requires accumulated evidence across many bins before switching state.
#'
#' @param bin_data   data.frame — one row per bin, output of estimate.ploidy.
#'                   Must contain: Count.Over, Count.Upstream, Total.Window.Size,
#'                   reads, ubh_tx, chromosome, start, end, use, mappability, bases.
#' @param min_length integer — effective overlap window size in bp, equal to
#'   \code{min_length - tag_overlap}. Read automatically from
#'   \code{metadata(sbird_sce)$overlap_window_size} by \code{ploidy_correction_hmm}.
#' @param cn_states  integer vector — CN states to model. Default 1:6.
#' @param p_stay     numeric — initial self-transition probability. Default 0.999.
#' @param max_iter   integer — maximum Baum-Welch iterations. Default 50.
#' @param tol        numeric — convergence tolerance on transition matrix. Default 1e-4.
#'
#' @return A list with:
#'   \item{cn_sequence}{integer vector of decoded CN state per bin (Viterbi)}
#'   \item{mean_ploidy}{numeric — length-weighted mean CN across good bins}
#'   \item{gamma}{matrix (n_bins x K) of posterior state probabilities}
#'   \item{transition}{K x K learned transition matrix}
#'   \item{bin_data}{input bin_data with cn_state and posterior columns appended}
#' @export
hmm_ploidy <- function(bin_data,
                       min_length = 50,
                       cn_states  = 1:6,
                       p_stay     = 0.999,
                       max_iter   = 50,
                       tol        = 1e-4) {

  K <- length(cn_states)

  # ---- quality mask --------------------------------------------------------
  good <- !is.na(bin_data$Count.Over) &
          !is.na(bin_data$Count.Upstream) &
          !is.na(bin_data$Total.Window.Size) &
          bin_data$Count.Upstream > 0 &
          bin_data$Total.Window.Size > 0 &
          is.finite(bin_data$Count.Upstream) &
          is.finite(bin_data$Total.Window.Size) &
          !is.na(bin_data$mappability) &
          bin_data$mappability >= 90 &
          !is.na(bin_data$bases) &
          bin_data$bases >= 90 &
          bin_data$use

  n_good <- sum(good)
  if (n_good < 10) {
    warning("Fewer than 10 usable bins — returning NA ploidy")
    bin_data$cn_state <- NA_integer_
    return(list(cn_sequence  = rep(NA_integer_, nrow(bin_data)),
                mean_ploidy  = NA_real_,
                gamma        = NULL,
                transition   = NULL,
                bin_data     = bin_data))
  }

  bd <- bin_data[good, ]

  # ---- sort to genomic order -----------------------------------------------
  has_prefix    <- any(grepl('^chr', bd$chr))
  chr_nums      <- if (has_prefix) paste0('chr', c(1:22, 'X', 'Y'))
                   else             c(as.character(1:22), 'X', 'Y')
  chr_factor    <- factor(bd$chr, levels = chr_nums[chr_nums %in% unique(bd$chr)])
  genomic_order <- order(chr_factor, bd$start)
  already_sorted <- identical(genomic_order, seq_len(nrow(bd)))
  if (!already_sorted) bd <- bd[genomic_order, ]
  restore_order <- if (!already_sorted) order(genomic_order) else seq_len(nrow(bd))

  # ---- transition matrix (shared across chromosomes) -----------------------
  A <- .build_trans(p_stay, K)

  # ---- run HMM independently per chromosome --------------------------------
  chrs      <- unique(bd$chr)
  cn_good   <- integer(nrow(bd))
  gamma     <- matrix(1/K, nrow(bd), K)

  for (this_chr in chrs) {
    idx <- which(bd$chr == this_chr)
    if (length(idx) < 3) {
      cn_good[idx] <- cn_states[round(K / 2)]
      next
    }
    bdc      <- bd[idx, ]
    n        <- nrow(bdc)

    # Doublet emission
    ratios   <- (cn_states - 1) / cn_states
    base_exp <- (bdc$Count.Upstream / bdc$Total.Window.Size) * min_length * bdc$Bin.Reads
    mu       <- outer(base_exp, ratios)
    mu[mu <= 0 | !is.finite(mu)] <- 1e-10

    log_emit <- matrix(NA_real_, n, K)
    for (k in seq_len(K))
      log_emit[, k] <- dpois(bdc$Count.Over, lambda = mu[, k], log = TRUE)

    # Wavelet warm start (chromosome-level)
    init_cn <- .wavelet_to_cn(bdc$reads, bdc$ubh_tx, cn_states, chr = bdc$chr)
    pi_init <- tabulate(match(init_cn, cn_states), nbins = K) / n
    pi_init <- pmax(pi_init, 1e-6); pi_init <- pi_init / sum(pi_init)

    res <- .hmm_bw_viterbi(log_emit, A, pi_init, K, n, max_iter, tol)
    cn_good[idx] <- cn_states[res$path]
    gamma[idx, ] <- res$gamma
  }

  # ---- restore original order ----------------------------------------------
  cn_good <- cn_good[restore_order]
  gamma   <- gamma[restore_order, , drop = FALSE]

  # ---- broadcast back to all bins ------------------------------------------
  cn_full            <- rep(NA_integer_, nrow(bin_data))
  cn_full[good]      <- cn_good
  gamma_full         <- matrix(NA_real_, nrow(bin_data), K)
  gamma_full[good, ] <- gamma
  colnames(gamma_full) <- paste0("P.CN", cn_states)

  bin_data$cn_state <- cn_full
  bin_data          <- cbind(bin_data, gamma_full)

  list(
    cn_sequence  = cn_full,
    mean_ploidy  = mean(cn_good, na.rm = TRUE),
    gamma        = gamma_full,
    transition   = A,
    bin_data     = bin_data
  )
}


# ---- helpers ----------------------------------------------------------------

#' Run Baum-Welch EM and Viterbi decoding on a pre-computed log emission matrix
#'
#' Shared by hmm_ploidy and hmm_ploidy_joint. Operates on a single chromosome.
#'
#' @param log_emit  n x K log emission matrix
#' @param A         K x K transition matrix (initial, updated by BW in place)
#' @param pi_init   K-vector initial state distribution (for BW only)
#' @param K         number of states
#' @param n         number of bins
#' @param max_iter  max Baum-Welch iterations (0 = Viterbi only)
#' @param tol       BW convergence tolerance
#' @return list(path = integer n-vector of state indices, gamma = n x K matrix)
#' @noRd
.hmm_bw_viterbi <- function(log_emit, A, pi_init, K, n, max_iter, tol) {

  log_pi <- log(pi_init)
  log_A  <- log(A)
  gamma  <- matrix(1/K, n, K)

  # ---- Baum-Welch ----------------------------------------------------------
  for (iter in seq_len(max_iter)) {
    log_A <- log(A)

    # Forward
    log_alpha      <- matrix(-Inf, n, K)
    log_alpha[1, ] <- log_pi + log_emit[1, ]
    for (t in 2:n)
      for (k in seq_len(K))
        log_alpha[t, k] <- log_emit[t, k] +
          .log_sum_exp(log_alpha[t - 1, ] + log_A[, k])

    # Backward
    log_beta <- matrix(0, n, K)
    for (t in (n - 1):1)
      for (k in seq_len(K))
        log_beta[t, k] <- .log_sum_exp(
          log_A[k, ] + log_emit[t + 1, ] + log_beta[t + 1, ])

    # Gamma
    log_gamma <- log_alpha + log_beta
    log_gamma <- log_gamma - apply(log_gamma, 1, .log_sum_exp)
    gamma     <- exp(log_gamma)

    # Xi -> A update
    log_lik <- .log_sum_exp(log_alpha[n, ])
    A_num   <- matrix(0, K, K)
    for (t in seq_len(n - 1))
      for (j in seq_len(K))
        for (k in seq_len(K))
          A_num[j, k] <- A_num[j, k] + exp(
            log_alpha[t, j] + log_A[j, k] +
            log_emit[t + 1, k] + log_beta[t + 1, k] - log_lik)

    A_new <- A_num / rowSums(A_num)
    A_new <- pmax(A_new, 1e-10)
    A_new <- A_new / rowSums(A_new)

    converged <- max(abs(A_new - A), na.rm = TRUE) < tol
    A <- A_new
    log_A <- log(A)
    if (converged) break
  }

  # ---- Viterbi (uniform start) ---------------------------------------------
  log_uniform    <- rep(log(1/K), K)
  log_delta      <- matrix(-Inf, n, K)
  psi            <- matrix(0L,   n, K)
  log_delta[1, ] <- log_uniform + log_emit[1, ]

  for (t in 2:n)
    for (k in seq_len(K)) {
      scores          <- log_delta[t - 1, ] + log_A[, k]
      best            <- which.max(scores)
      psi[t, k]       <- best
      log_delta[t, k] <- scores[best] + log_emit[t, k]
    }

  path    <- integer(n)
  path[n] <- which.max(log_delta[n, ])
  for (t in (n - 1):1)
    path[t] <- psi[t + 1, path[t + 1]]

  list(path = path, gamma = gamma, A = A)
}

# ---- helpers ----------------------------------------------------------------

#' Map wavelet segmented signal to nearest CN state
#'
#' When \code{chr} is supplied, uses chromosome-level medians rather than
#' bin-level wavelet levels to initialise CN states.  This produces a coarser
#' but more stable warm start that is less affected by wavelet over-segmentation
#' at very low coverage.  Without \code{chr}, falls back to bin-level mapping.
#'
#' @param reads   numeric vector of normalised reads per bin
#' @param ubh_tx  numeric vector of wavelet-segmented read levels per bin
#' @param cn_states integer vector of allowed CN states
#' @param chr     character vector of chromosome names per bin, or NULL for
#'                bin-level mapping (original behaviour)
#' @noRd
.wavelet_to_cn <- function(reads, ubh_tx, cn_states, chr = NULL) {

  n        <- length(ubh_tx)
  finite   <- is.finite(ubh_tx) & ubh_tx > 0

  # Genome-wide scale factor: modal wavelet level -> modal CN state
  # (same in both modes — provides the absolute scale)
  if (sum(finite) == 0) return(rep(cn_states[round(length(cn_states) / 2)], n))

  modal_cn    <- cn_states[round(length(cn_states) / 2)]
  modal_level <- stats::median(ubh_tx[finite], na.rm = TRUE)
  if (!is.finite(modal_level) || modal_level <= 0) modal_level <- 1
  scale <- modal_cn / modal_level

  if (!is.null(chr)) {
    # ---- chromosome-level initialisation ------------------------------------
    # Compute per-chromosome median wavelet level, scale to CN, then assign
    # every bin on that chromosome to the single rounded CN value.
    # Much coarser than bin-level but immune to wavelet over-segmentation.
    chr_levels <- tapply(ubh_tx, chr, function(x) {
      v <- x[is.finite(x) & x > 0]
      if (length(v) == 0) NA_real_ else stats::median(v, na.rm = TRUE)
    })

    # Map each bin to its chromosome median level
    bin_levels <- chr_levels[chr]
    raw_cn     <- as.numeric(bin_levels) * scale

  } else {
    # ---- bin-level initialisation (original behaviour) ----------------------
    raw_cn <- ubh_tx * scale
  }

  # Round and clamp to allowed CN states
  raw_cn  <- ifelse(is.finite(raw_cn), raw_cn, modal_cn)
  clamped <- pmax(min(cn_states), pmin(max(cn_states), round(raw_cn)))
  cn      <- cn_states[vapply(as.integer(clamped),
                              function(x) which.min(abs(cn_states - x)),
                              integer(1))]
  as.integer(cn)
}

#' Build a K x K transition matrix with p_stay on the diagonal
#' @noRd
.build_trans <- function(p_stay, K) {
  A        <- matrix((1 - p_stay) / (K - 1), K, K)
  diag(A)  <- p_stay
  A
}

#' Numerically stable log-sum-exp
#' @noRd
.log_sum_exp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(-Inf)
  m + log(sum(exp(x - m)))
}
