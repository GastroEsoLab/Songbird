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
#' @param min_length integer — overlap window size in bp (tag_overlap from process.batch).
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
  # Use the same bin selection as estimate.ploidy: mappability + bases + use flag
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

  # Work on the subset of good bins; results will be broadcast back at the end
  bd <- bin_data[good, ]
  n  <- nrow(bd)

  # ---- emission: expected overlap count per bin ----------------------------
  # ratio_k = (k-1)/k  (fraction of cuts that are doublets at CN=k)
  # Expected overlap count = background_density x overlap_window_size
  #   = (Count.Upstream / Total.Window.Size) x min_length
  # Using the raw Count.Upstream and mean window size avoids the asymmetry
  # introduced by dividing at fragment level when background windows (e.g.
  # 50,000 bp symmetric) are much larger than the overlap window (min_length bp).
  ratios   <- (cn_states - 1) / cn_states
  base_exp <- (bd$Count.Upstream / bd$Total.Window.Size) * min_length * bd$Bin.Reads  # n vector

  # mu[b, k] = expected Count.Over for bin b under CN state k
  mu        <- outer(base_exp, ratios)          # n x K
  mu[mu <= 0 | !is.finite(mu)] <- 1e-10

  # log emission matrix — computed once, never changes (background is fixed)
  log_emit <- matrix(NA_real_, n, K)
  for (k in seq_len(K)) {
    log_emit[, k] <- dpois(bd$Count.Over, lambda = mu[, k], log = TRUE)
  }

  # ---- wavelet initialisation of state sequence ----------------------------
  # ubh_tx is the wavelet-segmented read signal. We map its piecewise-constant
  # levels to CN states by rounding (reads / median_reads) * modal_CN.
  # This gives Viterbi a warm start rather than a uniform prior.
  init_cn <- .wavelet_to_cn(bd$reads, bd$ubh_tx, cn_states)

  # ---- transition matrix ---------------------------------------------------
  A <- .build_trans(p_stay, K)

  # Initial state distribution from wavelet assignment
  pi_init        <- tabulate(match(init_cn, cn_states), nbins = K) / n
  pi_init        <- pmax(pi_init, 1e-6)
  pi_init        <- pi_init / sum(pi_init)

  # ---- Baum-Welch ----------------------------------------------------------
  log_pi <- log(pi_init)

  for (iter in seq_len(max_iter)) {

    log_A <- log(A)

    # Forward pass
    log_alpha        <- matrix(-Inf, n, K)
    log_alpha[1, ]   <- log_pi + log_emit[1, ]
    for (t in 2:n) {
      for (k in seq_len(K)) {
        log_alpha[t, k] <- log_emit[t, k] +
          .log_sum_exp(log_alpha[t - 1, ] + log_A[, k])
      }
    }

    # Backward pass
    log_beta      <- matrix(0, n, K)
    for (t in (n - 1):1) {
      for (k in seq_len(K)) {
        log_beta[t, k] <- .log_sum_exp(
          log_A[k, ] + log_emit[t + 1, ] + log_beta[t + 1, ]
        )
      }
    }

    # Posterior state probabilities (gamma)
    log_gamma <- log_alpha + log_beta
    log_gamma <- log_gamma - apply(log_gamma, 1, .log_sum_exp)
    gamma     <- exp(log_gamma)

    # Posterior transition counts (xi), accumulated directly into A_num
    log_lik_total <- .log_sum_exp(log_alpha[n, ])
    A_num         <- matrix(0, K, K)
    for (t in seq_len(n - 1)) {
      for (j in seq_len(K)) {
        for (k in seq_len(K)) {
          A_num[j, k] <- A_num[j, k] + exp(
            log_alpha[t, j] + log_A[j, k] +
            log_emit[t + 1, k] + log_beta[t + 1, k] - log_lik_total
          )
        }
      }
    }

    # M step
    A_new        <- A_num / rowSums(A_num)
    A_new        <- pmax(A_new, 1e-10)
    A_new        <- A_new / rowSums(A_new)

    if (max(abs(A_new - A), na.rm = TRUE) < tol) {
      A <- A_new
      break
    }
    A <- A_new
  }

  # ---- Viterbi decoding (warm-started) ------------------------------------
  log_delta        <- matrix(-Inf, n, K)
  psi              <- matrix(0L,   n, K)
  log_delta[1, ]   <- log(pi_init) + log_emit[1, ]

  for (t in 2:n) {
    for (k in seq_len(K)) {
      scores          <- log_delta[t - 1, ] + log_A[, k]
      best            <- which.max(scores)
      psi[t, k]       <- best
      log_delta[t, k] <- scores[best] + log_emit[t, k]
    }
  }

  path    <- integer(n)
  path[n] <- which.max(log_delta[n, ])
  for (t in (n - 1):1) {
    path[t] <- psi[t + 1, path[t + 1]]
  }

  cn_good <- cn_states[path]

  # ---- broadcast back to all bins (NA for filtered bins) ------------------
  cn_full <- rep(NA_integer_, nrow(bin_data))
  cn_full[good] <- cn_good

  gamma_full <- matrix(NA_real_, nrow(bin_data), K)
  gamma_full[good, ] <- gamma
  colnames(gamma_full) <- paste0("P.CN", cn_states)

  bin_data$cn_state   <- cn_full
  bin_data            <- cbind(bin_data, gamma_full)

  mean_ploidy <- mean(cn_good, na.rm = TRUE)

  list(
    cn_sequence  = cn_full,
    mean_ploidy  = mean_ploidy,
    gamma        = gamma_full,
    transition   = A,
    bin_data     = bin_data
  )
}


# ---- helpers ----------------------------------------------------------------

#' Map wavelet segmented signal to nearest CN state
#' @noRd
.wavelet_to_cn <- function(reads, ubh_tx, cn_states) {
  # Identify the modal (most common) wavelet level — this is the dominant CN
  levels     <- unique(ubh_tx[is.finite(ubh_tx) & ubh_tx > 0])
  if (length(levels) == 0) return(rep(2L, length(reads)))

  level_counts <- vapply(levels, function(l) sum(ubh_tx == l, na.rm = TRUE),
                         integer(1))
  modal_level  <- levels[which.max(level_counts)]

  # Assume the modal level corresponds to the modal CN state (typically 2)
  modal_cn     <- cn_states[round(length(cn_states) / 2)]
  scale        <- modal_cn / modal_level

  # Map each bin's wavelet level to the nearest allowed CN state
  raw_cn  <- round(ubh_tx * scale)
  clamped <- pmax(min(cn_states), pmin(max(cn_states), raw_cn))
  cn      <- cn_states[vapply(clamped, function(x) which.min(abs(cn_states - x)),
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
