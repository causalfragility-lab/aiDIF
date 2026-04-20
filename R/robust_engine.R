#' @importFrom stats pnorm qnorm median var
#' @importFrom Matrix diag t bdiag
NULL

# ============================================================
# SECTION 1 — IRT SCALING FUNCTIONS
# ============================================================
# These compute item-level DIF statistics from paired group
# parameter estimates.  The estimand is a location parameter
# theta that represents the "typical" across-group difference
# after robust down-weighting of biased items.
# ============================================================

#' Compute item-level IRT scaling functions
#'
#' For each item, computes a standardised difference between
#' group parameter estimates.  The result is a vector \code{y}
#' whose robust location is estimated by
#' \code{\link{estimate_robust_scale}}.
#'
#' @param mle A validated \code{mle} list (output of
#'   \code{\link{read_ai_scored}} or constructed manually).
#'   Must contain \code{est$group.1}, \code{est$group.2} and
#'   matching \code{var.cov} matrices.
#' @param type One of \code{"intercept"} (default) or
#'   \code{"slope"}.  Determines which parameters are compared.
#' @param scale_by One of \code{"pooled"} (default),
#'   \code{"ref"}, or \code{"focal"}.  Controls the
#'   denominator used to standardise intercept differences:
#'   \code{"pooled"} uses \eqn{\sqrt{(a_1^2+a_2^2)/2}};
#'   \code{"ref"} uses \eqn{a_1}; \code{"focal"} uses
#'   \eqn{a_2}.  Ignored when \code{type = "slope"}.
#'
#' @return A named numeric vector of scaling-function values,
#'   one entry per item threshold (or per item for slopes).
#'
#' @keywords internal
compute_scaling_fn <- function(mle,
                               type     = "intercept",
                               scale_by = "pooled") {
  
  g1 <- mle$est$group.1
  g2 <- mle$est$group.2
  n_items      <- nrow(g1)
  n_thresholds <- ncol(g1) - 1L   # first column is slope
  
  if (type == "slope") {
    y <- g2$a1 / g1$a1
    names(y) <- paste0(rownames(g1), "_a")
    return(y)
  }
  
  # --- intercept differences -------------------------------------------
  num <- as.matrix(g2[, -1, drop = FALSE]) -
    as.matrix(g1[, -1, drop = FALSE])
  
  denom <- switch(scale_by,
                  pooled = sqrt((g1$a1^2 + g2$a1^2) / 2),
                  ref    = g1$a1,
                  focal  = g2$a1,
                  stop("`scale_by` must be 'pooled', 'ref', or 'focal'.", call. = FALSE)
  )
  
  y <- c(t(num / denom))
  names(y) <- paste0(
    rep(rownames(g1), each = n_thresholds),
    "_d",
    rep(seq_len(n_thresholds), times = n_items)
  )
  y
}


#' Gradient of the intercept scaling function
#'
#' Computes the Jacobian of \code{\link{compute_scaling_fn}}
#' with respect to all item parameters in both groups,
#' organised to be conformable with the block-diagonal
#' covariance matrix built by \code{\link{build_joint_vcov}}.
#'
#' @param mle A validated mle list.
#' @param theta Optional scalar: if supplied, item-specific
#'   scaling-function values are replaced by \code{theta} in
#'   the gradient (used when evaluating under H0).
#' @param scale_by Passed from \code{\link{compute_scaling_fn}}.
#'
#' @return A matrix with \code{n_items * n_thresholds} columns,
#'   each being the gradient vector of one scaling-function
#'   entry with respect to the full parameter vector.
#'
#' @keywords internal
grad_intercept_fn <- function(mle, theta = NULL, scale_by = "pooled") {
  
  g1 <- mle$est$group.1
  g2 <- mle$est$group.2
  n_items      <- nrow(g1)
  n_thresholds <- ncol(g1) - 1L
  n_pars       <- n_items * (n_thresholds + 1L) * 2L
  
  a1 <- rep(g1$a1, each = n_thresholds)
  a2 <- rep(g2$a1, each = n_thresholds)
  
  if (is.null(theta)) {
    theta_vec <- compute_scaling_fn(mle, "intercept", scale_by)
  } else {
    theta_vec <- rep(theta[1L], n_items * n_thresholds)
  }
  
  # Four gradient components per scaling-function entry,
  # ordered as (a1_i, d1_i, a2_i, d2_i):
  grad_block <- matrix(0, nrow = 4L, ncol = n_items * n_thresholds)
  if (scale_by == "pooled") {
    denom <- sqrt((a1^2 + a2^2) / 2)
    grad_block[1L, ] <- -a1 / (a1^2 + a2^2) * theta_vec
    grad_block[2L, ] <- -1 / denom
    grad_block[3L, ] <- -a2 / (a1^2 + a2^2) * theta_vec
    grad_block[4L, ] <-  1 / denom
  } else if (scale_by == "ref") {
    grad_block[1L, ] <- -theta_vec / a1
    grad_block[2L, ] <- -1 / a1
    grad_block[3L, ] <-  0
    grad_block[4L, ] <-  1 / a1
  } else {  # focal
    grad_block[1L, ] <-  0
    grad_block[2L, ] <- -1 / a2
    grad_block[3L, ] <- -theta_vec / a2
    grad_block[4L, ] <-  1 / a2
  }
  
  # Expand each column into a full n_pars-length vector
  half <- n_pars / 2L
  k <- 1L
  grad_list <- vector("list", n_items * n_thresholds)
  template  <- numeric(n_pars)
  for (i in seq_len(n_items)) {
    for (j in seq_len(n_thresholds)) {
      pos_a1 <- (i - 1L) * (n_thresholds + 1L) + 1L
      pos_d1 <- pos_a1 + j
      pos_a2 <- pos_a1 + half
      pos_d2 <- pos_d1 + half
      template[c(pos_a1, pos_d1, pos_a2, pos_d2)] <- grad_block[, k]
      grad_list[[k]] <- template
      template[] <- 0
      k <- k + 1L
    }
  }
  do.call(cbind, grad_list)
}


#' Block-diagonal joint covariance matrix for both groups
#'
#' @param mle A validated mle list.
#' @return A \code{Matrix::bdiag} sparse block-diagonal matrix.
#' @keywords internal
build_joint_vcov <- function(mle) {
  Matrix::bdiag(mle$var.cov$group.1, mle$var.cov$group.2)
}


#' Delta-method covariance matrix of the scaling functions
#'
#' @param mle  A validated mle list.
#' @param theta Optional scalar evaluated under H0.
#' @param scale_by Passed to gradient function.
#'
#' @return A square matrix of size
#'   \code{n_items * n_thresholds}.
#'
#' @keywords internal
vcov_scaling_fn <- function(mle, theta = NULL, scale_by = "pooled") {
  G   <- grad_intercept_fn(mle, theta, scale_by)
  Sig <- build_joint_vcov(mle)
  Matrix::t(G) %*% Sig %*% G
}


# ============================================================
# SECTION 2 — BI-SQUARE M-ESTIMATION (IRLS)
# ============================================================

#' Bi-square weight function
#' @param u Numeric vector of standardised residuals.
#' @param k Tuning parameter (default 1.96).
#' @return Numeric vector of weights in [0, 1].
#' @keywords internal
bisq_weight <- function(u, k = 1.96) {
  w <- (u / k)^2
  out <- (1 - w)^2
  out[abs(u) > k] <- 0
  out
}

#' Bi-square rho (objective) function
#' @keywords internal
bisq_rho <- function(u, k = 1.96) {
  w   <- (u / k)^2
  out <- 1 - (1 - w)^3
  out[abs(u) > k] <- 1
  out
}

#' Bi-square psi (influence) function
#' @keywords internal
bisq_psi <- function(u, k = 1.96) {
  w   <- (u / k)^2
  out <- u * (1 - w)^2
  out[abs(u) > k] <- 0
  out
}

#' Derivative of bi-square psi
#' @keywords internal
bisq_psi_prime <- function(u, k = 1.96) {
  w   <- (u / k)^2
  out <- (1 - w)^2 - 4 * w * (1 - w)
  out[abs(u) > k] <- 0
  out
}

#' Least trimmed squares estimate of location (starting value)
#' @param y Numeric vector.
#' @param trim Proportion to trim (default 0.5).
#' @keywords internal
lts_location <- function(y, trim = 0.5) {
  n <- length(y)
  h <- max(1L, floor(n * (1 - trim)))
  if (h >= n) return(mean(y))
  ys <- sort(y)
  wins <- n - h + 1L
  vars <- vapply(seq_len(wins),
                 function(j) var(ys[j:(j + h - 1L)]),
                 numeric(1))
  j_min <- which.min(vars)
  mean(ys[j_min:(j_min + h - 1L)])
}

#' Grid search for bi-square objective minimum (starting value)
#' @keywords internal
grid_rho_search <- function(y, var_fn, k, width = 0.01) {
  theta_grid <- seq(max(min(y), -2), min(max(y), 2), by = width)
  rho_vals <- vapply(theta_grid, function(th) {
    var_y <- var_fn(th)
    u     <- (y - th) / sqrt(var_y)
    sum(bisq_rho(u, k))
  }, numeric(1))
  theta_grid[which.min(rho_vals)]
}


#' Robust DIF scale estimation via IRLS
#'
#' Estimates a robust location parameter for the vector of
#' IRT scaling functions using iteratively re-weighted least
#' squares (IRLS) with the bi-square loss.  This is the core
#' estimation engine of \pkg{aiDIF}.
#'
#' @param mle   A validated mle list.
#' @param alpha Significance level controlling the bi-square
#'   tuning parameter \eqn{k = z_{1-\alpha/2}}. Default 0.05.
#' @param scale_by Scaling denominator; passed to
#'   \code{\link{compute_scaling_fn}}.  Default \code{"pooled"}.
#' @param tol   Convergence tolerance. Default \code{1e-7}.
#' @param maxit Maximum IRLS iterations. Default \code{100}.
#'
#' @return A list of class \code{rdif_fit} with elements:
#' \describe{
#'   \item{\code{est}}{Estimated robust scale parameter.}
#'   \item{\code{weights}}{Bi-square item weights.}
#'   \item{\code{rho_value}}{Value of objective at solution.}
#'   \item{\code{n_iter}}{Number of iterations used.}
#'   \item{\code{k}}{Tuning parameter used.}
#'   \item{\code{y}}{Raw scaling function values.}
#'   \item{\code{vcov_est}}{Covariance matrix of \code{y} at solution.}
#'   \item{\code{dif_test}}{Wald item-level DIF test (data.frame).}
#'   \item{\code{dtf_test}}{Wald test of differential test functioning.}
#' }
#'
#' @examples
#' dat <- simulate_aidif_data(n_items = 5, seed = 1)
#' fit <- estimate_robust_scale(dat$human)
#' print(fit$est)
#'
#' @export
estimate_robust_scale <- function(mle,
                                  alpha    = 0.05,
                                  scale_by = "pooled",
                                  tol      = 1e-7,
                                  maxit    = 100L) {
  check_aidif_mle(mle)
  
  k <- qnorm(1 - alpha / 2)
  y <- compute_scaling_fn(mle, "intercept", scale_by)
  
  # Variance function evaluated at a given theta
  var_at <- function(th) Matrix::diag(vcov_scaling_fn(mle, th, scale_by))
  
  # Three starting values: median, LTS, grid minimum
  starts <- c(
    median(y),
    lts_location(y),
    grid_rho_search(y, var_at, k)
  )
  
  run_irls <- function(start) {
    theta <- start
    nit   <- 0L
    conv  <- 1
    while (nit < maxit && conv > tol) {
      var_y   <- var_at(theta)
      u       <- (y - theta) / sqrt(var_y)
      w_star  <- bisq_weight(u, k)
      w       <- (w_star / var_y) / sum(w_star / var_y)
      new_th  <- sum(w * y)
      conv    <- abs(theta - new_th)
      theta   <- new_th
      nit     <- nit + 1L
    }
    var_y <- var_at(theta)
    u     <- (y - theta) / sqrt(var_y)
    list(est       = theta,
         weights   = bisq_weight(u, k),
         rho_value = sum(bisq_rho(u, k)),
         n_iter    = nit)
  }
  
  sols     <- lapply(starts, run_irls)
  rho_vals <- sapply(sols, `[[`, "rho_value")
  best     <- sols[[which.min(rho_vals)]]
  
  # Covariance at solution
  Vcov <- vcov_scaling_fn(mle, best$est, scale_by)
  
  best$y        <- y
  best$k        <- k
  best$vcov_est <- Vcov
  best$dif_test <- wald_dif_test(y, best$est, Vcov)
  best$dtf_test <- wald_dtf_test(y, best$est, best$weights,
                                 Vcov, vcov_scaling_fn(mle, NULL, scale_by), k)
  class(best)   <- "rdif"
  best
}


# ============================================================
# SECTION 3 — WALD TESTS
# ============================================================

#' Wald item-level DIF test
#'
#' Tests H0: y_i = theta for each item, using the projected
#' variance that accounts for estimation of theta itself.
#'
#' @param y      Scaling function values.
#' @param theta  Estimated robust scale parameter.
#' @param Vcov   Covariance matrix of y (at theta under H0).
#'
#' @return A data.frame with columns \code{delta}, \code{se},
#'   \code{z}, \code{p_val}.
#'
#' @keywords internal
wald_dif_test <- function(y, theta, Vcov) {
  n     <- length(y)
  delta <- y - theta
  var_y <- Matrix::diag(Vcov)
  
  # Projection removes the contribution of theta estimation
  P     <- matrix(1 / var_y / sum(1 / var_y), nrow = n, ncol = n)
  I_min_P <- diag(1, n) - P
  se    <- sqrt(pmax(Matrix::diag(Matrix::t(I_min_P) %*% Vcov %*% I_min_P), 0))
  z     <- ifelse(se > 0, delta / se, NA_real_)
  p     <- 2 * (1 - pnorm(abs(z)))
  
  data.frame(delta = delta, se = se, z = z, p_val = p,
             row.names = names(y))
}


#' Wald test of differential test functioning (DTF)
#'
#' Tests H0: mean(y) - theta = 0, i.e. whether the robust
#' scale estimate differs significantly from the naive mean.
#' A significant result indicates meaningful DTF.
#'
#' @param y          Scaling function values.
#' @param theta      Robust scale estimate.
#' @param weights    Bi-square weights from IRLS.
#' @param Vcov_H0    Covariance of y under H0 (theta plugged in).
#' @param Vcov_raw   Covariance of y (no theta substitution).
#' @param k          Bi-square tuning parameter.
#'
#' @return A one-row data.frame.
#' @keywords internal
wald_dtf_test <- function(y, theta, weights, Vcov_H0, Vcov_raw, k) {
  n     <- length(y)
  ybar  <- mean(y)
  delta <- ybar - theta
  
  var_y   <- Matrix::diag(Vcov_H0)
  u       <- (y - theta) / sqrt(var_y)
  pp      <- bisq_psi_prime(u, k)
  
  v_mean <- rep(1 / n, n)
  v_rob  <- (pmax(pp, 0) / var_y) / sum(pmax(pp, 0) / var_y)
  v_diff <- v_mean - v_rob
  
  se_mean  <- sqrt(as.numeric(Matrix::t(v_mean) %*% Vcov_raw %*% v_mean))
  se_rob   <- sqrt(as.numeric(Matrix::t(v_rob)  %*% Vcov_raw %*% v_rob))
  se_delta <- sqrt(pmax(as.numeric(Matrix::t(v_diff) %*% Vcov_raw %*% v_diff), 0))
  
  z   <- if (se_delta > 0) delta / se_delta else NA_real_
  p   <- 2 * (1 - pnorm(abs(z)))
  
  data.frame(naive_est = ybar, naive_se = se_mean,
             robust_est = theta, robust_se = se_rob,
             delta = delta, delta_se = se_delta,
             z = z, p_val = p)
}