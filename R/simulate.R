# -----------------------------------------------------------------------
#' Simulate item parameter estimates for the AI-DIF model.
#'
#' Generates a synthetic \code{aidif_data}-compatible list suitable for
#' benchmarking and method evaluation.  The data-generating model contains:
#' classical DIF in the human scoring condition (controlled via
#' \code{dif_items} and \code{dif_mag}), differential AI scoring bias
#' (controlled via \code{dasb_items} and \code{dasb_mag}), and a latent
#' group mean difference (\code{impact}).
#'
#' Rather than simulating item responses and refitting IRT models (which
#' requires additional dependencies), this function directly simulates
#' maximum-likelihood estimates and their asymptotic covariance matrices,
#' consistent with a 2PL model fitted to \code{n_obs} observations per
#' group.
#'
#' @param n_items   Integer. Number of items. Default: \code{10}.
#' @param n_obs     Integer. Approximate number of observations per group,
#'   used to scale the covariance matrices. Default: \code{500}.
#' @param impact    Numeric. Latent mean difference (group 2 minus group 1)
#'   in SD units. Default: \code{0.5}.
#' @param dif_items Integer vector. Indices of items with DIF in the human
#'   scoring condition (intercept shift added to group 2). Default:
#'   \code{1}.
#' @param dif_mag   Numeric. Magnitude of the intercept DIF effect (in
#'   IRT metric). Default: \code{0.5}.
#' @param dasb_items Integer vector. Indices of items where AI scoring
#'   introduces differential bias (intercept shift added to group 2 in the
#'   AI condition only). Default: \code{3}.
#' @param dasb_mag  Numeric. Magnitude of the DASB effect. Default:
#'   \code{0.4}.
#' @param ai_drift  Numeric. Uniform intercept shift applied to ALL items in
#'   BOTH groups under AI scoring (simulates AI calibration offset).
#'   Default: \code{0.1}.
#' @param seed      Integer seed for reproducibility, or \code{NULL}.
#'   Default: \code{42}.
#'
#' @return A list with elements \code{human} and \code{ai}, each formatted
#'   identically to the output of
#'   \code{\link{read_ai_scored}}.  Can be passed directly to
#'   \code{\link{fit_aidif}}.
#'
#' @examples
#' dat <- simulate_aidif_data(
#'   n_items   = 8,
#'   n_obs     = 600,
#'   dif_items = c(1, 2),
#'   dasb_items = 5
#' )
#' mod <- fit_aidif(dat$human, dat$ai)
#' summary(mod)
#'
#' @export
# -----------------------------------------------------------------------

simulate_aidif_data <- function(n_items    = 10L,
                                n_obs      = 500L,
                                impact     = 0.5,
                                dif_items  = 1L,
                                dif_mag    = 0.5,
                                dasb_items = 3L,
                                dasb_mag   = 0.4,
                                ai_drift   = 0.1,
                                seed       = 42L) {

  if (!is.null(seed)) set.seed(seed)

  # ---- True item parameters (group 1 / human) ---------------------------
  a_true <- stats::runif(n_items, 0.6, 1.8)          # slopes
  d_true <- stats::rnorm(n_items, 0,   0.6)           # intercepts (difficulty)

  # ---- Sampling variability (approximate asymptotic SEs) ----------------
  # For a 2PL item with a ~ 1, d ~ 0, n ~ 500:
  # SE(a) ~ 0.06,  SE(d) ~ 0.07 (rough approximation)
  se_a <- 0.06 * sqrt(500 / n_obs)
  se_d <- 0.07 * sqrt(500 / n_obs)

  make_est_df <- function(a_vec, d_vec) {
    df <- data.frame(a1 = a_vec, d1 = d_vec)
    rownames(df) <- paste0("item", seq_len(n_items))
    df
  }

  make_vcov <- function(n_items, se_a, se_d) {
    n_par <- n_items * 2L
    vcov  <- diag(c(rbind(se_a^2, se_d^2)), n_par)
    # Add small off-diagonals within item (a-d covariance ~ -0.3 * se_a * se_d)
    for (i in seq_len(n_items)) {
      idx_a <- (i - 1L) * 2L + 1L
      idx_d <- idx_a + 1L
      cov_ad <- -0.3 * se_a * se_d
      vcov[idx_a, idx_d] <- cov_ad
      vcov[idx_d, idx_a] <- cov_ad
    }
    vcov
  }

  # ---- Human scoring parameters -----------------------------------------
  # Group 1 (reference): true parameters + sampling noise
  a_h1 <- a_true + stats::rnorm(n_items, 0, se_a)
  d_h1 <- d_true + stats::rnorm(n_items, 0, se_d)

  # Group 2 (focal): shift intercepts by -impact (higher latent trait =>
  # items appear easier, i.e. d shifts negatively), then add DIF
  d_h2_base <- d_true - impact + stats::rnorm(n_items, 0, se_d)
  a_h2      <- a_true + stats::rnorm(n_items, 0, se_a)

  # Add human-scoring DIF to focal group
  d_h2 <- d_h2_base
  d_h2[dif_items] <- d_h2[dif_items] + dif_mag

  # ---- AI scoring parameters --------------------------------------------
  # AI scores shift ALL items uniformly by ai_drift (calibration bias)
  a_a1 <- a_h1 + stats::rnorm(n_items, 0, se_a * 0.5)
  d_a1 <- d_h1 + ai_drift + stats::rnorm(n_items, 0, se_d * 0.5)

  a_a2 <- a_h2 + stats::rnorm(n_items, 0, se_a * 0.5)
  d_a2 <- d_h2 + ai_drift + stats::rnorm(n_items, 0, se_d * 0.5)

  # Differential AI Scoring Bias: shift focal group only at dasb_items
  d_a2[dasb_items] <- d_a2[dasb_items] + dasb_mag

  # ---- Assemble mle lists -----------------------------------------------
  vcov1 <- make_vcov(n_items, se_a, se_d)
  vcov2 <- make_vcov(n_items, se_a, se_d)

  human <- list(
    par.names = list(
      internal = make_par_names(n_items),
      original = make_par_names(n_items)
    ),
    est = list(
      group.1 = make_est_df(a_h1, d_h1),
      group.2 = make_est_df(a_h2, d_h2)
    ),
    var.cov = list(group.1 = vcov1, group.2 = vcov2)
  )

  ai <- list(
    par.names = list(
      internal = make_par_names(n_items),
      original = make_par_names(n_items)
    ),
    est = list(
      group.1 = make_est_df(a_a1, d_a1),
      group.2 = make_est_df(a_a2, d_a2)
    ),
    var.cov = list(group.1 = vcov1, group.2 = vcov2)
  )

  list(human = human, ai = ai)
}

# -----------------------------------------------------------------------
#' Built-in example dataset for aiDIF
#'
#' Constructs and returns the built-in example dataset: paired
#' human and AI item parameter estimates for 6 items in two groups,
#' with known DIF and DASB planted at specific items.
#'
#' The data-generating model includes:
#' \itemize{
#'   \item \strong{Item 1}: DIF under human scoring (intercept +0.5
#'     in focal group).
#'   \item \strong{Item 3}: Differential AI Scoring Bias (DASB) —
#'     AI scoring adds +0.4 to the focal-group intercept only.
#'   \item \strong{Impact}: 0.5 SD (focal group higher on latent trait).
#'   \item \strong{AI drift}: uniform +0.1 calibration offset on all
#'     items in both groups.
#' }
#'
#' @return A list with elements \code{human} and \code{ai}, each a
#'   validated mle list (see \code{\link{simulate_aidif_data}} for
#'   format details).
#'
#' @examples
#' eg  <- make_aidif_eg()
#' mod <- fit_aidif(eg$human, eg$ai)
#' summary(mod)
#'
#' @seealso \code{\link{simulate_aidif_data}}, \code{\link{fit_aidif}}
#' @export
# -----------------------------------------------------------------------

make_aidif_eg <- function() {
  se_a <- 0.06; se_d <- 0.07; n <- 6L
  a_true <- c(1.2, 0.9, 1.4, 1.0, 1.3, 0.8)
  d_true <- c(0.3, -0.5, 0.1, -0.2, 0.4, -0.1)

  mk_df <- function(a, d) {
    df <- data.frame(a1 = round(a, 4), d1 = round(d, 4))
    rownames(df) <- paste0("item", seq_len(n)); df
  }
  mk_vc <- function() {
    v <- diag(c(rbind(se_a^2, se_d^2)), n * 2L)
    for (i in seq_len(n)) {
      ia <- (i - 1L)*2L + 1L; id <- ia + 1L
      v[ia, id] <- v[id, ia] <- -0.3 * se_a * se_d
    }; v
  }

  a_h1 <- a_true + c( .05,-.03, .04,-.02, .03,-.04)
  d_h1 <- d_true + c( .04,-.06, .03, .05,-.04, .02)
  a_h2 <- a_true + c(-.03, .04,-.02, .05,-.03, .02)
  d_h2 <- d_true - 0.5 + c(.03,-.05,.04,-.02,.06,-.03)
  d_h2[1L] <- d_h2[1L] + 0.5       # item 1: human DIF

  a_a1 <- a_h1 + c( .02,-.01, .03,-.01, .02,-.02)
  d_a1 <- d_h1 + 0.1 + c(.03,-.02,.01,.02,-.03,.01)
  a_a2 <- a_h2 + c(-.02, .03,-.01, .02,-.02, .01)
  d_a2 <- d_h2 + 0.1 + c(.02,-.03,.04,-.01,.03,-.02)
  d_a2[3L] <- d_a2[3L] + 0.4       # item 3: DASB

  pn  <- make_par_names(n)
  vc  <- mk_vc()
  mle <- function(g1, g2) list(
    par.names = list(internal = pn, original = pn),
    est       = list(group.1 = g1, group.2 = g2),
    var.cov   = list(group.1 = vc, group.2 = vc)
  )
  list(human = mle(mk_df(a_h1,d_h1), mk_df(a_h2,d_h2)),
       ai    = mle(mk_df(a_a1,d_a1), mk_df(a_a2,d_a2)))
}

# -----------------------------------------------------------------------
# Internal helper: build parameter name vectors for mle lists.
# -----------------------------------------------------------------------

make_par_names <- function(n_items) {
  unlist(lapply(seq_len(n_items), function(i) {
    c(paste0("item", i, ".a1"), paste0("item", i, ".d1"))
  }))
}
