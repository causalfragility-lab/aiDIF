# -----------------------------------------------------------------------
#' Fit the AI-DIF model
#'
#' The primary estimation function of \pkg{aiDIF}.  Runs the robust DIF
#' procedure under both human and AI scoring using the built-in IRLS engine
#' (\code{\link{estimate_robust_scale}}), then tests for Differential AI
#' Scoring Bias (DASB).
#'
#' @param human_mle A validated mle list for human-scored data.
#' @param ai_mle    A validated mle list for AI-scored data, or \code{NULL}.
#' @param alpha     Significance level.  Default \code{0.05}.
#' @param scale_by  Denominator for standardising intercept differences:
#'   \code{"pooled"} (default), \code{"ref"}, or \code{"focal"}.
#' @param tol       IRLS convergence tolerance.  Default \code{1e-7}.
#' @param maxit     Maximum IRLS iterations.  Default \code{100}.
#'
#' @return An object of class \code{"aidif"}.
#'
#' @examples
#' dat <- simulate_aidif_data(n_items = 6, seed = 1)
#' mod <- fit_aidif(dat$human, dat$ai)
#' print(mod)
#' summary(mod)
#'
#' @seealso \code{\link{estimate_robust_scale}},
#'   \code{\link{scoring_bias_test}}, \code{\link{simulate_aidif_data}}
#' @export
# -----------------------------------------------------------------------

fit_aidif <- function(human_mle,
                      ai_mle   = NULL,
                      alpha    = 0.05,
                      scale_by = "pooled",
                      tol      = 1e-7,
                      maxit    = 100L) {

  check_aidif_mle(human_mle, label = "human_mle")
  if (!is.null(ai_mle)) {
    check_aidif_mle(ai_mle, label = "ai_mle")
    check_compatible_mles(human_mle, ai_mle)
  }

  human_fit <- estimate_robust_scale(human_mle, alpha, scale_by, tol, maxit)
  dif_human <- human_fit$dif_test

  ai_fit <- NULL; dif_ai <- NULL
  scoring_bias <- NULL; ai_effect <- NULL

  if (!is.null(ai_mle)) {
    ai_fit       <- estimate_robust_scale(ai_mle, alpha, scale_by, tol, maxit)
    dif_ai       <- ai_fit$dif_test
    scoring_bias <- scoring_bias_test(human_mle, ai_mle)
    ai_effect    <- ai_effect_summary(dif_human, dif_ai, alpha)
  }

  structure(
    list(human_fit    = human_fit,
         ai_fit       = ai_fit,
         dif_human    = dif_human,
         dif_ai       = dif_ai,
         scoring_bias = scoring_bias,
         ai_effect    = ai_effect,
         data         = list(human = human_mle, ai = ai_mle),
         alpha        = alpha,
         scale_by     = scale_by),
    class = "aidif"
  )
}
