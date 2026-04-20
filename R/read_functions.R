#' @importFrom stats pnorm qnorm median var
NULL

# -----------------------------------------------------------------------
#' Validate and bundle paired human/AI parameter estimates
#'
#' Takes two mle lists (one per scoring condition) and returns a validated
#' \code{aidif_data} object for use in \code{\link{fit_aidif}}.
#'
#' @param human_mle An mle list for human-scored data.  Must contain
#'   \code{est} (a named list \code{group.1}, \code{group.2} of
#'   data.frames with columns \code{a1}, \code{d1}) and \code{var.cov}
#'   (matching list of covariance matrices).
#' @param ai_mle An mle list for AI-scored data in the same format.
#'
#' @return A list of class \code{"aidif_data"} with elements
#'   \code{human} and \code{ai}.
#'
#' @seealso \code{\link{fit_aidif}}, \code{\link{make_aidif_eg}},
#'   \code{\link{simulate_aidif_data}}
#' @export
read_ai_scored <- function(human_mle, ai_mle) {
  check_aidif_mle(human_mle, label = "human_mle")
  check_aidif_mle(ai_mle,    label = "ai_mle")
  check_compatible_mles(human_mle, ai_mle)
  structure(list(human = human_mle, ai = ai_mle), class = "aidif_data")
}
