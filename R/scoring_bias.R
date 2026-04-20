# -----------------------------------------------------------------------
#' Differential AI Scoring Bias (DASB) test.
#'
#' For each item, computes the change in item intercept from human to AI
#' scoring within each group, then tests whether this scoring shift differs
#' significantly across groups.  A significant result indicates the AI
#' scoring engine introduces a \emph{group-dependent} parameter distortion —
#' i.e., the AI does not merely re-scale all items uniformly but disfavours
#' (or favours) one group at specific items.
#'
#' @param human_mle Output of \code{\link{simulate_aidif_data}} for
#'   human-scored data.
#' @param ai_mle Output of \code{\link{simulate_aidif_data}} for
#'   AI-scored data.  Must have the same item/group structure.
#' @param fun Scaling function (passed to
#'   the internal scaling function) to use when normalising shifts.
#'   Default: \code{"d_fun3"}.
#'
#' @return A \code{data.frame} with one row per item (per threshold for
#'   polytomous items) and columns:
#' \describe{
#'   \item{\code{shift_g1}}{Scoring shift \eqn{\delta_{i1} = d_{i1}^{AI} - d_{i1}^{H}}.}
#'   \item{\code{shift_g2}}{Scoring shift \eqn{\delta_{i2} = d_{i2}^{AI} - d_{i2}^{H}}.}
#'   \item{\code{DASB}}{Differential AI Scoring Bias: \eqn{\delta_{i2} - \delta_{i1}}.}
#'   \item{\code{se}}{Standard error of DASB under the delta method.}
#'   \item{\code{z}}{Wald z-statistic.}
#'   \item{\code{p_val}}{Two-tailed p-value.}
#' }
#'
#' @details
#' \strong{Estimand.} Define the scoring shift in group \eqn{g} for item
#' \eqn{i} threshold \eqn{j} as:
#' \deqn{\delta_{igj} = d_{igj}^{\text{AI}} - d_{igj}^{\text{Human}}}
#' The DASB is \eqn{\delta_{i2j} - \delta_{i1j}}.  Under
#' \eqn{H_0: \text{DASB}_{ij} = 0} and independence across scoring
#' conditions and groups,
#' \deqn{\widehat{\mathrm{Var}}(\text{DASB}_{ij}) =
#'   \sigma^2_{i1j}^{H} + \sigma^2_{i2j}^{H} +
#'   \sigma^2_{i1j}^{AI} + \sigma^2_{i2j}^{AI}}
#' where each \eqn{\sigma^2} is the diagonal element of the corresponding
#' group-specific covariance matrix.
#'
#' @examples
#' eg <- make_aidif_eg()
#' scoring_bias_test(eg$human, eg$ai)
#'
#' @seealso \code{\link{fit_aidif}}, \code{\link{ai_effect_summary}}
#' @export
# -----------------------------------------------------------------------

scoring_bias_test <- function(human_mle, ai_mle, fun = "d_fun3") {
  
  check_aidif_mle(human_mle, label = "human_mle")
  check_aidif_mle(ai_mle,    label = "ai_mle")
  check_compatible_mles(human_mle, ai_mle)
  check_aidif_fun(fun)
  
  # Extract intercept (d) parameters for each group from each condition.
  # est is a list: group.1, group.2 ; each a data.frame with cols a1, d1, ...
  hg1 <- human_mle$est$group.1
  hg2 <- human_mle$est$group.2
  ag1 <- ai_mle$est$group.1
  ag2 <- ai_mle$est$group.2
  
  n_items      <- nrow(hg1)
  n_thresholds <- ncol(hg1) - 1L   # subtract the slope column
  
  # Variances: diagonal of the group vcov matrix.
  # Parameters are ordered: a1, d1, [d2, ...] per item.
  var_h1 <- Matrix::diag(human_mle$var.cov$group.1)
  var_h2 <- Matrix::diag(human_mle$var.cov$group.2)
  var_a1 <- Matrix::diag(ai_mle$var.cov$group.1)
  var_a2 <- Matrix::diag(ai_mle$var.cov$group.2)
  
  # Build output vectors.
  shift_g1 <- numeric(n_items * n_thresholds)
  shift_g2 <- numeric(n_items * n_thresholds)
  var_shift <- numeric(n_items * n_thresholds)
  row_names <- character(n_items * n_thresholds)
  
  k <- 1L
  for (i in seq_len(n_items)) {
    for (j in seq_len(n_thresholds)) {
      col_name <- paste0("d", j)
      item_name <- rownames(hg1)[i]
      
      # Scoring shift within each group (AI - Human intercept)
      shift_g1[k] <- ag1[i, col_name] - hg1[i, col_name]
      shift_g2[k] <- ag2[i, col_name] - hg2[i, col_name]
      
      # Parameter position in the vcov: (i-1)*(n_thresholds+1) + 1 = slope,
      # +j = threshold j.
      par_pos <- (i - 1L) * (n_thresholds + 1L) + 1L + j
      
      # Variance of DASB under independence of groups and scoring conditions
      var_shift[k] <- var_h1[par_pos] + var_h2[par_pos] +
        var_a1[par_pos] + var_a2[par_pos]
      
      row_names[k] <- if (n_thresholds == 1L) item_name else
        paste0(item_name, "_d", j)
      k <- k + 1L
    }
  }
  
  DASB <- shift_g2 - shift_g1
  se   <- sqrt(var_shift)
  z    <- DASB / se
  p    <- 2 * (1 - pnorm(abs(z)))
  
  out <- data.frame(
    shift_g1 = round(shift_g1, 4),
    shift_g2 = round(shift_g2, 4),
    DASB     = round(DASB, 4),
    se       = round(se, 4),
    z        = round(z, 3),
    p_val    = round(p, 4),
    row.names = row_names
  )
  out
}

# -----------------------------------------------------------------------
#' Summarise the effect of AI scoring on DIF flagging.
#'
#' Compares the DIF flagging patterns from human and AI scoring conditions
#' and classifies each item as: \code{"stable_clean"} (not flagged in
#' either), \code{"stable_dif"} (flagged in both), \code{"introduced"}
#' (flagged only under AI), \code{"masked"} (flagged only under human), or
#' \code{"new_direction"} (flagged in both but bias reverses sign).
#'
#' @param dif_human A \code{data.frame} returned by
#'   \code{\link{fit_aidif}} for the human scoring condition.
#' @param dif_ai A \code{data.frame} returned by
#'   \code{\link{fit_aidif}} for the AI scoring condition.
#' @param alpha Significance threshold for flagging. Default: \code{0.05}.
#'
#' @return A \code{data.frame} with one row per item/threshold and columns:
#' \describe{
#'   \item{\code{human_delta}}{Estimated DIF effect under human scoring.}
#'   \item{\code{ai_delta}}{Estimated DIF effect under AI scoring.}
#'   \item{\code{human_flag}}{Logical: flagged under human scoring?}
#'   \item{\code{ai_flag}}{Logical: flagged under AI scoring?}
#'   \item{\code{status}}{Classification (see Description).}
#' }
#'
#' @examples
#' eg <- make_aidif_eg()
#' mod <- fit_aidif(eg$human, eg$ai)
#' ai_effect_summary(mod$dif_human, mod$dif_ai)
#'
#' @seealso \code{\link{scoring_bias_test}}, \code{\link{fit_aidif}}
#' @export
# -----------------------------------------------------------------------

ai_effect_summary <- function(dif_human, dif_ai, alpha = 0.05) {
  # Accept full model objects as a convenience shorthand
  if (inherits(dif_human, "aidif")) dif_human <- dif_human$dif_human
  if (inherits(dif_ai,    "aidif")) dif_ai    <- dif_ai$dif_ai
  if (inherits(dif_human, "rdif"))  dif_human <- dif_human$dif_test
  if (inherits(dif_ai,    "rdif"))  dif_ai    <- dif_ai$dif_test
  
  # Normalise p-value column: accept both p_val (internal) and p.val (hand-built)
  .pval_col <- function(df, label) {
    if ("p_val" %in% names(df)) return("p_val")
    if ("p.val" %in% names(df)) return("p.val")
    stop(sprintf("`%s` must be a data.frame with columns 'delta' and 'p_val'.", label),
         call. = FALSE)
  }
  if (!is.data.frame(dif_human) || !"delta" %in% names(dif_human))
    stop("`dif_human` must be a data.frame from fit_aidif() or wald_dif_test().", call. = FALSE)
  if (!is.data.frame(dif_ai) || !"delta" %in% names(dif_ai))
    stop("`dif_ai` must be a data.frame from fit_aidif() or wald_dif_test().", call. = FALSE)
  if (nrow(dif_human) != nrow(dif_ai))
    stop("`dif_human` and `dif_ai` must have the same number of rows.", call. = FALSE)
  
  pcol_h <- .pval_col(dif_human, "dif_human")
  pcol_a <- .pval_col(dif_ai,    "dif_ai")
  flag_h <- dif_human[[pcol_h]] < alpha
  flag_a <- dif_ai[[pcol_a]]    < alpha
  
  status <- ifelse(!flag_h & !flag_a, "stable_clean",
                   ifelse( flag_h &  flag_a &
                             sign(dif_human$delta) == sign(dif_ai$delta),
                           "stable_dif",
                           ifelse( flag_h &  flag_a &
                                     sign(dif_human$delta) != sign(dif_ai$delta),
                                   "new_direction",
                                   ifelse(!flag_h &  flag_a, "introduced",
                                          "masked"))))
  
  data.frame(
    human_delta = round(dif_human$delta, 4),
    ai_delta    = round(dif_ai$delta, 4),
    human_flag  = flag_h,
    ai_flag     = flag_a,
    status      = status,
    row.names   = rownames(dif_human)
  )
}

# -----------------------------------------------------------------------
#' Anchor item weights from the robust AI-DIF procedure.
#'
#' Returns the bi-square weights assigned to each item under each scoring
#' condition.  Items with weight near zero are effectively excluded from
#' the robust scaling estimate, indicating likely DIF contamination.
#'
#' @param object An \code{aidif} object from \code{\link{fit_aidif}}.
#'
#' @return A \code{data.frame} with columns \code{human_weight} and
#'   (if AI data were provided) \code{ai_weight}.  Higher weight means
#'   the item is contributing more to the robust scale estimate.
#'
#' @examples
#' eg <- make_aidif_eg()
#' mod <- fit_aidif(eg$human, eg$ai)
#' anchor_weights(mod)
#'
#' @export
# -----------------------------------------------------------------------

anchor_weights <- function(object) {
  if (!inherits(object, "aidif"))
    stop("`object` must be of class 'aidif'.", call. = FALSE)
  
  w_human <- object$human_fit$weights
  item_names <- rownames(object$dif_human)
  
  if (!is.null(object$ai_fit)) {
    w_ai <- object$ai_fit$weights
    out  <- data.frame(
      human_weight = round(w_human, 4),
      ai_weight    = round(w_ai, 4),
      row.names    = item_names
    )
  } else {
    out <- data.frame(
      human_weight = round(w_human, 4),
      row.names    = item_names
    )
  }
  out
}