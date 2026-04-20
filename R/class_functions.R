# -----------------------------------------------------------------------
#' S3 print method for class \code{"aidif"}.
#'
#' Prints a compact summary of the estimated robust scaling parameters and,
#' when available, the number of items flagged for DIF and DASB.
#'
#' @param x   An object of class \code{"aidif"}.
#' @param ... Further arguments (currently ignored).
#' @return \code{x}, invisibly.
#' @export
print.aidif <- function(x, ...) {
  if (!inherits(x, "aidif"))
    stop("Object is not of class 'aidif'.", call. = FALSE)
  
  cat("AI-DIF Analysis\n")
  cat(rep("-", 40), "\n", sep = "")
  
  # Human scoring
  h_est <- round(x$human_fit$est, 4)
  h_se  <- round(x$human_fit$dtf_test$robust_se, 4)
  cat(sprintf("Human scoring  \u2014 robust scale est: %6.4f  (SE: %6.4f)\n", h_est, h_se))
  
  n_items <- nrow(x$dif_human)
  n_dif_h <- sum(x$dif_human$p_val < x$alpha)
  cat(sprintf("               \u2014 DIF items flagged: %d / %d\n", n_dif_h, n_items))
  
  if (!is.null(x$ai_fit)) {
    a_est <- round(x$ai_fit$est, 4)
    a_se  <- round(x$ai_fit$dtf_test$robust_se, 4)
    cat(sprintf("AI scoring     \u2014 robust scale est: %6.4f  (SE: %6.4f)\n", a_est, a_se))
    n_dif_a <- sum(x$dif_ai$p_val < x$alpha)
    cat(sprintf("               \u2014 DIF items flagged: %d / %d\n", n_dif_a, n_items))
    
    n_dasb <- sum(x$scoring_bias$p_val < x$alpha)
    cat(sprintf("DASB test      \u2014 items with differential AI bias: %d / %d\n",
                n_dasb, n_items))
  }
  invisible(x)
}

# -----------------------------------------------------------------------
#' S3 summary method for class \code{"aidif"}.
#'
#' Prints a detailed report including DIF test tables for each scoring
#' condition, the DASB table, and the AI-effect classification.
#'
#' @param object An object of class \code{"aidif"}.
#' @param ... Further arguments (currently ignored).
#' @return \code{NULL}, invisibly.
#' @export
summary.aidif <- function(object, ...) {
  if (!inherits(object, "aidif"))
    stop("Object is not of class 'aidif'.", call. = FALSE)
  
  n_items <- nrow(object$dif_human)
  fun     <- object$fun
  alpha   <- object$alpha
  
  cat("=============================================================\n")
  cat(" AI Differential Item Functioning Analysis (aiDIF)\n")
  cat("=============================================================\n\n")
  
  cat(sprintf("Items: %d   Scaling function: %s   Alpha: %.2f\n\n",
              n_items, fun, alpha))
  
  # ---- Human scoring section -------------------------------------------
  cat("--- Human Scoring DIF ----------------------------------------\n")
  h <- object$human_fit
  cat(sprintf("  Robust scale estimate:  %.4f  (SE: %.4f)\n",
              h$est, h$dtf_test$robust_se))
  cat(sprintf("  Iterations: %d   Multiple solutions: %s\n\n",
              h$n.iter, h$multiple.solutions))
  cat("  Wald DIF tests:\n")
  print(round(object$dif_human, 4))
  cat("\n")
  
  # ---- AI scoring section ----------------------------------------------
  if (!is.null(object$ai_fit)) {
    cat("--- AI Scoring DIF -------------------------------------------\n")
    a <- object$ai_fit
    cat(sprintf("  Robust scale estimate:  %.4f  (SE: %.4f)\n",
                a$est, a$dtf_test$robust_se))
    cat(sprintf("  Iterations: %d   Multiple solutions: %s\n\n",
                a$n.iter, a$multiple.solutions))
    cat("  Wald DIF tests:\n")
    print(round(object$dif_ai, 4))
    cat("\n")
    
    # ---- DASB table ------------------------------------------------------
    cat("--- Differential AI Scoring Bias (DASB) ---------------------\n")
    cat("  H0: AI scoring shift does not differ across groups\n")
    cat("  (Positive DASB => AI scoring disadvantages focal group)\n\n")
    print(object$scoring_bias)
    cat("\n")
    
    # ---- AI-effect summary -----------------------------------------------
    cat("--- AI-Effect Classification ---------------------------------\n")
    cat("  stable_clean  : not flagged in either condition\n")
    cat("  stable_dif    : flagged in both (same direction)\n")
    cat("  introduced    : flagged only under AI scoring\n")
    cat("  masked        : flagged only under human scoring\n")
    cat("  new_direction : flagged in both, opposite direction\n\n")
    print(object$ai_effect)
    cat("\n")
    
    # ---- Summary counts --------------------------------------------------
    tbl <- table(object$ai_effect$status)
    cat("  Status counts:\n")
    print(tbl)
  }
  
  invisible(NULL)
}

# -----------------------------------------------------------------------
#' S3 plot method for class \code{"aidif"}.
#'
#' Produces one of several diagnostic plots depending on \code{type}.
#'
#' @param x    An object of class \code{"aidif"}.
#' @param type Character. One of:
#'   \describe{
#'     \item{\code{"dif_forest"}}{Forest plot of DIF estimates with 95\%
#'       confidence intervals for both scoring conditions (default).}
#'     \item{\code{"dasb"}}{Bar chart of DASB estimates with error bars.}
#'     \item{\code{"weights"}}{Dot plot of bi-square anchor weights.}
#'     \item{\code{"rho"}}{Bi-square objective function for human scoring.}
#'   }
#' @param ... Additional graphical parameters passed to low-level plot
#'   functions.
#' @return \code{x}, invisibly.
#' @export
plot.aidif <- function(x, type = "dif_forest", ...) {
  if (!inherits(x, "aidif"))
    stop("Object is not of class 'aidif'.", call. = FALSE)
  
  type <- match.arg(type, c("dif_forest", "dasb", "weights", "rho"))
  
  switch(type,
         dif_forest = .plot_dif_forest(x, ...),
         dasb       = .plot_dasb(x, ...),
         weights    = .plot_weights(x, ...),
         rho        = .plot_rho(x, ...)
  )
  invisible(x)
}

# -----------------------------------------------------------------------
# Internal plot helpers
# -----------------------------------------------------------------------

.plot_dif_forest <- function(x, ...) {
  dif_h <- x$dif_human
  n     <- nrow(dif_h)
  items <- rownames(dif_h)
  z_crit <- qnorm(1 - x$alpha / 2)
  
  has_ai <- !is.null(x$dif_ai)
  dif_a  <- x$dif_ai
  
  # Compute CI limits
  lo_h <- dif_h$delta - z_crit * dif_h$se
  hi_h <- dif_h$delta + z_crit * dif_h$se
  
  xlim <- if (has_ai) {
    lo_a <- dif_a$delta - z_crit * dif_a$se
    hi_a <- dif_a$delta + z_crit * dif_a$se
    range(c(lo_h, hi_h, lo_a, hi_a), na.rm = TRUE)
  } else {
    range(c(lo_h, hi_h), na.rm = TRUE)
  }
  xlim <- xlim + c(-0.1, 0.1) * diff(xlim)
  
  old_par <- par(mar = c(4, 6, 3, 2))
  on.exit(par(old_par))
  
  y_h <- if (has_ai) seq(n, 1) + 0.15 else seq(n, 1)
  y_a <- seq(n, 1) - 0.15
  
  plot(dif_h$delta, y_h,
       xlim = xlim, ylim = c(0.5, n + 0.5),
       yaxt = "n", xlab = "DIF estimate", ylab = "",
       pch = 16, col = "#2166ac",
       main = "DIF Forest Plot: Human vs AI Scoring",
       ...)
  axis(2, at = seq(n, 1), labels = items, las = 1, cex.axis = 0.85)
  abline(v = 0, lty = 2, col = "grey50")
  
  # CIs for human
  segments(lo_h, y_h, hi_h, y_h, col = "#2166ac", lwd = 2)
  
  if (has_ai) {
    points(dif_a$delta, y_a, pch = 17, col = "#d6604d")
    segments(lo_a, y_a, hi_a, y_a, col = "#d6604d", lwd = 2)
    legend("topright", legend = c("Human", "AI"),
           pch = c(16, 17), col = c("#2166ac", "#d6604d"),
           lwd = 2, bty = "n", cex = 0.9)
  }
}

.plot_dasb <- function(x, ...) {
  if (is.null(x$scoring_bias))
    stop("No DASB results available. Re-run fit_aidif() with `ai_mle` provided.",
         call. = FALSE)
  
  sb    <- x$scoring_bias
  n     <- nrow(sb)
  items <- rownames(sb)
  z_crit <- qnorm(1 - x$alpha / 2)
  lo    <- sb$DASB - z_crit * sb$se
  hi    <- sb$DASB + z_crit * sb$se
  sig   <- sb$p_val < x$alpha
  col_  <- ifelse(sig, "#d6604d", "#4393c3")
  
  old_par <- par(mar = c(4, 6, 3, 2))
  on.exit(par(old_par))
  
  barplot(sb$DASB, horiz = TRUE, names.arg = items, las = 1,
          col = col_, border = NA,
          xlab = "DASB (group 2 minus group 1 AI shift)",
          main = "Differential AI Scoring Bias (DASB)",
          ...)
  abline(v = 0, col = "grey30")
  
  legend("topright",
         legend = c(sprintf("Significant (alpha=%.2f)", x$alpha), "Non-significant"),
         fill   = c("#d6604d", "#4393c3"),
         bty    = "n", cex = 0.9)
}

.plot_weights <- function(x, ...) {
  w <- anchor_weights(x)
  n     <- nrow(w)
  items <- rownames(w)
  
  old_par <- par(mar = c(4, 6, 3, 2))
  on.exit(par(old_par))
  
  plot(w$human_weight, seq_len(n),
       xlim = c(0, 1), ylim = c(0.5, n + 0.5),
       yaxt = "n", xlab = "Bi-square weight", ylab = "",
       pch  = 16, col = "#2166ac",
       main = "Anchor Weights: Human vs AI Scoring",
       ...)
  axis(2, at = seq_len(n), labels = items, las = 1, cex.axis = 0.85)
  abline(v = c(0, 1), lty = 3, col = "grey70")
  
  if (!is.null(x$ai_fit)) {
    points(w$ai_weight, seq_len(n) - 0.2, pch = 17, col = "#d6604d")
    legend("bottomright", legend = c("Human", "AI"),
           pch = c(16, 17), col = c("#2166ac", "#d6604d"),
           bty = "n", cex = 0.9)
  }
}

.plot_rho <- function(x, ...) {
  rp <- x$human_fit$rho.plot
  if (is.null(rp))
    stop("No rho grid available in `human_fit`.", call. = FALSE)
  
  plot(rp$theta, rp$rho, type = "l",
       xlab = expression(theta), ylab = expression(rho(theta)),
       main = "Bi-square Objective (Human Scoring)",
       col  = "#2166ac", lwd = 2, ...)
  abline(v = x$human_fit$est, lty = 2, col = "#d6604d")
  legend("topright",
         legend = sprintf("Est = %.4f", x$human_fit$est),
         lty = 2, col = "#d6604d", bty = "n", cex = 0.9)
}

# S3 methods are registered via NAMESPACE (S3method entries) and
# the @export tags on each function above \u2014 no .onLoad hook needed.
#' @importFrom graphics abline axis barplot legend par points segments
NULL