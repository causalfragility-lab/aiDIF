# Internal validation helpers.

check_aidif_mle <- function(mle, label = "mle") {
  if (!is.list(mle))
    stop(sprintf("`%s` must be a list (structured mle list (see simulate_aidif_data())).", label),
         call. = FALSE)
  if (!all(c("est", "var.cov") %in% names(mle)))
    stop(sprintf("`%s` must contain `est` and `var.cov` components.", label),
         call. = FALSE)
  if (!is.list(mle$est) || !is.list(mle$var.cov) || length(mle$est) < 2L)
    stop(sprintf("`%s$est` and `%s$var.cov` must each be lists with >= 2 groups.",
                 label, label), call. = FALSE)
  if (length(mle$est) != length(mle$var.cov))
    stop(sprintf("`%s$est` and `%s$var.cov` must have the same number of groups.",
                 label, label), call. = FALSE)
  invisible(TRUE)
}

check_compatible_mles <- function(human_mle, ai_mle) {
  n_groups_h <- length(human_mle$est)
  n_groups_a <- length(ai_mle$est)
  if (n_groups_h != n_groups_a)
    stop(sprintf("`human_mle` has %d groups but `ai_mle` has %d groups.",
                 n_groups_h, n_groups_a), call. = FALSE)

  n_items_h <- nrow(human_mle$est[[1]])
  n_items_a <- nrow(ai_mle$est[[1]])
  if (n_items_h != n_items_a)
    stop(sprintf("`human_mle` has %d items but `ai_mle` has %d items.",
                 n_items_h, n_items_a), call. = FALSE)

  n_par_h <- ncol(human_mle$est[[1]])
  n_par_a <- ncol(ai_mle$est[[1]])
  if (n_par_h != n_par_a)
    stop("Item parameter structure (number of columns in `est`) differs between `human_mle` and `ai_mle`.",
         call. = FALSE)

  invisible(TRUE)
}

check_aidif_fun <- function(fun) {
  allowed <- c("a_fun1", "a_fun2", "d_fun1", "d_fun2", "d_fun3")
  if (!is.character(fun) || length(fun) != 1L || !fun %in% allowed)
    stop(sprintf("`fun` must be one of: %s.", paste(allowed, collapse = ", ")),
         call. = FALSE)
  invisible(TRUE)
}

check_aidif_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("`alpha` must be a single numeric value in (0, 1).", call. = FALSE)
  invisible(TRUE)
}
