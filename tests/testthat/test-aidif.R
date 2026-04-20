library(testthat)
library(aiDIF)

# -----------------------------------------------------------------------
# Helper: build a minimal valid mle list (2 items, 2 groups, 2PL)
# -----------------------------------------------------------------------

make_mle <- function(d_shift_g2 = 0, n_items = 4L, se_a = 0.06, se_d = 0.07) {
  set.seed(1)
  a <- runif(n_items, 0.8, 1.5)
  d <- rnorm(n_items, 0, 0.5)

  make_df <- function(av, dv) {
    df <- data.frame(a1 = av, d1 = dv)
    rownames(df) <- paste0("item", seq_len(n_items))
    df
  }
  make_vc <- function() {
    v <- diag(c(rbind(se_a^2, se_d^2)), n_items * 2)
    for (i in seq_len(n_items)) {
      ia <- (i-1)*2 + 1; id <- ia + 1
      v[ia, id] <- v[id, ia] <- -0.3 * se_a * se_d
    }
    v
  }

  d2 <- d - 0.5 + d_shift_g2 * c(1, rep(0, n_items - 1))
  list(
    par.names = list(
      internal = c(rbind(paste0("item", 1:n_items, ".a1"), paste0("item", 1:n_items, ".d1"))),
      original = c(rbind(paste0("item", 1:n_items, ".a1"), paste0("item", 1:n_items, ".d1")))
    ),
    est = list(group.1 = make_df(a, d), group.2 = make_df(a, d2)),
    var.cov = list(group.1 = make_vc(), group.2 = make_vc())
  )
}

# -----------------------------------------------------------------------
# Input validation
# -----------------------------------------------------------------------

test_that("check_aidif_mle rejects non-list", {
  expect_error(check_aidif_mle("not_a_list"), "must be a list")
})

test_that("check_aidif_mle rejects missing est", {
  bad <- list(est = list(g1 = data.frame(a1=1, d1=0)))
  expect_error(check_aidif_mle(bad), "var.cov")
})

test_that("check_compatible_mles detects item mismatch", {
  m4 <- make_mle(n_items = 4)
  m6 <- make_mle(n_items = 6)
  expect_error(check_compatible_mles(m4, m6), "items")
})

test_that("read_ai_scored returns aidif_data class", {
  m <- make_mle()
  out <- read_ai_scored(m, m)
  expect_s3_class(out, "aidif_data")
  expect_named(out, c("human", "ai"))
})

# -----------------------------------------------------------------------
# scoring_bias_test
# -----------------------------------------------------------------------

test_that("scoring_bias_test returns correct structure", {
  human <- make_mle()
  ai    <- make_mle(d_shift_g2 = 0.4)   # DASB at item 1
  out   <- scoring_bias_test(human, ai)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("shift_g1", "shift_g2", "DASB", "se", "z", "p_val"))
  expect_equal(nrow(out), 4L)   # 4 items
})

test_that("scoring_bias_test detects planted DASB", {
  human <- make_mle(d_shift_g2 = 0)
  ai    <- make_mle(d_shift_g2 = 1.5)   # large planted DASB at item 1
  out   <- scoring_bias_test(human, ai)
  # Item 1 should be significant
  expect_true(out["item1", "p_val"] < 0.05)
})

test_that("scoring_bias_test DASB is near zero when no differential bias", {
  human <- make_mle()
  # AI adds same drift to both groups: DASB should be near 0
  ai <- human
  ai$est$group.1$d1 <- human$est$group.1$d1 + 0.1
  ai$est$group.2$d1 <- human$est$group.2$d1 + 0.1
  out <- scoring_bias_test(human, ai)
  expect_true(all(abs(out$DASB) < 0.01))
})

# -----------------------------------------------------------------------
# ai_effect_summary
# -----------------------------------------------------------------------

test_that("ai_effect_summary classifies correctly", {
  dif_human <- data.frame(delta = c(0.5, 0.0, -0.4, 0.0),
                          se    = c(0.1, 0.1,  0.1, 0.1),
                          p.val = c(0.01, 0.8, 0.02, 0.9))
  dif_ai    <- data.frame(delta = c(0.5, 0.4,  0.0, 0.0),
                          se    = c(0.1, 0.1,  0.1, 0.1),
                          p.val = c(0.01, 0.02, 0.8, 0.9))
  rownames(dif_human) <- rownames(dif_ai) <- paste0("item", 1:4)

  out <- ai_effect_summary(dif_human, dif_ai)
  expect_equal(out["item1", "status"], "stable_dif")
  expect_equal(out["item2", "status"], "introduced")
  expect_equal(out["item3", "status"], "masked")
  expect_equal(out["item4", "status"], "stable_clean")
})

# -----------------------------------------------------------------------
# fit_aidif
# -----------------------------------------------------------------------

test_that("fit_aidif runs with human_mle only", {
  human <- make_mle(d_shift_g2 = 0.5)
  mod   <- fit_aidif(human_mle = human)
  expect_s3_class(mod, "aidif")
  expect_null(mod$ai_fit)
  expect_null(mod$scoring_bias)
  expect_s3_class(mod$human_fit, "rdif")
})

test_that("fit_aidif runs with both scoring conditions", {
  human <- make_mle(d_shift_g2 = 0.5)
  ai    <- make_mle(d_shift_g2 = 0.8)
  mod   <- fit_aidif(human_mle = human, ai_mle = ai)
  expect_s3_class(mod, "aidif")
  expect_s3_class(mod$ai_fit, "rdif")
  expect_s3_class(mod$scoring_bias, "data.frame")
  expect_s3_class(mod$ai_effect, "data.frame")
})

test_that("fit_aidif rejects incompatible mles", {
  m4 <- make_mle(n_items = 4)
  m6 <- make_mle(n_items = 6)
  expect_error(fit_aidif(m4, m6), "items")
})

# -----------------------------------------------------------------------
# simulate_aidif_data
# -----------------------------------------------------------------------

test_that("simulate_aidif_data returns valid structure", {
  dat <- simulate_aidif_data(n_items = 6, seed = 99)
  expect_named(dat, c("human", "ai"))
  expect_named(dat$human$est, c("group.1", "group.2"))
  expect_equal(nrow(dat$human$est$group.1), 6L)
})

test_that("simulate_aidif_data can be passed to fit_aidif", {
  dat <- simulate_aidif_data(n_items = 5, dif_items = 1, dasb_items = 3, seed = 7)
  expect_no_error(fit_aidif(dat$human, dat$ai))
})

# -----------------------------------------------------------------------
# anchor_weights
# -----------------------------------------------------------------------

test_that("anchor_weights returns correct columns with ai data", {
  human <- make_mle(d_shift_g2 = 0.5)
  ai    <- make_mle(d_shift_g2 = 0.5)
  mod   <- fit_aidif(human, ai)
  w     <- anchor_weights(mod)
  expect_named(w, c("human_weight", "ai_weight"))
  expect_true(all(w$human_weight >= 0 & w$human_weight <= 1))
})

test_that("anchor_weights errors on non-aidif object", {
  expect_error(anchor_weights(list()), "class 'aidif'")
})
