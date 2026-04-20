## ============================================================
## RUN THIS ONCE in RStudio before devtools::document()
## Set your working directory to the package root, then:
##   source("data-raw/make_aidif_eg.R")
## ============================================================

set.seed(2024)
n_items <- 6L; se_a <- 0.06; se_d <- 0.07
a_true  <- c(1.2, 0.9, 1.4, 1.0, 1.3, 0.8)
d_true  <- c(0.3, -0.5, 0.1, -0.2, 0.4, -0.1)
impact  <- 0.5

make_est_df <- function(a, d) {
  df <- data.frame(a1 = round(a, 4), d1 = round(d, 4))
  rownames(df) <- paste0("item", seq_along(a)); df
}
make_vcov <- function(n) {
  v <- diag(c(rbind(se_a^2, se_d^2)), n * 2L)
  for (i in seq_len(n)) {
    ia <- (i-1L)*2L+1L; id <- ia+1L
    v[ia,id] <- v[id,ia] <- -0.3*se_a*se_d
  }; v
}

a_h1 <- a_true + c( .05,-.03, .04,-.02, .03,-.04)
d_h1 <- d_true + c( .04,-.06, .03, .05,-.04, .02)
a_h2 <- a_true + c(-.03, .04,-.02, .05,-.03, .02)
d_h2 <- d_true - impact + c(.03,-.05,.04,-.02,.06,-.03)
d_h2[1] <- d_h2[1] + 0.5   # Item 1 human DIF

a_a1 <- a_h1 + c( .02,-.01,.03,-.01,.02,-.02)
d_a1 <- d_h1 + 0.1 + c(.03,-.02,.01,.02,-.03,.01)
a_a2 <- a_h2 + c(-.02,.03,-.01,.02,-.02,.01)
d_a2 <- d_h2 + 0.1 + c(.02,-.03,.04,-.01,.03,-.02)
d_a2[3] <- d_a2[3] + 0.4   # Item 3 DASB

par_names <- c(rbind(paste0("item",1:n_items,".a1"),
                     paste0("item",1:n_items,".d1")))
vcov_r <- make_vcov(n_items)

aidif.eg <- list(
  human = list(
    par.names = list(internal=par_names, original=par_names),
    est     = list(group.1=make_est_df(a_h1,d_h1), group.2=make_est_df(a_h2,d_h2)),
    var.cov = list(group.1=vcov_r, group.2=vcov_r)
  ),
  ai = list(
    par.names = list(internal=par_names, original=par_names),
    est     = list(group.1=make_est_df(a_a1,d_a1), group.2=make_est_df(a_a2,d_a2)),
    var.cov = list(group.1=vcov_r, group.2=vcov_r)
  )
)

if (!dir.exists("data")) dir.create("data")
save(aidif.eg, file="data/aidif.eg.rda", compress="bzip2")
cat("Saved: data/aidif.eg.rda\n")

desc <- readLines("DESCRIPTION")
desc <- gsub("^LazyData: false", "LazyData: true", desc)
writeLines(desc, "DESCRIPTION")
cat("DESCRIPTION: LazyData set to true\n")
cat("\nNow run:\n  devtools::document()\n  devtools::install()\n")
