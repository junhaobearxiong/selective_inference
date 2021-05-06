# simulation for normals with correlated errors across features
# p-values can be computed for both F-test and jackstraw
# as well as the KS-test based joint null criterion described in the jackstraw paper
# remove comments for the desired portion to generate plots 

library(MASS)
library(ggplot2)
library(jackstraw)
library(gridExtra)

# construct covariance matrix 
gen_cov <- function(ngenes, pct.corr, corr=0.1) {
  Sigma <- diag(ngenes)
  Sigma[1:round(pct.corr * ngenes), 1:round(pct.corr * ngenes)] <- corr
  diag(Sigma) <- 1
  return(Sigma)
}

# generate expression data from normal distribution
gen_data <- function(ngenes, pct.corr, ncells=500, corr=0.1) {
  if (pct.corr == 0) {
    expr <- matrix(rnorm(ncells * ngenes), ncol = ngenes)
  } else {
    expr <- mvrnorm(n = ncells, mu = rep(0, ngenes), Sigma = gen_cov(ngenes, pct.corr, corr))
  } 
  return(expr)
}

# get F-test p-values distribution under the null (across all genes)
get_ftest_pval_dist <- function(expr) {
  ngenes <- ncol(expr)
  pvals <- matrix(NA, nrow = ngenes)
  pt <- as.vector(svd(expr, nu = 1, nv = 1)$u) # pseudotime
  for (i in 1:ngenes) {
    df <- data.frame(expr = expr[, i], pt = pt)
    # F-test 
    mod <- lm(expr ~ pt, data = df)
    pvals[i] <- summary(mod)$coefficients[2, 4]
  }
  return(pvals)
}

# get jackstraw p-values distribution under the null 
get_jackstraw_pval_dist <- function(expr) {
  # input for jackstraw has cells on columns and genes on rows
  test <- jackstraw_pca(dat = t(expr), r = 1)
  return(test$p.value)
}

# get p-values for KS test (across replications of simulated data)
get_kstest_pval <- function(test, ngenes, pct.corr, nreps=500, ncells=500, corr=0.1) {
  ks.pval <- matrix(NA, nrow = nreps)
  for (i in 1:nreps) {
    if (i %% 10 == 0) {
      print(i)
    }
    if (test == 'ftest') {
      pval <- get_ftest_pval_dist(ngenes = ngenes, pct.corr = pct.corr, ncells = ncells, corr = corr)
    } else if (test == 'jackstraw') {
      pval <- get_jackstraw_pval_dist(ngenes = ngenes, pct.corr = pct.corr, ncells = ncells, corr = corr)
    }
    kstest <- ks.test(x = pval, y = punif)
    ks.pval[i] <- kstest$p.value
  }
  return(ks.pval)
}

# pval <- get_jackstraw_pval_dist(ngenes = 1000, pct.corr = 0.5)
# pval <- get_kstest_pval(ngenes = 100, pct.corr = 0.0)
# pval1 <- get_ftest_pval_dist(ngenes = 100, pct.corr = 0.2)
# pval2 <- get_jackstraw_pval_dist(ngenes = 100, pct.corr = 0.2)
# ggplot() +
#   stat_qq(aes(sample = pval1), distribution = stats::qunif, colour = "green") + 
#   stat_qq(aes(sample = pval2), distribution = stats::qunif, colour = "red") + 
#   geom_abline(aes(slope = 1, intercept = 0), linetype = 2) + 
#   labs(title = 'combined')

#############################p-value over genes (F-test)##############################################
# ngenes.list <- c(100, 500, 1000, 2000)
# pct.corr.list <- c(0.0, 0.2, 0.5)
# 
# plots.ftest <- list()
# count <- 1
# for (pct.corr in pct.corr.list) {
#   for (ngenes in ngenes.list) {
#     print(paste(ngenes, pct.corr))
#     expr <- gen_data(ngenes = ngenes, pct.corr = pct.corr)
#     # F-test
#     pval.ftest <- get_ftest_pval_dist(expr)
#     print("f-test done")
#     df <- data.frame(ftest = pval.ftest)
#     plots.ftest[[count]] <- ggplot(data = df) +
#       stat_qq(aes(sample = ftest), distribution = stats::qunif) + 
#       geom_abline(aes(slope = 1, intercept = 0), color = "red") + 
#       labs(title = paste("ngenes", ngenes,", pct.corr", pct.corr))
#     count <- count + 1
#   }
# }
# pl.ftest <- do.call("grid.arrange", c(plots.ftest, ncol=4))
# ggsave(file="figures/corr_error_ftest.png", pl.ftest, width = 12, height = 9)

#############################p-value over genes (F-test + jackstraw)##############################################
# ngenes.list <- c(100, 500, 1000, 2000)
# pct.corr.list <- c(0.0, 0.2, 0.5)
# 
# plots.combined <- list()
# count <- 1
# for (pct.corr in pct.corr.list) {
#   for (ngenes in ngenes.list) {
#     print(paste(ngenes, pct.corr))
#     expr <- gen_data(ngenes = ngenes, pct.corr = pct.corr)
#     # F-test
#     pval.ftest <- get_ftest_pval_dist(expr)
#     print("f-test done")
#     # jackstraw
#     pval.js <- get_jackstraw_pval_dist(expr)
#     print("jackstraw done")
#     df <- data.frame(ftest = pval.ftest, js = pval.js)
#     plots.combined[[count]] <- ggplot(data = df) +
#       stat_qq(aes(sample = ftest), distribution = stats::qunif, colour = "red") + 
#       stat_qq(aes(sample = js), distribution = stats::qunif, colour = "green") + 
#       geom_abline(aes(slope = 1, intercept = 0), linetype = 2) + 
#       labs(title = paste("ngenes", ngenes,", pct.corr", pct.corr))
#     count <- count + 1
#   }
# }
# pl.combined <- do.call("grid.arrange", c(plots.combined, ncol=4))
# ggsave(file="figures/corr_error_combined.png", pl.combined, width = 12, height = 9)

##############################p-value over replications (KS test)########################################
# ngenes.list <- c(100, 500)
# pct.corr.list <- c(0.0, 0.2, 0.5)
# 
# plots.ks.js <- list()
# count <- 1
# for (pct.corr in pct.corr.list) {
#   for (ngenes in ngenes.list) {
#     print(paste(ngenes, pct.corr))
#     # jackstraw
#     pval.js <- get_kstest_pval(test = 'jackstraw', ngenes = ngenes, pct.corr = pct.corr, nreps = 100)
#     plots.ks.js[[count]] <- ggplot(data.frame(pval = pval.js), aes(sample = pval)) +
#       geom_abline(intercept = 0, slope = 1, alpha = 0.5, color='red') +
#       stat_qq(distribution = stats::qunif) + 
#       labs(title = paste("ngenes", ngenes,", pct.corr", pct.corr))
#     print("jackstraw done")
#     count <- count + 1
#   }
# }
# pl.ks.js <- do.call("grid.arrange", c(plots.ks.js, ncol=2))
# ggsave(file="figures/corr_error_ks_jackstraw.png", pl.ks.js, width = 6, height = 9)

