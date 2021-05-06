# simulation with cell maturity variable and correlated error
# p-values can be computed for both F-test and jackstraw
# F-test can use both PC1 and slingshot as pseudotime 
# remove comments for the desired portion to generate plots 

library(MASS)
library(ggplot2)
library(jackstraw)
library(gridExtra)
library(slingshot)
library(SingleCellExperiment)

# construct covariance matrix 
gen_cov <- function(ngenes, pct.corr, corr=0.1) {
  Sigma <- diag(ngenes)
  Sigma[1:round(pct.corr * ngenes), 1:round(pct.corr * ngenes)] <- corr
  diag(Sigma) <- 1
  return(Sigma)
}

# generate expression data from normal distribution
gen_data <- function(ngenes, pct.corr, pct.de, ncells=500, corr=0.1, maturity.effect=0.05) {
  if (pct.corr == 0) {
    expr <- matrix(rnorm(ncells * ngenes), ncol = ngenes)
  } else {
    expr <- mvrnorm(n = ncells, mu = rep(0, ngenes), Sigma = gen_cov(ngenes, pct.corr, corr))
  } 
  # indices of DE genes 
  genes.de <- sort(sample(c(1:ngenes), size = pct.de * ngenes))
  # indices of non-DE genes
  genes.nonde <- setdiff(c(1:ngenes), genes.de)
  # true "maturity" of cells (should be captured by pseudotime) 
  maturity <- rep(1:4, each = ncells / 4)
  # DE genes is affected by maturity 
  expr[, genes.de] <- expr[, genes.de] + maturity * maturity.effect
  # center expression matrix
  expr <- scale(expr, center = TRUE, scale = FALSE)
  return(list("expr" = expr, "genes.de" = genes.de, "genes.nonde" = genes.nonde))
}

# get pseudotime from slingshot
get_slingshot_pseudotime <- function(expr) {
  sce <- SingleCellExperiment(list(expr=t(expr)))
  pcs <- prcomp(expr)$x[, 1:5]
  reducedDims(sce) <- SimpleList(PCA = pcs)
  colData(sce)$cluster <- 1  # assume a single lineage
  sce <- slingshot(sce, reducedDim = "PCA", clusterLabels = "cluster")
  return(sce$slingPseudotime_1)
}

# get F-test p-values of each gene, given expression matrix
get_ftest_pval_dist <- function(expr, pt.method='pca') {
  ngenes <- ncol(expr)
  pvals <- matrix(NA, nrow = ngenes)
  # get pseudotime
  if (pt.method == 'pca') {
    pt <- as.vector(svd(expr, nu = 1, nv = 1)$u) 
  } else if (pt.method == 'slingshot') {
    pt <- get_slingshot_pseudotime(expr)
  }
  for (i in 1:ngenes) {
    df <- data.frame(expr = expr[, i], pt = pt)
    # F-test 
    mod <- lm(expr ~ pt, data = df)
    pvals[i] <- summary(mod)$coefficients[2, 4]
  }
  return(pvals)
}

# get jackstraw p-values of each gene, given expression matrix
get_jackstraw_pval_dist <- function(expr) {
  # input for jackstraw has cells on columns and genes on rows
  test <- jackstraw_pca(dat = t(expr), r = 1)
  return(test$p.value)
}

# pval1 <- get_ftest_pval_dist(data$expr)
# pval2 <- get_jackstraw_pval_dist(data$expr)
# 
# df <- data.frame(ftest = pval1[data$genes.nonde], js = pval2[data$genes.nonde])
# ggplot(data = df) +
#   stat_qq(aes(sample = ftest), distribution = stats::qunif, colour = "red") + 
#   stat_qq(aes(sample = js), distribution = stats::qunif, colour = "green") + 
#   geom_abline(aes(slope = 1, intercept = 0), linetype = 2) + 
#   labs(title = "")

################ pct.de vs. maturity effect (F-test: pca / slingshot) ########################

# pct.de.list <- c(0.2, 0.5, 0.8)
# maturity.effect.list <- c(0.05, 0.2, 1)
# plots <- list()
# count <- 1
# for (pct.de in pct.de.list) {
#   for (maturity.effect in maturity.effect.list) {
#     print(paste(pct.de, maturity.effect))
#     data <- gen_data(ngenes = 500, pct.corr = 0.3, pct.de = pct.de, maturity.effect = maturity.effect)
#     pval1 <- get_ftest_pval_dist(data$expr, pt.method = "pca")
#     pval2 <- get_ftest_pval_dist(data$expr, pt.method = "slingshot")
#     # only plot the p-values for non-de genes
#     df <- data.frame(pval1 = pval1[data$genes.nonde], pval2 = pval2[data$genes.nonde])
#     plots[[count]] <- ggplot(data = df) +
#       stat_qq(aes(sample = pval1), distribution = stats::qunif, colour = "red") +
#       stat_qq(aes(sample = pval2), distribution = stats::qunif, colour = "green") +
#       geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
#       labs(title = paste("pct.de", pct.de,", maturity.effect", maturity.effect))
#     count <- count + 1
#   }
# }
# pl <- do.call("grid.arrange", c(plots, ncol=3))
# ggsave(file="figures/corr_error_de_pctde_mateffect_slingshot.png", pl, width = 10, height = 9)

################ pct.de vs. maturity effect (F-test / js) ########################

# pct.de.list <- c(0.2, 0.5, 0.8)
# maturity.effect.list <- c(0.05, 0.2, 1)
# plots <- list()
# count <- 1
# for (pct.de in pct.de.list) {
#   for (maturity.effect in maturity.effect.list) {
#     print(paste(pct.de, maturity.effect))
#     data <- gen_data(ngenes = 500, pct.corr = 0.3, pct.de = pct.de, maturity.effect = maturity.effect)
#     pval.ftest <- get_ftest_pval_dist(data$expr)
#     pval.js <- get_jackstraw_pval_dist(data$expr)
#     # only plot the p-values for non-de genes
#     df <- data.frame(ftest = pval.ftest[data$genes.nonde], js = pval.js[data$genes.nonde])
#     plots[[count]] <- ggplot(data = df) +
#       stat_qq(aes(sample = ftest), distribution = stats::qunif, colour = "red") + 
#       stat_qq(aes(sample = js), distribution = stats::qunif, colour = "green") + 
#       geom_abline(aes(slope = 1, intercept = 0), linetype = 2) + 
#       labs(title = paste("pct.de", pct.de,", maturity.effect", maturity.effect))
#     count <- count + 1
#   }
# } 
# pl <- do.call("grid.arrange", c(plots, ncol=3))
# ggsave(file="figures/corr_error_de_pctde_mateffect.png", pl, width = 10, height = 9)

############### pct.de vs. corr (F-test / js) #########################
# corr.list <- c(0.2, 0.5)
# maturity.effect.list <- c(0.05, 0.2, 1, 5)
# plots <- list()
# count <- 1
# for (corr in corr.list) {
#   for (maturity.effect in maturity.effect.list) {
#     print(paste(corr, maturity.effect))
#     data <- gen_data(ngenes = 500, corr = corr, maturity.effect = maturity.effect, pct.corr = 0.3, pct.de = 0.2)
#     pval.ftest <- get_ftest_pval_dist(data$expr)
#     pval.js <- get_jackstraw_pval_dist(data$expr)
#     # only plot the p-values for non-de genes
#     df <- data.frame(ftest = pval.ftest[data$genes.nonde], js = pval.js[data$genes.nonde])
#     plots[[count]] <- ggplot(data = df) +
#       stat_qq(aes(sample = ftest), distribution = stats::qunif, colour = "red") + 
#       stat_qq(aes(sample = js), distribution = stats::qunif, colour = "green") + 
#       geom_abline(aes(slope = 1, intercept = 0), linetype = 2) + 
#       labs(title = paste("corr", corr, ", maturity.effect", maturity.effect))
#     count <- count + 1
#   }
# } 
# pl <- do.call("grid.arrange", c(plots, ncol=4))
# ggsave(file="figures/corr_error_de_corr_mateffect.png", pl, width = 14, height = 6)

############### pct.de vs. corr (F-test: pca / slingshot) #########################
# corr.list <- c(0.2, 0.5)
# maturity.effect.list <- c(0.05, 0.2, 1, 5)
# plots <- list()
# count <- 1
# for (corr in corr.list) {
#   for (maturity.effect in maturity.effect.list) {
#     print(paste(corr, maturity.effect))
#     data <- gen_data(ngenes = 500, corr = corr, maturity.effect = maturity.effect, pct.corr = 0.3, pct.de = 0.2)
#     pval1 <- get_ftest_pval_dist(data$expr, pt.method = "pca")
#     pval2 <- get_ftest_pval_dist(data$expr, pt.method = "slingshot")
#     # only plot the p-values for non-de genes
#     df <- data.frame(pval1 = pval1[data$genes.nonde], pval2 = pval2[data$genes.nonde])
#     plots[[count]] <- ggplot(data = df) +
#       stat_qq(aes(sample = pval1), distribution = stats::qunif, colour = "red") +
#       stat_qq(aes(sample = pval2), distribution = stats::qunif, colour = "green") +
#       geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
#       labs(title = paste("corr", corr, ", maturity.effect", maturity.effect))
#       count <- count + 1
#   }
# }
# pl <- do.call("grid.arrange", c(plots, ncol=4))
# ggsave(file="figures/corr_error_de_corr_mateffect_slingshot.png", pl, width = 14, height = 6)

