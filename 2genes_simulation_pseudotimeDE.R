# simulate 2 genes from negative-binomial
# compare F-test and pseudotimeDE p-values
# this simulation is NOT fleshed out

library(dplyr)
library(ggplot2)
library(PseudotimeDE)
library(SingleCellExperiment)
library(tibble)

ncells <- 500
ngenes <- 2 
nreps <- 500 # num replications for null distribution of p-values
nsubsamples <- 100 # number of subsamples for pseudotimeDE
  
pval <- matrix(NA, nrow=nreps, ncol=2) # the empirical p-values distribution 
colnames(pval) <- c("F-test", "PseudotimeDE")
for (i in 1:nreps) {
  if (i %% 10 == 0) print(i)

  # generate data from NB(2, 0.5)
  gene.exp <- matrix(rnbinom(ncells * ngenes, size = 1, prob = 0.5), ncol = ngenes)
  # gene.exp <- matrix(rnorm(ncells * ngenes), ncol=ngenes)
  # pseudotime is PC 1
  pt <- as.vector(svd(gene.exp, nu=1, nv=1)$u)
  df <- data.frame(gene1 = gene.exp[, 1], gene2 = gene.exp[, 2], pt = pt)
  if (i == 1) {
    # example data scatter plot
    pl.scatter <- ggplot(df, aes(x = gene1, y = gene2, color = pt)) + 
      geom_point()
  }

  # 1. f-test using linear regression 
  mod <- lm(gene1 ~ pt, data=df)
  pval[i, 1] <- summary(mod)$coefficients[2, 4]
  
  # 2. pseudotimeDE test
  sce <- SingleCellExperiment(list(counts=t(gene.exp)))
  colnames(sce) <- c(1:ncells)
  rownames(sce) <- c(1:ngenes)
  ori.tbl <- tibble(cell = colnames(sce), pseudotime = pt)

  options(mc.cores = 8)
  # get cell indices for each subsample
  cell.index <- mclapply(seq_len(nsubsamples), function(x) {
    sample(x = c(1:dim(sce)[2]), size = 0.8*dim(sce)[2], replace = FALSE)
  })

  sub.tbl <- mclapply(cell.index, function(x, sce) {
    # repeat pseudotime estimation for each subsample
    sce.sub <- sce[, x]
    pt.sub <- as.vector(svd(t(assays(sce.sub)$counts), nu=1, nv=1)$u)
    tbl <- tibble(cell = colnames(sce.sub), pseudotime = pt.sub)
    tbl
  }, sce = sce)

  res <- PseudotimeDE::runPseudotimeDE(gene.vec = c("1"),
                                       ori.tbl = ori.tbl,
                                       sub.tbl = sub.tbl,
                                       sce = sce,
                                       model = "nb")
  pval[i, 2] <- res$para.pv
}

pl.pval.ftest <- ggplot(data.frame(pval = pval[, 1]), aes(sample = pval)) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5, color='red') +
  stat_qq(distribution = stats::qunif) + 
  labs(title = "Q-Q plot: F-test")

pl.pval.pde <- ggplot(data.frame(pval = pval[, 2]), aes(sample = pval)) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5, color='red') +
  stat_qq(distribution = stats::qunif) + 
  labs(title = "Q-Q plot: PseudotimeDE")
