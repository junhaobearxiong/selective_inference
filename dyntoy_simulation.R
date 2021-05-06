# simulation to (roughly) replicate the simulation results in pseudotimeDE paper
# expression data simulation is based on `dyntoy`
# pseudotime is estimated using PC1 rather than slingshot
# this simulation is NOT fleshed out

library(dplyr)
library(ggplot2)
library(PseudotimeDE)
library(SingleCellExperiment)
library(tibble)
library(dyntoy)

nsubsamples <- 20 # number of subsamples for pseudotimeDE
ngenes.list <- c(100, 500) # number of genes
diffexp.list <- c(0, 0.2, 0.5) # proportion of DE genes

for (ngenes in ngenes.list) {
  for (diffexp in diffexp.list) {
    print(paste(ngenes, diffexp))
    
    # generate dataset using dyntoy
    dataset <- generate_dataset(model = "linear", num_cells = 100, num_features = ngenes, differentially_expressed_rate = diffexp)
    # convert to sce
    sce <- SingleCellExperiment(list(counts=t(dataset$counts)))
    
    # estimate pseudotime on the count data using pc1
    pt <- as.vector(svd(dataset$counts, nu=1, nv=1)$u)
    ori.tbl <- tibble(cell = colnames(sce), pseudotime = pt)
    
    # get de information of genes
    de <- dataset$tde_overall
    nonde.genes <- (de %>% filter(differentially_expressed == FALSE))$feature_id
    
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
    
    print(system.time(res <- PseudotimeDE::runPseudotimeDE(gene.vec = nonde.genes,
                                                           ori.tbl = ori.tbl,
                                                           sub.tbl = sub.tbl,
                                                           sce = sce,
                                                           model = "nb")
    ))
    
    pval <- res$para.pv
    saveRDS(pval, file=paste0('outputs/dyntoy_sim_ngene', ngenes, '_diffexp', diffexp, '.rds'))
  }
}

#####################################

for (ngenes in ngenes.list) {
  for (diffexp in diffexp.list) {
    pval <- readRDS(paste0('outputs/dyntoy_sim_ngene', ngenes, '_diffexp', diffexp, '.rds'))
    pl <- ggplot(data.frame(pval = pval), aes(sample = pval)) +
      geom_abline(intercept = 0, slope = 1, alpha = 0.5, color='red') +
      stat_qq(distribution = stats::qunif) +
      labs(title = paste0("QQ plot: # genes=", ngenes, ", DE rate=", diffexp))
    ggsave(paste0('figures/dyntoy_sim_ngene', ngenes, '_diffexp', diffexp, '.png'))
  }
}

