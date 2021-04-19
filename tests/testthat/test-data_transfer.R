test_that("data transfer", {
    library(SeuratData)
    # InstallData("panc8") # if not installed, require memory of 118MB
    data("panc8")

    dat = panc8[, panc8$tech %in% c('fluidigmc1')]
    int.out = readRDS('../../vignettes/results/panc_ref_integration_output.rds')

    genes = intersect(rownames(int.out$W),rownames(dat))
    exprs.list = expect_error(preprocess_for_integration(list(X=t(as.matrix(dat@assays$RNA@counts))), genes,
                                                         scale.factor=10^4, scale=T, center=F), NA)

    # integration, run 5 iterations for testing
    transfer.out = expect_error(CFITTransfer(exprs.list[[1]], Wref=int.out$W, max.niter = 2,
                                             future.plan='sequential', verbose=F), NA)
    transfer.out = expect_error(CFITTransfer(exprs.list[[1]], Wref=int.out$W, max.niter = 2,
                                             future.plan='multicore', workers=4, verbose=F), NA)
    transfer.out = expect_error(CFITTransfer(exprs.list[[1]], Wref=int.out$W, max.niter = 2,
                                             future.plan='multisession', workers=4, verbose=F), NA)
})

