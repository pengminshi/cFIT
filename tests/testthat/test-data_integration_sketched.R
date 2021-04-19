

test_that("data integration sketched version", {
    library(SeuratData)
    # InstallData("panc8") # if not installed, require memory of 118MB
    data("panc8")

    dat = panc8[, panc8$tech %in% c('indrop','smartseq2')] # two largest
    data.list = expect_error(split_dataset_by_batch(X=t(as.matrix(dat@assays$RNA@counts)),
                                                    batch = dat$tech), NA)
    data.list = expect_error(split_dataset_by_batch(X=t(as.matrix(dat@assays$RNA@counts)),
                                                    batch = dat$tech,
                                                    labels = dat$celltype,
                                                    metadata = dat@meta.data,
                                                    dataset.name = 'panc:'), NA)
    # gene selection
    genes = expect_error(select_genes(data.list$X.list, ngenes=2000, verbose=F), NA)
    expect_equal(length(genes), 2000)

    # preprocess for integration
    exprs.list = expect_error(preprocess_for_integration(data.list$X.list, genes,
                                                         scale.factor=10^4, scale=T, center=F), NA)
    for(i in 1:length(exprs.list)){
        expect_true(all(exprs.list[[i]] >=0))
    }

    # integration, run 5 iterations for testing
    int.out = expect_error(CFITIntegrate_sketched(X.list=exprs.list, r=10, subsample.prop=0.2, max.niter = 2,
                                                  future.plan='sequential', verbose=F), NA)

    int.out = expect_error(CFITIntegrate_sketched(X.list=exprs.list, r=10, subsample.prop=0.2, max.niter = 2,
                                                  future.plan='multicore', workers=4, verbose=F), NA)

    int.out = expect_error(CFITIntegrate_sketched(X.list=exprs.list, r=10, subsample.prop=0.2, max.niter = 2,
                                                  future.plan='multisession', workers=4, verbose=F), NA) # very slow

})

test_that("data integration (sketched version) with large number of genes", {
    library(SeuratData)
    # InstallData("panc8") # if not installed, require memory of 118MB
    data("panc8")

    dat = panc8[, panc8$tech %in% c('indrop','smartseq2')] # two largest
    data.list = expect_error(split_dataset_by_batch(X=t(as.matrix(dat@assays$RNA@counts)),
                                                    batch = dat@meta.data$tech), NA)
    data.list = expect_error(split_dataset_by_batch(X=t(as.matrix(dat@assays$RNA@counts)),
                                                    batch = dat@meta.data$tech,
                                                    labels = dat@meta.data$celltype,
                                                    metadata = dat@meta.data,
                                                    dataset.name = 'panc:'), NA)
    # gene selection
    genes = select_genes(data.list$X.list, ngenes=8000, verbose=F)
    expect_equal(length(genes), 8000)

    # preprocess for integration
    exprs.list = expect_error(preprocess_for_integration(data.list$X.list, genes,
                                                         scale.factor=10^4, scale=T, center=F), NA)
    for(i in 1:length(exprs.list)){
        expect_true(all(exprs.list[[i]] >=0))
    }

    # integration, run 5 iterations for testing
    int.out = expect_error(CFITIntegrate_sketched(X.list=exprs.list, r=10, subsample.prop=0.2, max.niter = 2,
                                                  future.plan='sequential', verbose=F), NA)

    int.out = expect_error(CFITIntegrate_sketched(X.list=exprs.list, r=10, subsample.prop=0.2, max.niter = 2,
                                                  future.plan='multicore', workers=4, verbose=F), NA)

    int.out = expect_error(CFITIntegrate_sketched(X.list=exprs.list, r=10, subsample.prop=0.2, max.niter = 2,
                                                  future.plan='multisession', workers=4, verbose=F), NA) # very slow
})





