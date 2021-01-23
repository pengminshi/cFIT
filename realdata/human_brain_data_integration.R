library(cFIT)
devtools::load_all()
data.path = '../../gene/project/sc-integration-and-transfer/data/'


# geschwind
geschwind = readRDS(file = paste0(data.path, 'Geschwind/processed/geschwind.rds'))
genes = geschwind$select.genes[[1]] # 4074 genes
geschwind.split.out = split_dataset_by_batch(X=geschwind$counts[,colnames(geschwind$counts)%in% genes],
                                             labels=geschwind$cell.info$Cluster,
                                             batch = paste0(geschwind$cell.info$Library, geschwind$cell.info$Donor),
                                             dataset.name='Geschwind:')


# kreigstein
kreigstein = readRDS(file = paste0(data.path, 'Kreigstein/processed/kreigstein.rds'))
kreigstein.split.out = split_dataset_by_batch(X=kreigstein$counts[,colnames(kreigstein$counts) %in% genes],
                                              labels = kreigstein$cell.info$WGCNAcluster,
                                              batch = rep('',length(kreigstein$cell.info$WGCNAcluster)),
                                              dataset.name='Kreigstein:')


X.list = c( geschwind.split.out$X.list,  kreigstein.split.out$X.list)
genes = readRDS('results/geschwind_selected_genes.rds')
X.list = preprocess_for_integration(X.list=X.list, genes=genes)


set.seed(42)
r = 40
max.niter = 100
out = data_integrate(X.list=X.list, r=r, max.niter=max.niter, tol=1e-5)
saveRDS(out, paste0(results.path, 'human_brain_int_out.rds'))
