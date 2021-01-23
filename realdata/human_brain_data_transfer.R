
library(cFIT)
data.path = ''



################################################################
# load zhong count data matrix
zhong = readRDS(paste0(data.path, 'zhong/processed/zhong.rds'))
# str(zhong)

# load reference factors from intergration output of Geschwind and Kreigstein data
int.out =  readRDS(paste0(results.path, 'human_brain_int_out.rds'))
# str(int.out)


# find the intesection of gene s
genes = intersect(rownames(int.out$W), colnames(zhong$counts))
# str(genes)

Wref = int.out$W[genes,]; str(Wref) # reference factor matrix

exprs = zhong$counts[, genes]
exprs = log2(exprs/rowSums(exprs) * 10^4 + 1)
exprs = scale(exprs, center=F, scale=T)
# str(exprs)


# perform transfer
out = data_transfer(Xtarget=exprs, Wref=Wref, max.niter=100, verbose=T)
saveRDS(out, paste0(results.path, 'zhong_transfer_out.rds'))



################################################################
# camp data

# load camp data from SOUP package
camp = SOUP::camp

# load reference factors from intergration output of Geschwind and Kreigstein data
int.out =  readRDS(paste0(results.path, 'human_brain_int_out.rds'))
# str(int.out)

# find the intesection of gene s
genes = intersect(rownames(int.out$W), colnames(camp$counts))
str(genes)

Wref = int.out$W[genes,]; str(Wref) # reference factor matrix

exprs = camp$counts[, genes]
exprs = log2(exprs/rowSums(exprs) * 10^4 + 1)
exprs = scale(exprs, center=F, scale=T)
str(exprs)

out = data_transfer(Xtarget=exprs, Wref=Wref, max.niter=100, verbose=T)
saveRDS(out, file=paste0(results.path, 'camp_transfer_out.rds'))

################################################################

# darmanis data
int.out =  readRDS(paste0(results.path, 'human_brain_int_out.rds'))

# adult cells
dat = readRDS(paste0(data.path, 'darmanis/processed/darmanis_adult_332cells.rds'))

# find the intesection of gene s
genes = intersect(rownames(int.out$W), colnames(dat$counts))
# str(genes)

Wref = int.out$W[genes,]; str(Wref) # reference factor matrix

exprs = dat$counts[, genes]
exprs = log2(exprs/rowSums(exprs) * 10^4 + 1)
exprs = scale(exprs, center=F, scale=T)
str(exprs)


out = data_transfer(Xtarget=exprs, Wref=Wref, max.niter=100, seed=42, verbose=T)
saveRDS(out, paste0(results.path, 'darmanis_adult_transfer_out.rds'))


# fetal cells
dat = readRDS(paste0(data.path, 'darmanis/processed/darmanis_fetal_134cells.rds'))

# find the intesection of gene s
genes = intersect(rownames(int.out$W), colnames(dat$counts))
# str(genes)

Wref = int.out$W[genes,]; str(Wref) # reference factor matrix

exprs = dat$counts[, genes]
exprs = log2(exprs/rowSums(exprs) * 10^4 + 1)
exprs = scale(exprs, center=F, scale=T)
str(exprs)


out = data_transfer(Xtarget=exprs, Wref=Wref, max.niter=100, seed=42, verbose=T)
saveRDS(out, paste0(results.path, 'darmanis_fetal_transfer_out.rds'))


################################################################
# Sestan data

# load sestan count data matrix
dat = readRDS(paste0(data.path, 'Sestan/processed/sestan.rds'))
# load reference factors from intergration output of Geschwind and Kreigstein data
int.out =  readRDS(paste0(results.path, 'human_brain_int_out.rds'))

# find the intesection of gene
genes = intersect(rownames(int.out$W), colnames(dat$counts))
# str(genes)

Wref = int.out$W[genes,]; str(Wref) # reference factor matrix

exprs = dat$counts[, genes]
exprs = log2(exprs/rowSums(exprs) * 10^4 + 1)
exprs = scale(exprs, center=F, scale=T)
str(exprs)


out = data_transfer(Xtarget=exprs, Wref=Wref, max.niter=100, seed=42, verbose=T)
saveRDS(out, paste0(results.path, 'sestan_transfer_out.rds'))



