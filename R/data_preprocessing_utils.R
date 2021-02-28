#' Split the expression matrix into a list of per batch expressions
#'
#' @param X An ncells-by-ngenes gene expression matrix
#' @param batch categorical/characteristic vector as batch id
#' @param labels (optional) charateristic vector of cell labels (default NULL)
#' @param metadata (optinal) dataframe of ncells row saving the cell metadata
#' @param dataset.name string, prefix of list name
#'
#' @return a list containing  \describe{
#'  \item{X.list}{a list of gene expression matrix per batch}
#'  \item{labels.list}{a list of vectors, saving the cell labels per batch}
#'  \item{metadata.list}{a list of dataframe saving the metadata per batch}
#' }
#'
#' @import checkmate
#' @export
split_dataset_by_batch <- function(X, batch, labels = NULL, metadata = NULL, dataset.name = NULL) {

    checkmate::assert_true(nrow(X) == length(batch))
    batch = as.character(batch)

    batch.unique = unique(batch)
    nbatch = length(batch.unique)

    # initialize
    X.list = vector(mode = "list", length = nbatch)
    labels.list = vector(mode = "list", length = nbatch)
    metadata.list = vector(mode = "list", length = nbatch)

    # save the X, labels, metadata in batches
    for (i in 1:nbatch) {
        ind = which(batch == batch.unique[i])
        X.list[[i]] = X[ind, ]
        labels.list[[i]] = as.character(labels[ind])
        if (!is.null(rownames(X.list[[i]]))) {
            names(labels.list[[i]]) = rownames(X.list[[i]])
        }

        if (!is.null(metadata)) {
            metadata.list[[i]] = metadata[ind, ]
            if (!is.null(rownames(X.list[[i]]))) {
                rownames(metadata.list[[i]]) = rownames(X.list[[i]])
            }
        }
    }

    # name the batchs adding dataset.name
    names(X.list) = paste0(dataset.name, batch.unique)
    names(labels.list) = paste0(dataset.name, batch.unique)
    names(metadata.list) = paste0(dataset.name, batch.unique)

    return(list(X.list = X.list, labels.list = labels.list, metadata.list = metadata.list))
}




#' Select highly variabel genes shared among datasets
#'
#' wrapper for Seurat gene selection functions \code{FindVariableFeatures} and \code{SelectIntegrationFeatures}.
#'
#' @param X.list a list of m cells-by-genes gene expression matrices from m data sets
#' @param ngenes integer, number of highly variable genes to select
#' @param verbose boolean scalar, whether to show Seurat gene selection messages (default FALSE)
#'
#' @import Seurat
#' @export
select_genes <- function(X.list, ngenes = 3000, subset = NULL, verbose = F) {
    m = length(X.list)
    obj.list = lapply(1:m, function(j) {
        obj = Seurat::CreateSeuratObject(counts = t(X.list[[j]]))
        obj = Seurat::NormalizeData(obj, verbose = verbose)
        obj = Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = ngenes,
            verbose = verbose)
        obj
    })
    selected.genes = Seurat::SelectIntegrationFeatures(object.list = obj.list, nfeatures = ngenes,
        verbose = TRUE)
    return(selected.genes)
}



#' Preprocess expression matrix as input for data integration
#'
#' Preprocess, leave highly expressed genes, depth normalized, log transform, scale but not centered
#'
#' @param X.list A list of m ncells-by-ngenes, gene expression matrices from m data sets
#' @param genes character vector, a vector of gene names (selected highly variable genes)
#' @param scale.factor scalar, scaling factor for library size normalization (default 1e4)
#' @param scale boolean, whether to scale per gene expression (default TRUE)
#' @param center boolean, whether to center per gene expression (default FALSE)
#'
#' @return preprocessed expression list X.list
#'
#' @importFrom Seurat ScaleData
#' @export
preprocess_for_integration <- function(X.list, genes, scale.factor = 10000, scale = T,
    center = F, verbose = F) {
    m = length(X.list)
    datasets = names(X.list)
    for (i in 1:m) {
        genes = intersect(genes, colnames(X.list[[i]]))
    }

    X.list = lapply(1:m, function(j) {
        x = X.list[[j]][, genes]
        x = log(x/rowSums(x) * scale.factor + 1)
        x = t(Seurat::ScaleData(t(x), do.center = center, do.scale = scale, verbose = verbose))
        x
    })
    names(X.list) = datasets

    return(X.list)
}

#' Preprocess expression matrix as input for data transfer
#'
#' Subset to the shared gene sets between target data and reference factor matrix.
#' Depth normalized, log transform, scale but not centered
#'
#' @param counts ncells-by-ngenes raw count matrix from target data
#' @param Wref ngenes-by-r reference factor matrix
#' @param scale.factor scalar, scaling factor for library size normalization (default 1e4)
#' @param scale boolean, whether to scale per gene expression (default TRUE)
#' @param center boolean, whether to center per gene expression (default FALSE)
#' @param verbose boolean scalar, whether to show extensive program logs (default TRUE)
#'
#' @return A list containing \describe{
#'  \item{exprs}{ncell-by-ngenes preprocessed expresssion matrix}
#'  \item{Wref}{ngenes-by-r reference factor matrix, with genes filtered}
#' }
#'
#' @importFrom Seurat ScaleData
#' @export
preprocess_for_transfer <- function(counts, Wref, scale.factor = 10000, scale = T,
    center = F, verbose = F) {
    genes = intersect(rownames(Wref), colnames(counts))

    if (verbose)
        logmsg("Transfer on ", nrow(counts), ", ", length(genes), " genes ...")

    Wref = Wref[genes, ]
    exprs = counts[, genes]
    exprs = log(exprs/rowSums(exprs) * scale.factor + 1)
    exprs = t(Seurat::ScaleData(t(exprs), do.center = center, do.scale = scale, verbose = verbose))

    return(list(exprs = exprs, Wref = Wref))
}


#' Match the human gene expressions to the mouse gene expressions
#'
#' Match the expression of human genes to expression of mouse genes. If the case of multi-match, the expressions are averaged
#'
#' @param human.expres ncell-by-ngene expression matrix of human
#' @param mouse.genes.list a character vector of mouse gene symbols, indicating which genes are to match
#'
#' @return ncell-by-ngene' of mapped mouse expression matrix
#'
#' @export
match_human_expr_to_mouse_genes <- function(human.exprs, mouse.genes.list) {

    expr.human.genes = colnames(human.exprs)

    map.list = mouse_gene_to_human_genes(mouse.genes.list)  # find the mapped human gene of each mouse gene
    mapped.genes = names(map.list)
    nb.mapped.genes = length(mapped.genes)

    mapped.genes.exist = rep(T, nb.mapped.genes)
    gene.exprs = lapply(1:nb.mapped.genes, function(i) {
        gene = mapped.genes[i]

        human.genes = map.list[[gene]]
        human.genes = unique(human.genes[human.genes %in% expr.human.genes])

        nb.genes = length(human.genes)

        if (nb.genes >= 1) {
            # multiple match
            if (nb.genes > 1) {
                cat("mouse gene ", gene, "mapped to human genes", human.genes, "\n")
            }

            exprs.of.gene = human.exprs[, human.genes]
            if (nb.genes >= 2) {
                exprs.of.gene = rowMeans(exprs.of.gene)
            }
            return(exprs.of.gene)
        } else {
            mapped.genes.exist[i] <<- F
        }
        return(NULL)
    })
    gene.exprs = do.call(cbind, gene.exprs)
    colnames(gene.exprs) = mapped.genes[mapped.genes.exist]

    return(gene.exprs)
}

#' Find the mapped human gene of each mouse gene
#'
#' @param gene.list a character vector of mouse gene symbols, indicating which mouse genes are to match
#' @return a list of vectors, corresponding to the mapped human genes for each mouse gene
#' @export
mouse_gene_to_human_genes <- function(gene.list) {

    if (!require("biomaRt", quietly = TRUE)) {
        message("install biomaRt package ...")
        BiocManager::install("biomaRt")
        require(biomaRt, quietly = TRUE)
    }

    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

    gene.map = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = gene.list,
        mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows = T)

    tmp = aggregate(as.character(gene.map[, 2]), by = list(gene.map[, 1]), FUN = function(x) list(as.character(x)))

    map.list = tmp$x
    names(map.list) = tmp$Group.1
    return(map.list)
}

#' Match the mouse gene expressions to the human gene expressions
#'
#' Match the expression of mouse genes to expression of human genes.
#' If the case of multi-match, the expressions are averaged
#'
#' @param mouse.expres ncell-by-ngene expression matrix of mouse
#' @param human.genes.list a character vector of human gene symbols, indicating which genes are to match
#'
#' @return ncell-by-ngene' of mapped human expression matrix
#'
#' @export
match_mouse_expr_to_human_genes <- function(mouse.exprs, human.genes.list) {
    expr.mouse.genes = colnames(mouse.exprs)

    map.list = human_gene_to_mouse_genes(human.genes.list)
    mapped.genes = names(map.list)
    nb.mapped.genes = length(mapped.genes)

    mapped.genes.exist = rep(T, nb.mapped.genes)
    gene.exprs = lapply(1:nb.mapped.genes, function(i) {
        gene = mapped.genes[i]

        mouse.genes = map.list[[gene]]
        mouse.genes = unique(mouse.genes[mouse.genes %in% expr.mouse.genes])

        nb.genes = length(mouse.genes)

        if (nb.genes >= 1) {
            exprs.of.gene = mouse.exprs[, mouse.genes]
            if (nb.genes >= 2) {
                exprs.of.gene = rowMeans(exprs.of.gene)
            }
            return(exprs.of.gene)
        } else {
            mapped.genes.exist[i] <<- F
        }
        return(NULL)
    })
    gene.exprs = do.call(cbind, gene.exprs)
    colnames(gene.exprs) = mapped.genes[mapped.genes.exist]
    return(gene.exprs)
}

#' Find the mapped mouse gene of each human gene
#'
#' @param gene.list a character vector of human gene symbols, indicating which human genes are to match
#' @return a list of vectors, corresponding to the mapped mouse genes for each human gene
#' @export
human_gene_to_mouse_genes <- function(gene.list) {

    if (!require("biomaRt", quietly = TRUE)) {
        message("install biomaRt package ...")
        BiocManager::install("biomaRt")
        require("biomaRt")
    }

    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

    gene.map = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = gene.list,
        mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows = T)

    tmp = aggregate(as.character(gene.map[, 2]), by = list(gene.map[, 1]), FUN = function(x) list(as.character(x)))

    map.list = tmp$x
    names(map.list) = tmp$Group.1
    return(map.list)
}

