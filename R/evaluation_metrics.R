#' Calculate the alginment score
#'
#' Calculate the alginment score per sample, per cluster and per dataset
#'
#' @param X.list a list of gene expression matrix per batch
#' @param k integer, number of nearest neighbor to consider
#' @param balanced boolean, whether to subset to balance the size of baches
#' @param min.sample.per.dataset integer, minimum cells to sample per dataset (default 100)
#' @param dataset factor/characteristic vector that indicate the batch of cells if len(X.list)==1
#' @param labels factor/characteristic vector, cluster/cell type labels of the cells
#' @param max.k integer, maximum number of nearest neighbors to consider (cap for large dataset)
#' @param seed random seed
#' @param verbose boolean scalar, whether to show extensive program logs (default TRUE)
#' @param pca if pca!=F, perform pca first before calculating the pairwise distance
#'
#' @return A list containing \describe{
#'  \item{alignment.per.sample}{a numeric vector, alignment score per cell}
#'  \item{alignment.per.cluster}{a numeric vector, alignment score per cluster. Only calculated}
#'  with labels are provided
#'  \item{alignment.per.dataset}{a numeric vector, alignment score per dataset.}
#' }
#'
#' @importFrom FNN get.knn
#' @export
calculate_alignment <- function(X.list, k = NULL, balanced = F, min.sample.per.dataset = 100, 
    dataset = NULL, labels = NULL, max.k = 50, seed = 0, verbose = F, pca = F) {
    set.seed(seed)
    
    m = length(X.list)
    
    if (m == 1 & is.null(dataset)) {
        stop("Error: please provide more than one dataset!")
    } else if (m == 1) {
        unique.ds = unique(dataset)
        X.list = lapply(unique.ds, function(d) X.list[[1]][dataset == d, ])
        m = length(X.list)
    }
    
    if (verbose) 
        logmsg("Calculate alignment of ", m, " datasets ...")
    
    if (is.null(dataset)) {
        dataset = as.character(do.call(c, lapply(1:m, function(j) rep(j, nrow(X.list[[j]])))))
    } else {
        dataset = as.character(dataset)
    }
    
    if (balanced) {
        min.sample.per.dataset = max(min(sapply(X.list, nrow)), min.sample.per.dataset)
        if (verbose) 
            logmsg("Subsample tp ", min.sample.per.dataset, " points per dataset ...")
        
        tmp = 0
        ind.list = c()
        for (j in 1:m) {
            nj = nrow(X.list[[j]])
            ind = sample(1:nj, min.sample.per.dataset, replace = F)
            X.list[[j]] = X.list[[j]][ind, ]
            ind.list = c(ind.list, tmp + ind)
            
            tmp = tmp + nj
        }
        dataset = dataset[ind.list]
        if (!is.null(labels)) {
            labels = labels[ind.list]
        }
    }
    
    n = sum(sapply(X.list, nrow))
    X = do.call(rbind, X.list)
    if (pca) {
        X = prcomp(X, rank = max(pca, 10))$x
    }
    
    if (is.null(k)) {
        k = min(max.k, max(10, floor(0.01 * n)))
    }
    
    if (verbose) 
        logmsg("Find k nearest neighbor with k = ", k, "...")
    knn.out = FNN::get.knn(X, k = k)
    
    
    unique.ds = unique(dataset)
    dataset.ratio = sapply(unique.ds, function(ds) mean(dataset == ds))  # replace 1/m
    names(dataset.ratio) = unique.ds
    
    
    # n = nrow(knn.out$nn.index)
    if (verbose) 
        logmsg("Calculate alignment score")
    alignment.per.sample = sapply(1:n, function(i) {
        nb_same = sum(dataset[knn.out$nn.index[i, ]] == dataset[i])
        ratio = dataset.ratio[dataset[i]]
        1 - (nb_same/k - ratio)/(1 - ratio)  # normalization
    })
    
    if (!is.null(labels)) {
        if (verbose) 
            logmsg("Aggregate alignment score per cluster ...")
        agg.out = aggregate(alignment.per.sample, by = list(labels = labels), FUN = mean)
        alignment.per.cluster = agg.out$x
        names(alignment.per.cluster) = agg.out$labels
    } else {
        alignment.per.cluster = NULL
    }
    
    if (verbose) 
        logmsg("Aggregate alignment score per dataset ...")
    agg.out = aggregate(alignment.per.sample, by = list(dataset = dataset), FUN = mean)
    alignment.per.dataset = agg.out$x
    names(alignment.per.dataset) = agg.out$dataset
    
    return(list(alignment.per.sample = alignment.per.sample, alignment.per.cluster = alignment.per.cluster, 
        alignment.per.dataset = alignment.per.dataset))
}


#' Calculate the alginment score
#'
#' Calculate the alginment score per cluster, where the normalization is performed based on the
#' per-group batch population ratios
#'
#' @param X ncell-by-ngene expression matrix
#' @param dataset factor/characteristic vector that indicate the batch of cells
#' @param labels factor/characteristic vector, cluster/cell type labels of the cells
#' @param k integer, number of nearest neighbor to consider
#' @param max.k integer, maximum number of nearest neighbors to consider (cap for large dataset)
#' @param seed random seed
#' @param verbose boolean scalar, whether to show extensive program logs (default TRUE)
#' @param pca if pca!=F, perform pca first before calculating the pairwise distance
#'
#' @return A list containing \describe{
#'  \item{alignment.per.sample}{a numeric vector, alignment score per cell}
#'  \item{alignment.per.cluster}{a numeric vector, alignment score per cluster. Only calculated
#'  with labels are provided}
#'  \item{alignment.per.dataset}{a numeric vector, alignment score per dataset.}
#'  \item{alignment.per.dataset.cluster}{numeric vector, alignment score per dataset per cluster}
#' }
#'
#' @importFrom FNN get.knn
#' @export
calculate_alignment_per_type <- function(X, dataset = NULL, labels = NULL, k = NULL, 
    max.k = 50, seed = 0, verbose = F, pca = F) {
    
    set.seed(seed)
    n = nrow(X)
    
    if (pca) {
        X = prcomp(X, rank = max(pca, 10))$x
    }
    
    if (is.null(k)) {
        k = min(max.k, max(10, floor(0.01 * n)))
    }
    
    if (verbose) 
        logmsg("Find k nearest neighbor with k = ", k, "...")
    knn.out = FNN::get.knn(X, k = k)
    
    # n = nrow(knn.out$nn.index)
    if (verbose) 
        logmsg("Calculate alignment score")
    alignment.per.sample = sapply(1:n, function(i) {
        nb.neighbor.in.same.label = sum(labels[knn.out$nn.index[i, ]] == labels[i])
        if (nb.neighbor.in.same.label > 0) {
            nb_same = sum(labels[knn.out$nn.index[i, ]] == labels[i] & dataset[knn.out$nn.index[i, 
                ]] == dataset[i])
            ratio = mean(dataset == dataset[i] & labels == labels[i])/mean(labels == 
                labels[i])  # proportion of dataset conditioned on label
            return(1 - (nb_same/nb.neighbor.in.same.label - ratio)/(1 - ratio))  # # normalization
        } else {
            return(NA)
        }
    })
    
    if (verbose) 
        logmsg("Aggregate alignment score per cluster ...")
    agg.out = aggregate(alignment.per.sample, by = list(labels = labels), FUN = mean, 
        na.rm = T)
    alignment.per.cluster = agg.out$x
    names(alignment.per.cluster) = agg.out$labels
    
    if (verbose) 
        logmsg("Aggregate alignment score per dataset ...")
    agg.out = aggregate(alignment.per.sample, by = list(dataset = dataset), FUN = mean, 
        na.rm = T)
    alignment.per.dataset = agg.out$x
    names(alignment.per.dataset) = agg.out$dataset
    
    if (verbose) 
        logmsg("Aggregate alignment per cluster per dataset")
    alignment.per.dataset.cluster = aggregate(alignment.per.sample, by = list(dataset = dataset, 
        labels = labels), FUN = mean, na.rm = T)
    
    
    return(list(alignment.per.sample = alignment.per.sample, alignment.per.cluster = alignment.per.cluster, 
        alignment.per.dataset = alignment.per.dataset, alignment.per.dataset.cluster = alignment.per.dataset.cluster))
}


#' Map cell labels for each cell type for target data given the cell type labels of reference data
#'
#' Assign labels to each cell by the mode of neighboring reference cell types. Calculate the proportion
#' of mapped reference cell types for each target cell types.
#'
#' @param exprs.source ncells-by-ngenes expression matrix of source data
#' @param exprs.target ncells-by-ngenes expression matrix of target data
#' @param labels.source a vector of labels for source cells
#' @param labels.source a vector of labels for target cells
#' @param k integer, number of nearest neighor (by default \code{min(max.k,
#' max(10, floor(0.01* (n.source + n.target))))} is used)
#' @param max.k upper bound of k (default 50)
#' @param verbose boolean scalar, whether to show extensive program logs (default TRUE)
#'
#' @importFrom FNN get.knn
#' @export
get_cluster_match <- function(exprs.source, exprs.target, labels.source, labels.target, 
    k = NULL, max.k = 50, verbose = F) {
    
    X.list = list(exprs.source, exprs.target)
    n.source = nrow(exprs.source)
    n.target = nrow(exprs.target)
    
    dataset = c(rep("source", n.source), rep("target", n.target))
    
    if (is.null(k)) {
        k = min(max.k, max(10, floor(0.01 * (n.source + n.target))))
    }
    
    if (verbose) 
        logmsg("Find k nearest neighbor with k = ", k, "...")
    knn.out = FNN::get.knn(do.call(rbind, X.list), k = k)
    
    target.knn.index = knn.out$nn.index[-(1:n.source), ]
    target.knn.index[target.knn.index > n.source] = NA  # only consider index from source
    
    if (verbose) 
        logmsg("Calculate major neighbor type per target sample...")
    matched.type = unlist(lapply(1:nrow(target.knn.index), function(i) {
        ind = target.knn.index[i, ]
        if (all(is.na(ind))) 
            return(NA)
        return(getmode(labels.source[ind[!is.na(ind)]]))  # matched major source label for each target sample
    }))
    
    
    if (verbose) 
        logmsg("Calculate the neighbor composition of each class...")
    unique.labels.target = unique(labels.target)
    source.class.size = table(labels.source)
    matched.type.per.label = lapply(unique.labels.target, function(label) {
        ind = which(labels.target == label)
        tab = table(matched.type[ind])
        
        # tab = tab / source.class.size[names(tab)] # rewight by the inverse of class
        # size
        tab = sort(tab, decreasing = T)
        
        if (length(tab) > 0) {
            tab = tab/length(ind)  #sum(tab) # calculate the proportion
        }
        list(n = length(ind), matched.type = tab)
    })
    names(matched.type.per.label) = unique.labels.target
    
    # return(list(target.knn.index=target.knn.index, matched.type=matched.type,
    # unique.labels.target=unique.labels.target, labels.target=labels.target))
    return(matched.type.per.label)
}

#' Assign labels for each cell in for target data given the reference labels
#'
#' Assign labels to each cell by the mode of neighboring reference cell types
#'
#' @param exprs.source ncells-by-ngenes expression matrix of source data
#' @param exprs.target ncells-by-ngenes expression matrix of target data
#' @param labels.source a vector of labels for  source cells
#' @param k integer, number of nearest neighor (by default \code{min(max.k,
#' max(10, floor(0.01* (n.source + n.target))))} is used)
#' @param max.k upper bound of k (default 50)
#' @param verbose boolean scalar, whether to show extensive program logs (default TRUE)
#'
#' @importFrom FNN get.knn
#' @export
asign_labels <- function(exprs.source, exprs.target, labels.source, k = NULL, max.k = 50, 
    verbose = F) {
    
    X.list = list(exprs.source, exprs.target)
    n.source = nrow(exprs.source)
    n.target = nrow(exprs.target)
    
    dataset = c(rep("source", n.source), rep("target", n.target))
    
    if (is.null(k)) {
        k = min(max.k, max(10, floor(0.01 * (n.source + n.target))))
    }
    
    if (verbose) 
        logmsg("Find k nearest neighbor with k = ", k, "...")
    knn.out = FNN::get.knn(do.call(rbind, X.list), k = k)
    
    target.knn.index = knn.out$nn.index[-(1:n.source), ]
    target.knn.index[target.knn.index > n.source] = NA  # only consider index from source
    
    if (verbose) 
        logmsg("Calculate major neighbor type per target sample...")
    matched.type = unlist(lapply(1:nrow(target.knn.index), function(i) {
        ind = target.knn.index[i, ]
        if (all(is.na(ind))) 
            return(NA)
        return(getmode(labels.source[ind[!is.na(ind)]]))  # matched major source label for each target sample
    }))
    return(matched.type)
}



#' Adjusted Rand Index
#'
#' @param labels.ref Reference label vector
#' @param labels.est Estimated label vector
#'
#' @return calculated ARI
#'
#' @importFrom  mclust adjustedRandIndex
#' @export
ari <- function(labels.ref, labels.est) {
    mclust::adjustedRandIndex(labels.ref, labels.est)
}

#' Get the the mode in a categorical vector
#'
#' @param x a categorical vector
#'
#' @return mode in the vector
#' @export
getmode <- function(x) {
    tab = table(x)
    mode = names(tab)[order(tab, decreasing = T)[1]]
    return(mode)
}


#' Plot the confution matrix
#'
#' @param est_label vector of estimated labels
#' @param true_label vector of true labels
#' @param short.names (optional) rename the reference true labels for plots
#' @param y.ord (optional) a vector lenght of the number of rows of table, re-ordering the rows
#' @param x.ord (optional) a vector lenght of the number of columns of table, re-ordering the columns
#' @param xlab (optional) string, x axis label (default reference)
#' @param ylab (optional) string, y axis label (default empty)
#' @param threshold (optional) non-negative value, lower threshold, above which the values are shown (default NULL)
#'
#' @return the plot (ggplot object)
#'
#' @import ggplot2
#' @export
plotContTable <- function(est_label, true_label, short.names = NULL, y.ord = NULL, 
    x.ord = NULL, xlab = "Reference", ylab = "", threshold = NULL) {
    
    if (is.null(short.names)) {
        short.names = levels(factor(true_label))
    }
    
    cont.table <- table(true_label, est_label)
    if (!is.null(y.ord)) {
        if (length(y.ord) != ncol(cont.table)) {
            cat("wrong order", ncol(cont.table), "columns")
        } else {
            cont.table = cont.table[, y.ord]
        }
    }
    if (!is.null(x.ord)) {
        if (length(x.ord) != nrow(cont.table)) {
            cat("wrong order", nrow(cont.table), "rows")
        } else {
            cont.table = cont.table[x.ord, ]
            short.names = short.names[x.ord]
        }
    }
    K <- ncol(cont.table)
    
    if (!is.null(threshold)) {
        cont.table[cont.table <= threshold] = 0
    }
    
    sub.clusters <- paste0("", colnames(cont.table))
    
    cont.table <- apply(as.matrix(cont.table), 2, as.integer)
    cont.table <- data.frame(cont.table)
    cont.table$Reference = factor(short.names, levels = short.names)
    colnames(cont.table) <- c(sub.clusters, "Reference")
    
    dat3 <- reshape2::melt(cont.table, id.var = "Reference")
    grid.labels = as.character(dat3$value)
    grid.labels[grid.labels == "0"] = ""
    
    
    g <- ggplot2::ggplot(dat3, aes(Reference, variable)) + geom_tile(aes(fill = value)) + 
        geom_text(aes(label = grid.labels), size = 4.5) + scale_fill_gradient(low = "white", 
        high = "lightsteelblue") + labs(y = ylab, x = xlab) + theme(panel.background = element_blank(), 
        axis.line = element_blank(), axis.text.x = element_text(size = 10, angle = 90, 
            hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 10), axis.ticks = element_blank(), 
        axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12), 
        legend.position = "none")
    
    return(g)
}
