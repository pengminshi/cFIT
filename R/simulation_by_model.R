#' Simulate data from the cFIT model
#'
#' Simulate \deqn{X_j = H_jW\Lambda_j + 1_{n_j}b_j +E_j, j = 1,\ldots, ntask}. The nonnegative matrix
#' is generated as the cluster centers given number of clusters K. H_j is generated as the binary
#' membership matrix, where the labels are generated from a Dirichlet distribution with parameter
#' alpha. Distortion lambda and shift b are generated from truncated normal distribution, Noise
#' matrix E is generate with each entry from iid normal distribution.
#'
#' @param n number of data point per dataset
#' @param ntask number of batches
#' @param K number of clusters
#' @param p number of genes
#' @param alpha parameter for Dirichilet distribution used to generate the labels (default 10,
#' representing equal cluster size. smaller alpha corresponds to more unbalanced types)
#' @param sig with cluster variance
#' @param cl.sep Cluster center separation, the higher the clusters are well separated
#' @param batch.effect.sig batch effect variance, higher the large batch effects are
#'
#' @return A list of generated data, \describe{
#'  \item{X.list}{a list of n-by-p expression matrix}
#'  \item{H.list}{a list of n-by-K binary membership matrix}
#'  \item{lambda.list}{a list of length-p vectors of per-dataset scaling}
#'  \item{b.list}{a list of length-p vectors of per-dataset shift}
#'  \item{E.list}{a list of noise matrix}
#'  \item{label.list}{a list of length-n vectors of cluster labels}
#'  \item{W}{p-by-K common factor matrix}
#' }
#' @importFrom gtools rdirichlet
#' @export
generate_data <- function(n, ntask, K, p, alpha = NULL, sig = 1, cl.sep = 1, batch.effect.sig = 0.1) {

    if (is.null(alpha)) {
        alpha = rep(10, K)  # nearly 1/K for each cluster
    } else {
        if (length(alpha) == 1) {
            alpha = rep(alpha, K)
        }
    }

    # generate W: representing the cluster centerss
    centers.out = generate_centers_uniform(K = K, d = p, cl.sep = cl.sep)
    # W = t(abs(centers.out$centers)/rowSums(abs(centers.out$centers)) * 100) #
    # column sum is 100
    W = t(centers.out$centers)

    # generate labels
    label.list = lapply(1:ntask, function(l) {
        # generate proportions of gaussian mixture
        prob = gtools::rdirichlet(1, alpha = alpha)

        # generate labels
        sample(1:K, size = n, replace = T, prob = prob)
    })

    # generate Hj
    H.list = lapply(1:ntask, function(l) {

        # H, binary and row sum is 1
        label_to_membership(label.list[[l]], label.names = 1:K)  # n * K
    })

    # generate lambdaj, mean 1 sd = 0.5, trimed at 0
    lambda.list = lapply(1:ntask, function(l) {
        lambd = rnorm(p, mean = 1, sd = batch.effect.sig)
        lambd[lambd < 0] = 0
        return(lambd)
    })

    # generate bj, mean 0 abs(rnorm (0, sig))
    b.list = lapply(1:ntask, function(l) {
        abs(rnorm(p, mean = 0, sd = batch.effect.sig))
    })

    # generate Ej
    E.list = lapply(1:ntask, function(l) {
        matrix(rnorm(n * p, mean = 0, sd = sig), nrow = n)
    })

    # generate Ej to get Xj
    X.list = lapply(1:ntask, function(l) {
        X = H.list[[l]] %*% t(W * lambda.list[[l]]) + matrix(1, nrow = n, ncol = 1) %*%
            b.list[[l]] + E.list[[l]]
        X[X < 0] = 0

        colnames(X) = paste0("gene_", 1:ncol(X))
        rownames(X) = paste0("cell_", l, 1:nrow(X))
        X
    })

    return(list(X.list = X.list, H.list = H.list, lambda.list = lambda.list, b.list = b.list,
        E.list = E.list, label.list = label.list, W = W))

}

#' Generate gaussian centers with fixed minimum separation c
#'
#' Generate K centers in d dimenaional space such that min ||mu_i-mu_j|| = c(d * sig)^{1/2]}
#'
#' @param K number of clusters
#' @param d dimension of space where K centers are generated
#' @param cl.sep cluster separation
#'
#' @return A list containing \describe{
#'  \item{centers}{generated K-by-d corresponding to K centers in d-dim space}
#'  \item{status}{boolean indicating whether the generationg is successful}
#' }
#' @import checkmate
generate_centers_uniform <- function(K, d, cl.sep) {
    checkmate::assert_true(K >= 2)

    val_norm = sapply(runif(K * d), dist_norm, d, K)

    centers = matrix(val_norm * cl.sep, nrow = K)

    return(list(centers = centers, status = T))
}



dist_norm <- function(x, d, K) {
    x/gamma(1 + 1/d)/gamma(K) * gamma(1/d + K)
}




#' Convert label vector to membership matrix
#'
#' @param labels a vector of labels from K distint classes
#' @param label.names (optional) alternative label names, used for naming columns of membership matrix
#'
#' @return an n-by-K binary membership matrix
#' @import checkmate
#' @export
label_to_membership <- function(labels, label.names = NULL) {
    if (is.null(label.names)) {
        label.names = sort(unique(labels))
    } else {
        checkmate::assert_true(all(labels %in% label.names))
    }

    memb = t(sapply(labels, function(lab) as.numeric(label.names == lab)))

    return(memb)
}


