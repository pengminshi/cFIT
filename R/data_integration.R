#' Integration of multiple data source.
#'
#' Solve the model parameters through Iterative Nonnegative Matrix Factorization (iNMF),
#' by minimizing the objective function \deqn{1/N \sum_j||X_j -(H_JW^T\Lambda_j  + 1_nj b_j^T)||_F^2 +
#' gamma \sum_{l=1}^p(\sum_{j=1}^mnj/N \lambda_{jl}-1)^2}.
#'
#' @param X.list a list of m ncells-by-ngenes, gene expression matrices from m data sets
#' @param r scalar, dimension of common factor matrix, which can be chosen as the rough number of
#' identifiable cells types in the joint population (default 15).
#' @param max.niter integer, max number of iterations (default 100).
#' @param tol numeric scalar, tolerance used in stopping criteria (default 1e-5).
#' @param gamma numeric scalar, parameter for the penalty term. It is sufficient to chose a large
#' value (default 1e6)
#' @param nrep integer, number of repeated runs (to reduce effect of local optimum, default 1)
#' @param init a list of parameters for parameter initialization. The list either contains all
#' parameter sets: W,lambda.list, b.list, H.list, or only W will be used if provided (default NULL).
#' @param verbose boolean scalar, whether to show extensive program logs (default TRUE)
#' @param seed random seed used (default 0)
#' @param n.cores integer, number of cores used for parallel computation
#'
#' @return a list containing  \describe{
#'  \item{W}{ngenes-by-r numeric matrix, estimated common factor matrix}
#'  \item{H.list}{A list of m factor loading matrix of size ncells-by-r, estimated factor loading matrices}
#'  \item{b.list}{A list of estimated shift vector of size p (ngenes).}
#'  \item{lambda.list}{A list of estimated scaling vector of size p (ngenes).}
#'  \item{convergence}{boolean, whether the algorithm converge}
#'  \item{obj}{numeric scalar, value of the objective function at convergence or when maximum iteration achieved}
#'  \item{niter}{integer, the iteration at convergence (or maximum iteration if not converge)}
#'  \item{delta}{numeric scalar, the relative difference in last objective update}
#'  \item{params}{list of parameters used for the algorithm: gamma, max.iter, tol, nrep}
#' }
#'
#' @import checkmate parallel
#' @export
CFITIntegrate <- function(X.list, r = 15, max.niter = 100, tol = 1e-05, gamma = 1e+06,
    nrep = 1, init = NULL, verbose = T, seed = 0, n.cores = parallel::detectCores() -
        4) {

    message("Run cFIT with ", n.cores, " cores ...")
    m = length(X.list)

    # subset to the shared genes
    genes = colnames(X.list[[1]])
    for (j in 1:m) {
        if (is.null(colnames(X.list[[j]]))) {
            stop("gene symbols missing for import X.list")
        }
        genes = genes[genes %in% colnames(X.list[[j]])]
    }

    if (length(genes) < max(10, m)) {
        warning("Too few genes (", length(genes), "), check the data source")
    }

    p = length(genes)
    X.list = lapply(X.list, function(x) x[, genes])

    if (verbose)
        logmsg("Integrate ", m, " datasets (", p, " genes)")

    # total number of samples
    ntotal = sum(sapply(1:m, function(j) nrow(X.list[[j]])))

    # save the best results with minimum objective for output
    obj.best = Inf
    obj.history.best = NULL
    params.list.best = NULL
    niter.best = NULL
    delta.best = NULL
    delta.history.best = NULL
    deltaw.history.best = NULL  # diff w /w
    converge.best = NULL
    seed.best = NULL

    # run each repeat
    time.start = Sys.time()
    for (rep in 1:nrep) {
        set.seed(seed + rep - 1)  # set random seed

        # parameters known from previous iterations
        if (all(c("W", "lambda.list", "b.list", "H.list") %in% names(init))) {
            params.list = list(W = init$W, H.list = init$H.list, b.list = init$b.list,
                lambda.list = init$lambda.list)
        } else {
            # initialize the parameters, W can be supplied or initialized
            if (verbose)
                logmsg("Initialize W, H, b, Lambda ...")
            params.list = initialize_params(X.list = X.list, r = r, gamma = gamma,
                W = init$W, verbose = verbose, n.cores = n.cores)
        }

        obj = objective_func(X.list = X.list, W = params.list$W, H.list = params.list$H.list,
            lambda.list = params.list$lambda.list, b.list = params.list$b.list, gamma = gamma,
            n.cores = n.cores)
        if (verbose)
            logmsg("Objective for initialization = ", obj)

        # initialize
        delta = Inf
        converge = F
        obj.history = obj
        delta.history = NULL
        deltaw.history = NULL


        # solve 4 set of parameters iteratively
        for (iter in 1:max.niter) {

            obj.old = obj  # save or calculating the gap
            w.old = params.list$W

            # random permute update order to ensure convergence
            params.to.update.list = sample(c("W", "lambda", "b", "H"), 4, replace = F)
            if (verbose)
                logmsg("iter ", iter, ", update by: ", paste(params.to.update.list,
                  collapse = "->"))

            for (params.to.update in params.to.update.list) {

                params.list = solve_subproblem(params.to.update = params.to.update,
                  X.list = X.list, W = params.list$W, H.list = params.list$H.list,
                  b.list = params.list$b.list, lambda.list = params.list$lambda.list,
                  gamma = gamma, verbose = verbose, n.cores = n.cores)
            }

            obj = objective_func(X.list = X.list, W = params.list$W, H.list = params.list$H.list,
                lambda.list = params.list$lambda.list, b.list = params.list$b.list,
                gamma = gamma, n.cores = n.cores)
            obj.history = c(obj.history, obj)

            # relative difference of objective function
            delta = abs(obj - obj.old)/mean(c(obj, obj.old))
            delta.history = c(delta.history, delta)
            deltaw.history = c(deltaw.history, norm(w.old - params.list$W)/norm(w.old))

            if (verbose)
                logmsg("iter ", iter, ", objective=", obj, ", delta=diff/obj = ",
                  delta)

            # check if converge
            if (delta < tol) {
                logmsg("Converge at iter ", iter, ", obj delta = ", delta)
                converge = T
                break
            }
        }

        if (obj < obj.best) {
            # update to save the best results
            obj.best = obj
            obj.history.best = obj.history
            params.list.best = params.list
            niter.best = iter
            delta.best = delta
            delta.history.best = delta.history
            deltaw.history.best = deltaw.history
            converge.best = converge
            seed.best = seed + rep - 1
        }
    }
    time.elapsed = difftime(time1 = Sys.time(), time2 = time.start, units = "auto")

    if (verbose) {
        logmsg("Finised in ", time.elapsed, " ", units(time.elapsed), "\nBest result with seed ",
            seed.best, ":\nConvergence status: ", converge.best, " at ", niter.best,
            " iterations\nFinal objective delta:", delta.best)
    }

    if (!is.null(rownames(X.list[[1]]))) {
        params.list.best$H.list = lapply(1:m, function(j) {
            H = params.list.best$H.list[[j]]
            rownames(H) = rownames(X.list[[j]])
            H
        })
        rownames(params.list.best$W) = colnames(X.list[[1]])
    }

    return(list(H.list = params.list.best$H.list, W = params.list.best$W, b.list = params.list.best$b.list,
        lambda.list = params.list.best$lambda.list, convergence = converge.best,
        obj = obj.best, obj.history = obj.history.best, delta = delta.best, delta.history = delta.history.best,
        deltaw.history.best = deltaw.history.best, niter = niter.best, params = list(gamma = gamma,
            max.niter = max.niter, tol = tol, nrep = nrep, seed = seed, n.cores = n.cores)))
}


#' Calculate the objective function
#'
#' \deqn{1/N \sum_j||X_j -(H_JW^T\Lambda_j  + 1_nj b_j^T)||_F^2 +gamma \sum_{l=1}^p(\sum_{j=1}^mnj/N \lambda_{jl}-1)^2}
#'
#' @param X.list a list of ncells-by-ngenes gene expression matrix
#' @param W ngenes-by-r numeric matrix.
#' @param lambda.list A list of scaling vector of size p (ngenes).
#' @param H.list A list of factor loading matrix of size ncells-by-r
#' @param b.list A list of shift vector of size p (ngenes).
#' @param gamma numeric scalar, parameter for the penalty term.
#' @param n.cores number of cores used for parallel computing
#'
#' @return numeric scalar, the value of the objective function
#' @export
objective_func <- function(X.list, W, lambda.list, H.list, b.list, gamma, n.cores) {
    obj.list = lapply(1:length(X.list), function(j) {
        tmp1 = H.list[[j]] %*% t(W)
        tmp2 = matrix(1, nrow = nrow(X.list[[j]]), ncol = 1) %*% b.list[[j]]
        sum((X.list[[j]] - tmp1 * rep(lambda.list[[j]], rep.int(nrow(H.list[[j]]),
            nrow(W))) - tmp2)^2)
    })

    nvec = sapply(X.list, nrow)
    ntotal = sum(nvec)

    penalty = gamma * sum((colSums(nvec/ntotal * do.call(rbind, lambda.list)) - 1)^2)  # *ntotal /ntotal
    obj = sum(do.call(c, obj.list))/ntotal + penalty

    return(obj)
}


#' Solve subproblems via coordinate descent, given the parameter to update
#'
#' Solve the subproblem given which parameter set to update. For each subproblem,
#' the exact solution is obtained
#'
#' @param params.to.update a characteristic scalar, choice of ('W','lambda','b','H'),
#' specifying which set of parameters to update
#' @param X.list a list of ncells-by-ngenes gene expression matrix
#' @param W ngenes-by-r numeric matrix.
#' @param lambda.list A list of scaling vector of size p (ngenes).
#' @param H.list A list of factor loading matrix of size ncells-by-r
#' @param b.list A list of shift vector of size p (ngenes).
#' @param gamma numeric scalar, parameter for the penalty term.
#' @param verbose boolean scalar, whether to show extensive program logs (default TRUE)
#' @param n.cores integer number of cores used for parallel computing
#'
#' @return a list containing updated parameters: W, H.list, lambda.list,  b.list
#' @export
solve_subproblem <- function(params.to.update = c("W", "lambda", "b", "H"), X.list,
    W, H.list, b.list, lambda.list, gamma, verbose = T, n.cores) {
    params.to.update = match.arg(params.to.update)
    m = length(X.list)

    if (params.to.update == "W") {
        W = solve_W(X.list = X.list, H.list = H.list, lambda.list = lambda.list,
            b.list = b.list, n.cores = n.cores)
    } else if (params.to.update == "lambda") {
        lambda.list = solve_lambda_list(X.list = X.list, W = W, H.list = H.list,
            b.list = b.list, gamma = gamma, n.cores = n.cores)
    } else if (params.to.update == "b") {
        b.list = lapply(1:m, function(j) solve_b(X.list[[j]], W = W, H = H.list[[j]],
            lambd = lambda.list[[j]], n.cores = n.cores))
    } else {
        H.list = lapply(1:m, function(j) solve_H(X = X.list[[j]], W = W, lambd = lambda.list[[j]],
            b = b.list[[j]], n.cores = n.cores))
    }

    return(list(W = W, lambda.list = lambda.list, b.list = b.list, H.list = H.list))
}



#' Random paramter initialization
#'
#' @param X.list a list of ncells-by-ngenes gene expression matrix
#' @param r scalar, dimensional of common factor matrix, which can be chosen as the rough number of
#' identifiable cells types in the joint population.
#' @param seed random seed used for random generation (default 0).
#'
#' @return a list containing initialized parameters: W, H.list, lambda.list,  b.list
#' @export
initialize_params_random <- function(X.list, r, seed = 0) {
    set.seed(seed)
    m = length(X.list)
    p = ncol(X.list[[1]])
    n.list = sapply(1:m, function(j) nrow(X.list[[j]]))

    H.list = lapply(1:m, function(j) matrix(runif(n.list[j] * r, min = 0, max = 1),
        nrow = n.list[j], ncol = r))
    W = matrix(runif(p * r, min = 0, max = 2), nrow = p)
    lambda.list = lapply(1:m, function(i) rep(1, p))
    b.list = lapply(1:m, function(i) rep(0, p))

    return(list(W = W, H.list = H.list, lambda.list = lambda.list, b.list = b.list))
}



#' Initialize parameters for data integration or transfer
#'
#' Initialize the non-negative factor loading with the label encoding by performing k-means on the
#' centered and scaled data matrix concatenated. Then the comman factor matrix W is initilizated
#' sequentially by nonnegative matrix fatorization using scaled but not centered data matrix.
#'
#' @param X.list a list of ncells-by-ngenes gene expression matrix
#' @param r scalar, dimensional of common factor matrix, which can be chosen as the rough number of
#' identifiable cells types in the joint population (default 15).
#' @param gamma numeric scalar, parameter for the penalty term. It is sufficient to chose a large
#' value (default 1e6)
#' @param W ngenes-by-r numeric matrix. Supplied if parameter initialization is provided (default NULL).
#' @param verbose boolean scalar, whether to show extensive program logs (default TRUE)
#' @param n.cores number of cores used for parallel computing
#'
#' @return a list containing initialized parameters: W, H.list, lambda.list,  b.list
#' @import parallel checkmate
#' @export
initialize_params <- function(X.list, r, gamma, W = NULL, verbose = TRUE, n.cores = parallel::detectCores() -
    4) {
    m = length(X.list)
    p = ncol(X.list[[1]])

    # initialize W
    if (is.null(W)) {
        X.list.scale = lapply(X.list, function(X) {
            X = scale(X, center = T, scale = T)
            X[is.na(X)] = 0
            X
        })

        if (sum(do.call(c, lapply(X.list.scale, nrow))) > 10000) {
            nstart = 1
        } else {
            nstart = 3
        }

        kmeans.out = kmeans(do.call(rbind, X.list.scale), centers = r, nstart = nstart,
            iter.max = 10)
        W = t(kmeans.out$centers)  # there are negative values in W, this is used to initiate H

        # initialize H nonegative ensured
        H.list = lapply(1:m, function(j) solve_H(X.list.scale[[j]], W = W, lambd = rep(1,
            p), b = rep(0, p), n.cores = n.cores))

        # update W, nonnegative ensured
        b.list.init = lapply(1:m, function(j) rep(0, p))
        lambda.list.init = lapply(1:m, function(j) rep(1, p))
        W = solve_W(X.list = X.list, H.list = H.list, lambda.list = lambda.list.init,
            b.list = b.list.init, n.cores = n.cores)
    } else {

        # check nonegative
        checkmate::assert_true(all(W >= 0))

        # initialize H nonegative ensured
        H.list = lapply(1:m, function(j) solve_H(X.list[[j]], W = W, lambd = rep(1,
            p), b = rep(0, p), n.cores = n.cores))
    }

    # initialize lambda
    b.list.init = lapply(1:m, function(j) rep(0, p))
    lambda.list = solve_lambda_list(X.list = X.list, W = W, H.list = H.list, b.list = b.list.init,
        gamma = gamma, n.cores = n.cores)

    # initialize b
    b.list = lapply(1:m, function(j) solve_b(X = X.list[[j]], W = W, H = H.list[[j]],
        lambd = lambda.list[[j]], n.cores = n.cores))

    return(list(W = W, H.list = H.list, lambda.list = lambda.list, b.list = b.list))
}

#' Solve for factor loading paramters H
#'
#' \deqn{argmin_{H>=0} ||X- HW^T diag(lambd) - 1_n b^T||_F^2   s.t H1 = 1)}
#'
#' @param X ncells-by-ngenes gene expression matrix
#' @param W ngenes-by-r non-negative common factor matrix
#' @param lambd nonnegative numeric scalar, scaling associated with the dataset
#' @param b nonegative scalar, shift associated with tht edataset
#' @param n.cores number of cores used for the process
#'
#' @return ncells-by-r matrix, factor loading matrix H
#' @import checkmate parallel
#' @importFrom lsei pnnls
#' @export
solve_H <- function(X, W, lambd, b, n.cores) {
    # check size X n*p, W p*r, lambda p, b p
    checkmate::assert_true(all(c(ncol(X), nrow(W), length(lambd)) == rep(length(b),
        3)))

    A = lambd * W  #p*r
    H = do.call(rbind, lapply(1:nrow(X), function(i) {
        lsei::pnnls(a = A, b = X[i, ] - b, sum = 1)$x
    }))

    checkmate::assert_true(any(is.na(H)) == F)
    return(H)
}


#' Solve for nonnegative common factor matrix W
#'
#' \deqn{argmin_{W>=0} ||X- HW^T diag(lambd) - 1_n b^T||_F^2  }
#'
#' @param X.list A list of ncells-by-ngenes gene expression matrix.
#' @param H.list A list of factor loading matrix of size ncells-by-r
#' @param lambda.list A list of scaling vector of size p (ngenes).
#' @param b.list A list of shift vector of size p (ngenes).
#'
#' @return W ngenes-by-r common factor matrix shared among datasets
#' @import checkmate parallel
#' @importFrom lsei nnls
#'
#' @export
solve_W <- function(X.list, H.list, lambda.list, b.list, n.cores) {
    p = length(lambda.list[[1]])
    m = length(H.list)
    checkmate::assert_true(all(c(length(lambda.list), length(b.list)) == rep(m, 2)))

    nj.list = lapply(X.list, nrow)  # avoid repeated calculation

    W = do.call(rbind, parallel::mclapply(1:p, function(l) {
        A = do.call(rbind, lapply(1:m, function(j) lambda.list[[j]][l] * H.list[[j]]  # nj*r
))
        B = do.call(c, lapply(1:m, function(j) {
            X.list[[j]][, l] - b.list[[j]][l]  # nj*1
        }))

        lsei::nnls(a = A, b = B)$x
    }, mc.cores = n.cores))

    checkmate::assert_true(any(is.na(W)) == F)

    return(W)
}

#' Solve for dataset specific scalings lambda.list
#'
#' \deqn{argmin_{lambda_j>=0} ||X- HW^T diag(lambda_j) - 1_n b^T||_F^2 +
#' gamma * N * sum_{l=1}^p(\sum_{j=1}^mnj/N \lambda_{jl}-1)^2}
#'
#' @param X.list A list of ncells-by-ngenes gene expression matrix.
#' @param H.list A list of factor loading matrix of size ncells-by-r
#' @param W ngenes-by-r non-negative common factor matrix
#' @param b.list A list of shift vector of size p (ngenes).
#' @param gamma numeric scalar, parameter for the penalty term.
#' @param n.cores integer, number of cores used for parallel computation
#'
#' @return lambda.list A list of m scaling vector of size p (ngenes).
#' @import checkmate parallel
#' @importFrom lsei nnls
#'
#' @export
solve_lambda_list <- function(X.list, W, H.list, b.list, gamma, n.cores) {
    nvec = sapply(X.list, nrow)
    ntotal = sum(nvec)
    m = length(nvec)

    if (m > 1) {
        lambda.out = parallel::mclapply(1:nrow(W), function(l) {
            Ajl.list = lapply(1:m, function(j) {
                H.list[[j]] %*% W[l, ]  # nj * 1
            })
            Bjl.list = lapply(1:m, function(j) {
                X.list[[j]][, l] - b.list[[j]][l]
            })

            Amat <- diag(sapply(Ajl.list, function(Ajl) sum(Ajl^2))) + gamma * matrix(nvec,
                ncol = 1) %*% matrix(nvec, nrow = 1)/ntotal
            Bvec <- sapply(1:m, function(j) sum(Ajl.list[[j]] * Bjl.list[[j]])) +
                gamma * nvec

            lambd = lsei::nnls(a = Amat, b = Bvec)$x
            lambd[is.na(lambd)] = 0

            lambd
        }, mc.cores = n.cores)
        lambda.list = lapply(1:m, function(j) sapply(lambda.out, function(lambd) lambd[j]))
    } else {
        # only one dataset
        lambda.list = list(rep(1, nrow(W)))
    }
    return(lambda.list)
}


#' Solve for dataset specific shift b
#'
#' \deqn{argmin_{b_j>=0} ||X- HW^T diag(lambda_j) - 1_n b^T||_F^2 }
#'
#' @param X ncells-by-ngenes gene expression matrix
#' @param W ngenes-by-r non-negative common factor matrix
#' @param H ncells-by-r nonnegative factor loading matrix
#' @param lambd numeric scalar, scaling associated with the dataset
#' @param n.cores integer, number of cores used for parallel computation
#'
#' @return b shift vector of size p (ngenes).
#'
#' @import parallel
#' @export
solve_b <- function(X, W, H, lambd, n.cores) {

    b = do.call(c, parallel::mclapply(1:nrow(W), function(l) {
        max(0, mean(X[, l] - H %*% W[l, ] * lambd[l]))
    }, mc.cores = n.cores))

    return(b)
}

