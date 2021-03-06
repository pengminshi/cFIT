#' Integration of multiple data source.
#'
#' Solve the model parameters through Iterative Nonnegative Matrix Factorization (iNMF),
#' by minimizing the objective function \deqn{1/N \sum_j||X_j -(H_jW^T\Lambda_j  + 1_{n_j} b_j^T)||_F^2} with penalties.
#'
#' @param X.list a list of m ncells-by-ngenes, gene expression matrices from m data sets
#' @param r scalar, dimension of common factor matrix, which can be chosen as the rough number of
#' identifiable cells types in the joint population (default 15).
#' @param max.niter integer, max number of iterations (default 100).
#' @param tol numeric scalar, tolerance used in stopping criteria (default 1e-5).
#' @param nrep integer, number of repeated runs (to reduce effect of local optimum, default 1)
#' @param init a list of parameters for parameter initialization. The list either contains all
#' parameter sets: W,lambda.list, b.list, H.list, or only W will be used if provided (default NULL).
#' @param future.plan plan for future parallel computation, can be chosen from 'sequential','transparent','multicore','multisession' and 'cluster'. Default is 'sequential'. Note that Rstudio does not support 'multicore'.
#' @param workers additional parameter for \code{future::plan()}, in cases of 'multicore','multisession' and 'cluster'.
#' @param verbose boolean scalar, whether to show extensive program logs (default TRUE)
#' @param seed random seed used (default 0)
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
#'  \item{params}{list of parameters used for the algorithm: max.iter, tol, nrep}
#' }
#'
#' @import checkmate parallel
#' @importFrom future.apply future_apply
#' @importFrom future plan
#' @export
CFITIntegrate <- function(X.list, r = 15, max.niter = 100, tol = 1e-05,
    nrep = 1, init = NULL,
    future.plan=c('sequential','transparent','multicore','multisession','cluster'),
    workers = parallel::detectCores() - 1,
    verbose = T, seed = 0) {

    env.plan = future::plan()
    future.plan = match.arg(future.plan)
    future::plan(future.plan)

    if (verbose){
        logmsg("Run cFIT with ", future.plan,  " plan ...")
        if (future.plan %in% c('multicore', 'multisession', 'cluster')){
            logmsg('Worker: ', workers)
        }
    }

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
            params.list = initialize_params(X.list = X.list, r = r,
                W = init$W, verbose = verbose)
        }

        obj = objective_func(X.list = X.list, W = params.list$W, H.list = params.list$H.list,
            lambda.list = params.list$lambda.list, b.list = params.list$b.list)
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
            params.to.update.list = sample(c("W", "lambda", "H"), replace = F)
            if (verbose)
                logmsg("iter ", iter, ", update by: ", paste(params.to.update.list,
                  collapse = "->"))

            for (params.to.update in params.to.update.list) {

                params.list = solve_subproblem(params.to.update = params.to.update,
                  X.list = X.list, W = params.list$W, H.list = params.list$H.list,
                  b.list = params.list$b.list, lambda.list = params.list$lambda.list,
                  verbose = verbose)
            }

            obj = objective_func(X.list = X.list, W = params.list$W, H.list = params.list$H.list,
                lambda.list = params.list$lambda.list, b.list = params.list$b.list)
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

    future::plan(env.plan)

    return(list(H.list = params.list.best$H.list, W = params.list.best$W, b.list = params.list.best$b.list,
        lambda.list = params.list.best$lambda.list, convergence = converge.best,
        obj = obj.best, obj.history = obj.history.best, delta = delta.best, delta.history = delta.history.best,
        deltaw.history.best = deltaw.history.best, niter = niter.best, params = list(#gamma = gamma,
            max.niter = max.niter, tol = tol, nrep = nrep, seed = seed)))
}


#' Calculate the objective function
#'
#' \deqn{1/N \sum_j||X_j -(H_JW^T\Lambda_j  + 1_{n_j} b_j^T)||_F^2}
#'
#' @param X.list a list of ncells-by-ngenes gene expression matrix
#' @param W ngenes-by-r numeric matrix.
#' @param lambda.list A list of scaling vector of size p (ngenes).
#' @param H.list A list of factor loading matrix of size ncells-by-r
#' @param b.list A list of shift vector of size p (ngenes).
#'
#' @return numeric scalar, the value of the objective function
#' @export
objective_func <- function(X.list, W, lambda.list, H.list, b.list) {
    if (is.null(H.list)){
        H.list = lapply(1:length(X.list), function(j) solve_H(X = X.list[[j]], W = W, lambd = lambda.list[[j]],
                                                 b = b.list[[j]]))
    }
    obj.list = lapply(1:length(X.list), function(j) {
        tmp1 = H.list[[j]] %*% t(W)
        tmp2 = matrix(1, nrow = nrow(X.list[[j]]), ncol = 1) %*% b.list[[j]]
        sum((X.list[[j]] - tmp1 * rep(lambda.list[[j]], rep.int(nrow(H.list[[j]]),
            nrow(W))) - tmp2)^2)
    })

    nvec = sapply(X.list, nrow)
    ntotal = sum(nvec)

    obj = sum(do.call(c, obj.list))/ntotal

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
#' @param verbose boolean scalar, whether to show extensive program logs (default TRUE)
#'
#' @return a list containing updated parameters: W, H.list, lambda.list,  b.list
#' @export
solve_subproblem <- function(params.to.update = c("W", "lambda", "H"), X.list,
    W, H.list, b.list, lambda.list, verbose = T) {
    params.to.update = match.arg(params.to.update)
    m = length(X.list)

    if (params.to.update == "W") {
        W = solve_W(X.list = X.list, H.list = H.list, lambda.list = lambda.list,
            b.list = b.list)
        b.list = lapply(1:m, function(j) solve_b(X.list[[j]], W = W, H = H.list[[j]],
                                                 lambd = lambda.list[[j]]))
    } else if (params.to.update == "lambda") {
        lambda.list = solve_lambda_list(X.list = X.list, W = W, H.list = H.list,
            b.list = b.list)
    } else {
        H.list = lapply(1:m, function(j) solve_H(X = X.list[[j]], W = W, lambd = lambda.list[[j]],
            b = b.list[[j]]))
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
#' @param W ngenes-by-r numeric matrix. Supplied if parameter initialization is provided (default NULL).
#' @param verbose boolean scalar, whether to show extensive program logs (default TRUE)
#'
#' @return a list containing initialized parameters: W, H.list, lambda.list,  b.list
#' @import parallel checkmate
#' @export
initialize_params <- function(X.list, r, W = NULL, verbose = TRUE) {
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
            p), b = rep(0, p)))

        # update W, nonnegative ensured
        b.list.init = lapply(1:m, function(j) rep(0, p))
        lambda.list.init = lapply(1:m, function(j) rep(1, p))
        W = solve_W(X.list = X.list, H.list = H.list, lambda.list = lambda.list.init,
            b.list = b.list.init)
    } else {

        # check nonegative
        checkmate::assert_true(all(W >= 0))

        # initialize H nonegative ensured
        H.list = lapply(1:m, function(j) solve_H(X.list[[j]], W = W, lambd = rep(1,
            p), b = rep(0, p)))
    }

    # initialize lambda
    b.list.init = lapply(1:m, function(j) rep(0, p))
    lambda.list = solve_lambda_list(X.list = X.list, W = W, H.list = H.list, b.list = b.list.init)

    # initialize b
    b.list = lapply(1:m, function(j) solve_b(X = X.list[[j]], W = W, H = H.list[[j]],
        lambd = lambda.list[[j]]))

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
#'
#' @return ncells-by-r matrix, factor loading matrix H
#' @import checkmate parallel
#' @importFrom lsei pnnls
#' @export
solve_H <- function(X, W, lambd, b) {
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
solve_W <- function(X.list, H.list, lambda.list, b.list) {
    p = length(lambda.list[[1]])
    m = length(H.list)
    checkmate::assert_true(all(c(length(lambda.list), length(b.list)) == rep(m, 2)))

    nj.list = lapply(X.list, nrow)  # avoid repeated calculation

    W = do.call(rbind, future.apply::future_lapply(1:p, function(l) {
        A = do.call(rbind, lapply(1:m, function(j) lambda.list[[j]][l] * H.list[[j]]  # nj*r
))
        B = do.call(c, lapply(1:m, function(j) {
            X.list[[j]][, l] - b.list[[j]][l]  # nj*1
        }))

        lsei::nnls(a = A, b = B)$x
    }))

    checkmate::assert_true(any(is.na(W)) == F)

    return(W)
}

#' Solve for dataset specific scalings lambda.list
#'
#' \deqn{argmin_{lambda_j>=0} ||X- HW^T diag(lambda_j) - 1_n b^T||_F^2}
#'
#' @param X.list A list of ncells-by-ngenes gene expression matrix.
#' @param H.list A list of factor loading matrix of size ncells-by-r
#' @param W ngenes-by-r non-negative common factor matrix
#' @param b.list A list of shift vector of size p (ngenes).
#'
#' @return lambda.list A list of m scaling vector of size p (ngenes).
#' @import checkmate parallel
#' @importFrom lsei nnls
#'
#' @export
solve_lambda_list <- function(X.list, W, H.list, b.list) {
    nvec = sapply(X.list, nrow)
    ntotal = sum(nvec)
    m = length(nvec)
    p = nrow(W)

    if (m > 1) {
        lambda.list = lapply(1:m, function(j){
            lambd = future.apply::future_sapply(1:p, function(l) {
                y = X.list[[j]][, l] - b.list[[j]][l]
                x = H.list[[j]] %*% W[l, ]
                xx = sum(x * x)
                xy = sum(x * y)
                if (xx == 0 | xy < 0) {
                    return(0)
                }
                return(xy / xx)
            })
        })

        # calculate the scaling for each gene
        scale.per.gene = sapply(1:p, function(l) {
            lambdas = sapply(1:m, function(j) lambda.list[[j]][l])
            lambda.sums = sum(lambdas * nvec)
            if (lambda.sums == 0){
                return(1)
            }
            return(ntotal / lambda.sums)
        })

        # rescaled lambda.list
        lambda.list = lapply(lambda.list, function(lambd) {
            lambd = lambd * scale.per.gene
            # lambd[is.na(lambd)] = 1
            return(lambd)
        })
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
#' @param b.gamma tunning parameter for the L2 panelty on b.gamma.
#'
#' @return b shift vector of size p (ngenes).
#'
#' @import parallel
#' @export
solve_b <- function(X, W, H, lambd, b.gamma = 0.1) {

    b = do.call(c, future.apply::future_lapply(1:nrow(W), function(l) {
        mean(X[, l] - H %*% W[l, ] * lambd[l])/ (1+b.gamma)
    }))
    # b = rep(0, nrow(W))

    return(b)
}

