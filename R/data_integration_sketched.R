#' Integration of multiple data source by fast cFIT with sketching and Stochastic Proximal Point method (SPP).
#'
#' Solve the model parameters through Iterative Nonnegative Matrix Factorization (iNMF),
#' by minimizing the sketched objective function \deqn{1/\tilde{N} \sum_j||SX_j -(SH_JW^T\Lambda_j  + S1_nj b_j^T)||_F^2 +
#' gamma \sum_{l=1}^p(\sum_{j=1}^m\tilde{n}_j/\tilde{N} \lambda_{jl}-1)^2}, with additional penalty for SPP.
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
#' @param subsample.prop a scalar between 0 and 1. smaller proportion with results in fast computation but less
#' accurate results. By default the value is set to  \code{min(5*10^4/ntotal, 1)}
#' @param weight.list weights for performing weighted subsampling sketching. Note that the weight.list is a list of weights
#' per batch. The weights for each batch is a vector of nonnegative values of the same size as the number of cells in the batch.
#' To use the Statistical leverage score as weights, run
#' @param future.plan plan for future parallel computation, can be chosen from 'sequential','transparent','multicore','multisession' and 'cluster'. Note that Rstudio does not support 'multicore'.
#' @param workers additional parameter for \code{future::plan()}, in cases of 'multicore','multisession' and 'cluster'.
#' \code{weight.list = lapply(1:length(X.list), function(j) statistical_leverage_score(X.list[[j]], k=r))}
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
#'  \item{obj/history}{a numeric vector, value of the objective function per iteration}
#'  \item{deltaw}{numeric, the relative change in W (common factor matrix) measured by Frobenious norm}
#'  \item{deltaw.history}{a vector of numeric values, the relative change in W (common factor matrix) per iteration.}
#'  \item{niter}{integer, the iteration at convergence (or maximum iteration if not converge)}
#'  \item{params}{list of parameters used for the algorithm: gamma, max.iter, tol, nrep, subsample.prop, weight.list}
#' }
#'
#' @import checkmate parallel
#' @export
CFITIntegrate_sketched <- function(X.list, r = 15, max.niter = 100, tol = 1e-04,
    gamma = 1e+06, nrep = 1, init = NULL, subsample.prop = NULL, weight.list = NULL,
    future.plan=c('multicore','sequential','transparent','multisession','cluster'),
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

    # proportion to be subsampled for updates
    if (is.null(subsample.prop)) {
        subsample.prop = min(5 * 10^4/ntotal, 1)
        logmsg("Use subsample proportion ", subsample.prop)
    }

    # check the validity of weights
    if (!is.null(weight.list)) {
        checkmate::assert_true(length(weight.list) == length(X.list))
    }

    # save the best results with minimum objective for output
    obj.best = Inf
    obj.history.best = NULL
    deltaw.history.best = NULL
    params.list.best = NULL
    niter.best = NULL
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

            if (!is.null(weight.list)) {
                X.list.sub = subsample(X.list, subsample.prop, weight.list = weight.list)
            } else {
                X.list.sub = subsample(X.list, subsample.prop)
            }

            params.list = initialize_params(X.list = X.list.sub, r = r, gamma = gamma,
                W = init$W, verbose = verbose)
            # params.list = initialize_params_random(X.list = X.list, r = r) # random
            # initialization
        }

        obj = objective_func(X.list = X.list.sub, W = params.list$W, H.list = params.list$H.list,
            lambda.list = params.list$lambda.list, b.list = params.list$b.list, gamma = gamma)
        if (verbose)
            logmsg("Objective for initialization = ", obj)

        # initialize
        converge = F
        obj.history = obj
        deltaw.history = NULL

        # solve 4 set of parameters iteratively
        for (iter in 1:max.niter) {

            w.old = params.list$W
            params.list.last = params.list

            if (!is.null(weight.list)) {
                X.list.sub = subsample(X.list, subsample.prop, weight.list = weight.list)
            } else {
                X.list.sub = subsample(X.list, subsample.prop)
            }

            # estimate the membership matrix first random permute update order for the other
            # three sets of parameters
            params.to.update.list = c("H", sample(c("W", "lambda", "b"), replace = F))
            if (verbose)
                logmsg("iter ", iter, ", update by: ", paste(params.to.update.list,
                  collapse = "->"))

            for (params.to.update in params.to.update.list) {
                # logmsg('test ----- > str(param.list)', str(params.list))
                params.list = solve_subproblem_penalized(params.to.update = params.to.update,
                  X.list = X.list.sub, W = params.list$W, H.list = params.list$H.list,
                  b.list = params.list$b.list, lambda.list = params.list$lambda.list,
                  gamma = gamma, iter = iter, params.list.last = params.list.last,
                  verbose = verbose)
            }
            # logmsg('test----> str(params.list)', str(params.list))
            obj = objective_func(X.list = X.list.sub, W = params.list$W, H.list = params.list$H.list,
                lambda.list = params.list$lambda.list, b.list = params.list$b.list,
                gamma = gamma)
            obj.history = c(obj.history, obj)
            # relative difference of objective function
            deltaw = norm(params.list$W - w.old)/norm(w.old)
            deltaw.history = c(deltaw.history, deltaw)
            if (verbose)
                logmsg("iter ", iter, ", objective=", obj, ", delta_w=", deltaw)

            # check if converge
            if (deltaw < tol) {
                logmsg("Converge at iter ", iter, ", deltaw = ", deltaw)
                converge = T
                break
            }
        }
        if (verbose){
            logmsg("Estimate the H for all samples")
        }
        params.list = solve_subproblem(params.to.update = "H", X.list = X.list, W = params.list$W,
            H.list = params.list$H.list, b.list = params.list$b.list, lambda.list = params.list$lambda.list,
            gamma = gamma, verbose = verbose)

        if (!is.null(rownames(X.list[[1]]))) {
            params.list$H.list = lapply(1:m, function(j) {
                H = params.list$H.list[[j]]
                rownames(H) = rownames(X.list[[j]])
                H
            })
            rownames(params.list$W) = colnames(X.list[[1]])
        }
        if (verbose){
            logmsg("Objective using all samples")
        }
        obj = objective_func(X.list = X.list, W = params.list$W, H.list = params.list$H.list,
            lambda.list = params.list$lambda.list, b.list = params.list$b.list, gamma = gamma)

        if (obj < obj.best) {
            # update to save the best results
            obj.best = obj
            obj.history.best = obj.history
            params.list.best = params.list
            niter.best = iter
            deltaw.best = deltaw
            deltaw.history.best = deltaw.history
            converge.best = converge
            seed.best = seed + rep - 1
        }
    }

    time.elapsed = difftime(time1 = Sys.time(), time2 = time.start, units = "min")
    if (verbose) {
        logmsg("Finised in ", time.elapsed, " ", units(time.elapsed), "\nBest result with seed ",
            seed.best, ":\nConvergence status: ", converge.best, " at ", niter.best,
            " iterations\nFinal objective delta:", deltaw.best)
    }

    future::plan(env.plan)
    return(list(H.list = params.list.best$H.list, W = params.list.best$W, b.list = params.list.best$b.list,
        lambda.list = params.list.best$lambda.list, convergence = converge.best,
        obj = obj.best, obj.history = obj.history.best, deltaw = deltaw.best, deltaw.history = deltaw.history.best,
        niter = niter.best, params = list(gamma = gamma, max.niter = max.niter, tol = tol,
            nrep = nrep, subsample.prop = subsample.prop, weight.list = weight.list)))
}


#' Subsample from a list of data matrix
#'
#' Subsample the samples in each batch, eighter weighted with the input weight.list or unweighted.
#'
#' @param X.list a list of m ncells-by-ngenes, gene expression matrices from m data sets
#' @param subsample.prop a scalar between 0 and 1. smaller proportion with results in fast computation but less accurate results.
#' @param min.samples integer, minimum number of samples from each batch (20 by default).
#' @param weight.list weights for performing weighted subsampling sketching. Note that the weight.list is a list of weights
#' per batch. The weights for each batch is a vector of nonnegative values of the same size as the number of cells in the batch.
#' No weight by default (weight.list=NULL).
#'
#' @return a list of subsampled datasets
#' @export
subsample <- function(X.list, subsample.prop, min.samples = 20, weight.list = NULL) {
    X.list.sub = lapply(1:length(X.list), function(j) {

        x = X.list[[j]]
        n = nrow(x)
        if (is.null(weight.list)) {
            x[sample(1:n, round(max(n * subsample.prop, min(n, min.samples)))), ]
        } else {
            w = weight.list[[j]]
            w = w/sum(w)
            x[sample(1:n, round(max(n * subsample.prop, min(n, min.samples))), prob = w),
                ]
        }
    })

    return(X.list.sub)
}


#' Solve subproblems via coordinate descent with Stochastic proximal point method
#'
#' Solve the subproblem given which parameter set to update. For each subproblem,
#' the exact solution is obtained. The sketched objective function is solve w.r.t the set
#' of parameter with additional penalty in the form of \deqn{\frac{1}{\mu_t} \|w-w^{}t-1\|},
#' where the step size \deqn{mu_t = 0.1\times iter}.
#'
#' @param params.to.update a characteristic scalar, choice of ('W','lambda','b','H'),
#' specifying which set of parameters to update
#' @param X.list a list of ncells-by-ngenes gene expression matrix
#' @param W ngenes-by-r numeric matrix.
#' @param lambda.list A list of scaling vector of size p (ngenes).
#' @param H.list A list of factor loading matrix of size ncells-by-r
#' @param b.list A list of shift vector of size p (ngenes).
#' @param gamma numeric scalar, parameter for the penalty term.
#' @param iter integer, the current iteration
#' @param params.list.last The parameters from last iteration, for SPP
#' @param verbose boolean scalar, whether to show extensive program logs (default TRUE)
#'
#' @return a list containing updated parameters: W, H.list, lambda.list,  b.list
#'
#' @export
solve_subproblem_penalized <- function(params.to.update = c("W", "lambda", "b", "H"),
    X.list, W, H.list, b.list, lambda.list, gamma, iter, params.list.last, verbose = T) {
    params.to.update = match.arg(params.to.update)
    m = length(X.list)

    if (params.to.update == "W") {
        W = solve_W_penalized(X.list = X.list, H.list = H.list, lambda.list = lambda.list,
            b.list = b.list, W.old = params.list.last$W, iter = iter)
    } else if (params.to.update == "lambda") {
        lambda.list = solve_lambda_list_penalized(X.list = X.list, W = W, H.list = H.list,
            b.list = b.list, gamma = gamma, lambda.list.old = params.list.last$lambda.list,
            iter = iter)
    } else if (params.to.update == "b") {
        b.list = lapply(1:m, function(j) solve_b_penalized(X.list[[j]], W = W, H = H.list[[j]],
            lambd = lambda.list[[j]], b.old = params.list.last$b.list[[j]], iter = iter))
    } else {
        H.list = lapply(1:m, function(j) solve_H(X = X.list[[j]], W = W, lambd = lambda.list[[j]],
            b = b.list[[j]]))
    }

    return(list(W = W, lambda.list = lambda.list, b.list = b.list, H.list = H.list))
}


#' Solve for nonnegative common factor matrix W with Stochastic Porximal Point method
#'
#' \deqn{argmin_{W>=0} ||X- HW^T diag(lambd) - 1_n b^T||_F^2  += \frac{1}{2\mu_t}\|W-W^{t-1}\|_F/ }
#'
#' @param X.list A list of ncells-by-ngenes gene expression matrix.
#' @param H.list A list of factor loading matrix of size ncells-by-r
#' @param lambda.list A list of scaling vector of size p (ngenes).
#' @param b.list A list of shift vector of size p (ngenes).
#' @param W.old W from last iteration.
#' @param iter current iteration, for calculating the step size.
#'
#' @return W ngenes-by-r common factor matrix shared among datasets
#' @import checkmate parallel
#' @importFrom lsei nnls
#' @export
solve_W_penalized <- function(X.list, H.list, lambda.list, b.list, W.old, iter) {
    p = length(lambda.list[[1]])
    m = length(H.list)
    checkmate::assert_true(all(c(length(lambda.list), length(b.list)) == rep(m, 2)))

    nj.list = lapply(X.list, nrow)  # avoid repeated calculation
    n = sum(do.call(c, nj.list))
    r = ncol(W.old)

    W = do.call(rbind, future.apply::future_lapply(1:p, function(l) {
        A = do.call(rbind, lapply(1:m, function(j) lambda.list[[j]][l] * H.list[[j]]  # nj*r
))  # n * r
        B = do.call(c, lapply(1:m, function(j) {
            X.list[[j]][, l] - b.list[[j]][l]  # nj*1
        }))  # n * 1

        # add the panelty
        mu = iter * 0.1  # penalty parameter
        A = rbind(A * sqrt(1/n), sqrt(mu/r) * diag(rep(1, r)))
        B = c(B * sqrt(1/n), W.old[l, ] * sqrt(mu/r))

        lsei::nnls(a = A, b = B)$x
    }))

    checkmate::assert_true(any(is.na(W)) == F)

    return(W)
}


#' Solve for dataset specific scalings lambda.list with Stochastic Porximal Point method
#'
#' @param X.list A list of ncells-by-ngenes gene expression matrix.
#' @param H.list A list of factor loading matrix of size ncells-by-r
#' @param W ngenes-by-r non-negative common factor matrix
#' @param b.list A list of shift vector of size p (ngenes).
#' @param gamma numeric scalar, parameter for the penalty term.
#' @param lambda.list.old lambda.list from last iteration
#' @param iter current iteration, for calculating the step size.
#'
#' @return lambda.list A list of m scaling vector of size p (ngenes).
#' @import checkmate parallel
#' @importFrom lsei nnls
#' @export
solve_lambda_list_penalized <- function(X.list, W, H.list, b.list, gamma, lambda.list.old,
    iter) {
    nvec = sapply(X.list, nrow)
    ntotal = sum(nvec)
    m = length(nvec)
    lambda.old.mat = do.call(rbind, lambda.list.old)  # m * p

    if (m > 1) {
        lambda.out = future.apply::future_lapply(1:nrow(W), function(l) {
            Ajl.list = lapply(1:m, function(j) {
                H.list[[j]] %*% W[l, ]  # nj * 1
            })
            Bjl.list = lapply(1:m, function(j) {
                X.list[[j]][, l] - b.list[[j]][l]
            })

            Amat <- diag(sapply(Ajl.list, function(Ajl) sum(Ajl^2))) + gamma * matrix(nvec,
                ncol = 1) %*% matrix(nvec, nrow = 1)/ntotal  # m*m matrix
            Bvec <- sapply(1:m, function(j) sum(Ajl.list[[j]] * Bjl.list[[j]])) +
                gamma * nvec

            # add penalty
            mu = iter * 0.1
            Amat = rbind(Amat, mu * diag(rep(1, m)))
            Bvec = c(Bvec, mu * lambda.old.mat[, l])

            lambd = lsei::nnls(a = Amat, b = Bvec)$x
            lambd[is.na(lambd)] = 0

            lambd
        })
        lambda.list = lapply(1:m, function(j) sapply(lambda.out, function(lambd) lambd[j]))
    } else {
        # only one dataset
        lambda.list = list(rep(1, nrow(W)))
    }
    return(lambda.list)
}


#' Solve for dataset specific shift b with Stochastic Porximal Point method
#'
#' @param X ncells-by-ngenes gene expression matrix
#' @param W ngenes-by-r non-negative common factor matrix
#' @param H ncells-by-r nonnegative factor loading matrix
#' @param lambd numeric scalar, scaling associated with the dataset
#' @param b.old b vector  from last iteration
#' @param iter current iteration, for calculating the step size.
#'
#' @return b shift vector of size p (ngenes).
#' @import parallel
#' @importFrom lsei nnls
#' @export
solve_b_penalized <- function(X, W, H, lambd, b.old, iter) {

    n = nrow(X)
    b = do.call(c, future.apply::future_lapply(1:nrow(W), function(l) {
        A = matrix(1, nrow = n, ncol = 1)
        B = X[, l] - H %*% W[l, ] * lambd[l]

        # add penalty
        mu = iter * 0.1
        A = rbind(A/sqrt(n), 1 * mu)
        B = c(B/sqrt(n), b.old[l] * mu)

        lsei::nnls(a = A, b = B)$x
    }))

    return(b)
}

#' Calculate the statistical leverage score
#'
#' @param A data matrix
#' @param k number of singular vectors used for calculating scores
#'
#' @return a vector of scores per samples
#'
#' @importFrom RSpectra svds
#' @import checkmate
#' @export
statistical_leverage_score <- function(A, k = NULL) {
    if (is.null(k)) {
        u = svd(A)$u
    } else {
        checkmate::assert_true(k < ncol(A))
        u = RSpectra::svds(A, k = k)$u
    }
    score = rowSums(u^2)

    return(score)
}

