#' Transfer the learned gene expression signature onto a new dataset
#'
#' Estimate the data specific parameters through Iterative Nonnegative Matrix Factorization (iNMF),
#' by minimizing the objective function \deqn{1/n||X -(HW^T\Lambda  + 1_n b^T)||_F^2}.
#'
#' @param Xtarget An ncells-by-ngenes gene expression matrix
#' @param Wref An ngenes-by-r reference low dimensional factor matrix
#' @param max.niter integer, max number of iterations (default 100).
#' @param tol numeric scalar, tolerance used in stopping criteria (default 1e-5).
#' @param nrep integer, number of repeated runs (to reduce effect of local optimum, default 1)
#' @param init a list of parameters for parameter initialization: lambda, b, H
#' @param verbose boolean scalar, whether to show extensive program logs (default TRUE)
#' @param seed random seed used (default 0)
#'
#' @return a list containing  \describe{
#'  \item{H} estimated factor loading matrix of size ncells-by-r
#'  \item{b}  estimated shift vector of size p (ngenes)
#'  \item{lambda} estimated scaling vector of size p (ngenes)
#'  \item{convergence} boolean, whether the algorithm converge
#'  \item{obj} numeric scalar, value of the objective function at convergence or when maximum iteration achieved
#'  \item{niter} integer, the iteration at convergence (or maximum iteration if not converge)
#'  \item{delta} numeric scalar, the relative difference in last objective update
#'  \item{params} list of parameters used for the algorithm: max.iter, tol, nrep
#' }
#'
#' @import checkmate
#' @export
CFITTransfer <- function(Xtarget, Wref, max.niter = 100, tol = 1e-05, init = NULL, seed = 0,
                         verbose = T, n.cores = parallel::detectCores() - 4) {


    checkmate::assert_true(ncol(Xtarget) == nrow(Wref))

    if (!is.null(rownames(Wref)) & !is.null(colnames(Xtarget))) {
        checkmate::assert_true(all(rownames(Wref) == colnames(Xtarget)))
    }
    n = nrow(Xtarget)

    # remove the genes that has na/inf/nan
    has.na = apply(Xtarget, 2, function(x) any(is.na(x)))
    has.inf = apply(Xtarget, 2, function(x) any(is.infinite(x)))
    good.genes = which(!(has.na | has.inf))
    Xtarget = Xtarget[, good.genes]
    Wref = Wref[good.genes, ]
    p = nrow(Wref)

    set.seed(seed)
    time.start = Sys.time()

    # if initial values of params provide
    if (all(c("lambda", "b", "H") %in% names(init))) {
        params.list = list(H = init$H, b = init$b, lambda = init$lambda)
    } else {
        # initialize
        if (verbose)
            logmsg("Initialize target H, b, Lambda ...")
        params.list = list(b = rep(0, p), lambda = rep(1, p), H = solve_H(Xtarget, W = Wref, lambd = rep(1, p), b = rep(0, p)))
    }

    obj = transfer_objective_func(X = Xtarget, W = Wref, H = params.list$H, lambd = params.list$lambda, b = params.list$b, n.cores=n.cores)
    converge = F

    for (iter in 1:max.niter) {
        obj.old = obj

        params.to.update.list = sample(c("lambda", "b", "H"), 3, replace = F)
        if (verbose)
            logmsg("iter ", iter, ", update by: ", paste(params.to.update.list, collapse = "->"))

        for (params.to.update in params.to.update.list) {

            params.list = transfer_solve_subproblem(params.to.update = params.to.update, X = Xtarget, W = Wref,
                                                    H = params.list$H, b = params.list$b,
                                                    lambd = params.list$lambda, verbose = verbose, n.cores=n.cores)
        }

        obj = transfer_objective_func(X = Xtarget, W = Wref, H = params.list$H,
                                      lambd = params.list$lambda, b = params.list$b, n.cores=n.cores)

        delta = abs(obj - obj.old)/mean(c(obj, obj.old))
        if (verbose)
            logmsg("iter ", iter, ", objective=", obj, ", delta=diff/obj = ", delta)

        # check if converge
        if (delta < tol) {
            if (verbose)
                logmsg("Converge at iter ", iter, ", obj delta = ", delta)
            converge = T
            break
        }

    }
    time.elapsed = difftime(time1 = Sys.time(), time2 = time.start, units = "auto")

    if (verbose) {
        logmsg("Finised in ", time.elapsed, " ", units(time.elapsed), "\n", "Convergence status: ",
               converge, " at ", iter, " iterations\nFinal objective delta:", delta)
    }

    if (!is.null(rownames(Xtarget))) {
        rownames(params.list$H) = rownames(Xtarget)
    }

    if (!is.null(colnames(Xtarget))) {
        names(params.list$b) = colnames(Xtarget)
        names(params.list$lambda) = colnames(Xtarget)
    }

    return(list(H = params.list$H, b = params.list$b, lambda = params.list$lambda,
                convergence = converge, obj = obj, niter = iter, delta = delta,
                params = list(max.niter = max.niter, tol = tol, Wref = Wref)))

}

#' Calculate the objective function
#'
#' \deqn{G(H, \Lambda, b; Wref):= \|X -(H Wref^T\Lambda +  1 d^T)\|_F^2}
#'
#' @param X An ncells-by-ngenes gene expression matrix
#' @param W An ngenes-by-r reference low dimensional factor matrix
#' @param H A factor loading matrix of size ncells-by-r
#' @param lambd A numeric vecor of the scaling
#' @param b A numeric shift vector of size p (ngenes).
#'
#' @return numeric scalar, the value of the objective function
transfer_objective_func <- function(X, W, H, lambd, b, n.cores) {
    n = nrow(X)
    tmp = t(X) - lambd * W %*% t(H) - b  # p by n
    obj = sum(tmp^2)/n

    return(obj)
}


#' Solve subproblems via coordinate descent, given the parameter to update
#'
#' Solve the subproblem given which parameter set to update. For each subproblem,
#' the exact solution is obtained
#'
#' @param params.to.update a characteristic scalar, choice of (lambda','b','H'),
#' specifying which set of parameters to update
#' @param X An ncells-by-ngenes gene expression matrix
#' @param W An ngenes-by-r reference low dimensional factor matrix
#' @param H A factor loading matrix of size ncells-by-r
#' @param b A numeric shift vector of size p (ngenes).
#' @param lambda A numeric scaling vector of size p (ngenes)
#' @param verbose boolean scalar, whether to show extensive program logs (default TRUE)
#'
#' @return a list containing updated parameters: H, lambda,  b
transfer_solve_subproblem <- function(params.to.update = c("lambda", "b", "H"),
                                      X, W, H, b, lambd, verbose = T, n.cores) {
    params.to.update = match.arg(params.to.update)

    if (params.to.update == "lambda") {
        lambd = transfer_solve_lambda(X = X, W = W, H = H, b = b, n.cores=n.cores)
    } else if (params.to.update == "b") {
        b = solve_b(X = X, W = W, H = H, lambd = lambd, n.cores=n.cores)
    } else {
        H = solve_H(X = X, W = W, lambd = lambd, b = b, n.cores=n.cores)
    }
    return(list(lambda = lambd, b = b, H = H))
}


#' Solve for dataset specific scalings lambda for target data
#'
#' \deqn{argmin_{lambda>=0} ||X- HW^T diag(lambd) - 1_n b^T||_F^2}
#'
#' @param X An ncells-by-ngenes gene expression matrix
#' @param W An ngenes-by-r reference low dimensional factor matrix
#' @param H A factor loading matrix of size ncells-by-r
#' @param b A numeric shift vector of size p (ngenes).
#'
#' @return numeric vector, scaling of target data with respect to the reference factor matrix
#' @import parallel
transfer_solve_lambda <- function(X, W, H, b, n.cores) {

    n = nrow(X)

    A = H %*% t(W)
    xmean = mean(c(X))

    lambd = do.call(c, parallel::mclapply(1:ncol(X), function(l) {
        a_l = A[, l]
        b_l = X[, l] - b[l]

        max(0, (a_l %*% b_l)/(a_l %*% a_l))
    }, mc.cores = n.cores))

    return(lambd)
}


