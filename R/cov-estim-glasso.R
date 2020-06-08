#' Graphical Lasso (GLASSO) Covariance Estimation
#'
#' Computes the GLASSO estimator of the covariance matrix.
#'
#' @param data an nxp data matrix.
#' @param type a character, the type of matrix to be estimated. Possible values are c("cor", "cov"). Default value is "cor" for the correlation matrix.
#' @param rho a double or a sequence, the non-negative regularization parameter \eqn{\rho} for lasso. \eqn{\rho=0} means no regularization.
#' Can be a scalar (usual) or a symmetric p by p matrix, or a vector of length p. In the latter case, the penalty matrix has jkth element \eqn{\sqrt{\rho_j*\rho_k}}. Default value is NULL
#' and an optimal regularization is computed with the cross-validation (CV) procedure as in \insertCite{cvglassopackage;textual}{CovEstim}.
#' @param rho_len an integer, indicating the length for setting the rho sequence.
#' @param rho_seq_min_ratio a double, the minimum ratio for deviation between the rho values in the sequence.
#' @param nfolds an integer, indicating the number of folds for the CV. Default value is 5.
#' @param crit a character, indicating which selection criterion within the CV. Possible values are "loglik", "AIC" and "BIC". Default is set to "BIC".
#' @param pendiag_log a logical, indicating whether the diagonal of the sample covariance matrix is to be penalized (TRUE) or not (FALSE). Default value is FALSE.
#' @param start a character, specifying the start type of the glasso algorithm. Possible values are "warm" or "cold". Default value is "cold".
#' @param tol a double, indicating the tolerance for the glasso algorithm. Default value is set to 1e-05.
#' @param maxit an integer, indicating the maximum number of iterations for the glasso algorithm. Default value is set to 10000.
#' @param cores an integer, indicating how many cores should be used for the CV. Default value is 1. cores cannot be higher than the maximum number of cores of the processor in use.
#' @param seed an integer, the seed for the performed cross-validation. Default value is 1234.
#'
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the lasso penalty.
#' }
#'
#' @details The GLASSO estimator is elaborated in detail in \insertCite{friedman2008sparse;textual}{CovEstim}. More information on the functionality can be found in
#' \insertCite{glassopackage;textual}{CovEstim} and \insertCite{cvglassopackage;textual}{CovEstim}.
#'
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_glasso <- sigma_estim_glasso(sp_rets, type="cov", rho=0.0001)[[1]]
#'
#' @import glasso
#' @import CVglasso
#' @import doParallel
#' @import foreach
#' @importFrom Rdpack reprompt
#' @references
#'\insertAllCited
#'
#' @export sigma_estim_glasso
#'
sigma_estim_glasso <-
  function(data,
           rho = NULL,
           rho_len = 100,
           rho_seq_min_ratio = 0.001,
           type = "cor",
           nfolds = 5,
           crit = "loglik",
           pendiag_log = FALSE,
           start = "warm",
           tol = 1e-04,
           maxit = 10000,
           cores = 1,
           seed = 1234) {
    data <- as.matrix(data)

    if (type == "cov") {
      n <- dim(data)[1]
      centered <- apply(data, 2, function(x)
        x - mean(x))
      sigma_sample <- t(centered) %*% centered / (n - 1)

      if (is.null(rho)) {
        sigma_sample_diag0 <- sigma_sample
        diag(sigma_sample_diag0) <- 0
        rho_max <- max(abs(sigma_sample_diag0))
        rho_min <- rho_seq_min_ratio * rho_max
        rho <-
          10 ^ seq(log10(rho_min), log10(rho_max), length = rho_len)
      }
      results <-
        CVglasso::CVglasso(
          X = data,
          lam = rho,
          K = nfolds,
          crit.cv = crit,
          diagonal = pendiag_log,
          start = start,
          tol = tol,
          maxit = maxit,
          cores = cores,
          trace = "none"
        )
        rho <- results$Tuning[2]
        sigma_mat <- results$Sigma
        rownames(sigma_mat) <- colnames(data)
        colnames(sigma_mat) <- colnames(data)
        
    } else if (type == "cor") {
      results <-
        CVglasso_mod(
          X = data,
          lam = rho,
          nlam = rho_len,
          seed = seed,
          K = nfolds,
          crit.cv = crit,
          start = start,
          tol = tol,
          maxit = maxit,
          cores = cores,
          trace = "none"
        )
        
        rho <- results$Tuning[2]
        
        vola_mat <- diag(apply(data, 2, stats::sd, na.rm = TRUE))
        sigma_mat <- vola_mat%*%results$Sigma%*%vola_mat
        
        rownames(sigma_mat) <- colnames(data)
        colnames(sigma_mat) <- colnames(data)

    } else{
      stop("Type not recognized. Please, enter either cov or cor.")
    }


    return(list(sigma_mat, rho))

  }


CVglasso_mod = function(X = NULL,
                        lam = NULL,
                        nlam = 50,
                        seed = 1234,
                        path = FALSE,
                        tol = 1e-04,
                        maxit = 10000,
                        adjmaxit = NULL,
                        K = 5,
                        crit.cv = c("loglik", "AIC", "BIC"),
                        start = c("warm", "cold"),
                        cores = 1,
                        trace = c("progress", "print", "none")) {
  if (is.null(X)) {
    stop("Must provide entry for X!")
  }
  if (!all(lam > 0)) {
    stop("lam must be positive!")
  }
  if (!(all(c(tol, maxit, adjmaxit, K, cores, nlam) > 0))) {
    stop("Entry must be positive!")
  }
  if (!(all(sapply(
    c(tol, maxit, adjmaxit, K, cores, nlam),
    length
  ) <= 1))) {
    stop("Entry must be single value!")
  }
  if (all(c(maxit, adjmaxit, K, cores) %% 1 != 0)) {
    stop("Entry must be an integer!")
  }
  if (cores < 1) {
    stop("Number of cores must be positive!")
  }
  if (cores > 1 && path) {
    cat("\nParallelization not possible when producing solution path. Setting cores = 1...")
    cores = 1
  }

  K = ifelse(path, 1, K)

  if (cores > K) {
    cat("\nNumber of cores exceeds K... setting cores = K")
    cores = K
  }
  if (is.null(adjmaxit)) {
    adjmaxit = maxit
  }

  # match values
  crit.cv = match.arg(crit.cv)
  start = match.arg(start)
  trace = match.arg(trace)
  call = match.call()
  MIN.error = AVG.error = CV.error = NULL
  n = nrow(X)

  # compute sample correlation matrix
  C = stats::cor(X)
  Cminus = C
  diag(Cminus) = 0

  # compute grid of lam values, if necessary
  if (is.null(lam)) {
    lam = sort(exp(seq(log(1 / nlam), log(1), length.out = nlam)))

  } else {
    # sort lambda values
    lam = sort(lam)

  }

  # perform cross validation, if necessary
  if ((length(lam) > 1) & (!is.null(X) || path)) {
    # run CV in parallel?
    if (cores > 1) {
      # execute CVP
      GLASSO = CVP_mod(
        X = X,
        lam = lam,
        tol = tol,
        maxit = maxit,
        adjmaxit = adjmaxit,
        K = K,
        crit.cv = crit.cv,
        start = start,
        cores = cores,
        trace = trace
      )
      MIN.error = GLASSO$min.error
      AVG.error = GLASSO$avg.error
      CV.error = GLASSO$cv.error

    } else {
      GLASSO = CV_mod(
        X = X,
        lam = lam,
        path = path,
        tol = tol,
        maxit = maxit,
        adjmaxit = adjmaxit,
        K = K,
        crit.cv = crit.cv,
        start = start,
        trace = trace
      )
      MIN.error = GLASSO$min.error
      AVG.error = GLASSO$avg.error
      CV.error = GLASSO$cv.error
      Path = GLASSO$path

    }

    # print warning if lam on boundary
    if ((GLASSO$lam == lam[1]) && (length(lam) != 1) && !path) {
      cat(
        "\nOptimal tuning parameter on boundary... consider providing a smaller lam value or decreasing lam.min.ratio!"
      )
    }

    # provide estimate that is pd and dual feasible
    alpha = min(c(GLASSO$lam / max(abs(Cminus)), 1))
    init = (1 - alpha) * C
    diag(init) = diag(C)

    # compute final estimate at best tuning parameters
    lam_ = GLASSO$lam
    GLASSO = glasso::glasso(
      s = C,
      rho = lam_,
      thr = tol,
      maxit = maxit,
      penalize.diagonal = FALSE,
      start = "warm",
      w.init = init,
      wi.init = diag(ncol(C)),
      trace = FALSE
    )
    GLASSO$lam = lam_

  } else {
    # execute ADMM_sigmac
    if (length(lam) > 1) {
      stop("Must set specify X, set path = TRUE, or provide single value for lam.")
    }

    # provide estimate that is pd and dual feasible
    alpha = min(c(lam / max(abs(Cminus)), 1))
    init = (1 - alpha) * C
    diag(init) = diag(C)

    GLASSO = glasso::glasso(
      s = C,
      rho = lam,
      thr = tol,
      maxit = maxit,
      penalize.diagonal = FALSE,
      start = "warm",
      w.init = init,
      wi.init = diag(ncol(C)),
      trace = FALSE
    )
    GLASSO$lam = lam

  }

  diag_pen = 1 - diag(ncol(C))

  # compute penalized loglik
  loglik = sum(GLASSO$wi * C) - determinant(GLASSO$wi, logarithm = TRUE)$modulus[1] + GLASSO$lam * sum(abs(diag_pen * GLASSO$wi))

  # return values
  tuning = matrix(c(log10(GLASSO$lam), GLASSO$lam), ncol = 2)
  colnames(tuning) = c("log10(lam)", "lam")
  if (!path) {
    Path = NULL
  }

  returns = list(
    Call = call,
    Iterations = GLASSO$niter,
    Tuning = tuning,
    Lambdas = lam,
    maxit = maxit,
    Omega = GLASSO$wi,
    Sigma = GLASSO$w,
    Path = Path,
    Loglik = loglik,
    MIN.error = MIN.error,
    AVG.error = AVG.error,
    CV.error = CV.error
  )

  class(returns) = "CVglasso"
  return(returns)

}

CV_mod = function(X = NULL,
                  lam = exp(seq(log(1 / 100), log(1), length.out = 100)),
                  seed = 1234,
                  path = FALSE,
                  tol = 1e-04,
                  maxit = 10000,
                  adjmaxit = NULL,
                  K = 5,
                  crit.cv = c("loglik", "AIC", "BIC"),
                  start = c("warm", "cold"),
                  cores = 1,
                  trace = c("progress", "print", "none")) {
  set.seed(seed)

  # match values
  crit.cv = match.arg(crit.cv)
  start = match.arg(start)
  trace = match.arg(trace)
  lam = sort(lam)

  # initialize
  Path = NULL
  initmaxit = maxit

  # compute sample correlation matrix
  C = stats::cor(X)
  C.train = C.valid = C

  CV_errors = array(0, c(length(lam), K))

  # set progress bar
  if (trace == "progress") {
    progress = utils::txtProgressBar(max = K * length(lam), style = 3)
  }

  # no need to create folds if K = 1
  if (K == 1) {
    # initialize Path, if necessary
    if (path) {
      Path = array(0, c(ncol(C), ncol(C), length(lam)))
    }

  } else {
    # designate folds and shuffle -- ensures randomized folds
    n = nrow(X)
    ind = sample(n)

  }

  # parse data into folds and perform CV
  for (k in 1:K) {
    if (K > 1) {
      # training set
      leave.out = ind[(1 + floor((k - 1) * n / K)):floor(k *
                                                           n / K)]
      X.train = X[-leave.out, , drop = FALSE]
      X_bar = apply(X.train, 2, mean)
      X.train = scale(X.train, center = X_bar, scale = FALSE)

      # validation set
      X.valid = X[leave.out, , drop = FALSE]
      X.valid = scale(X.valid, center = X_bar, scale = FALSE)

      C.train = stats::cor(X.train)
      C.valid = stats::cor(X.valid)

    }


    # re-initialize values for each fold
    maxit = initmaxit
    init = C.train
    initOmega = diag(ncol(C.train))

    # provide estimate that is pd and dual feasible
    Cminus = C.train
    diag(Cminus) = 0
    alpha = min(c(lam[1] / max(abs(Cminus)), 1))
    init = (1 - alpha) * C.train
    diag(init) = diag(C.train)


    # loop over all tuning parameters
    for (i in 1:length(lam)) {
      # set temporary tuning parameter
      lam_ = lam[i]

      # compute the penalized likelihood precision matrix estimator
      GLASSO = glasso::glasso(
        s = C.train,
        rho = lam_,
        thr = tol,
        maxit = maxit,
        penalize.diagonal = FALSE,
        start = "warm",
        w.init = init,
        wi.init = initOmega,
        trace = FALSE
      )

      if (start == "warm") {
        # option to save initial values for warm starts
        init = GLASSO$w
        initOmega = GLASSO$wi
        maxit = adjmaxit

      }

      CV_errors[i, k] = sum(GLASSO$wi * C.valid) - determinant(GLASSO$wi, logarithm = TRUE)$modulus[1]

      # update for crit.cv, if necessary
      if (crit.cv == "AIC") {
        CV_errors[i, k] = CV_errors[i, k] + sum(GLASSO$wi !=
                                                  0)
      }
      if (crit.cv == "BIC") {
        CV_errors[i, k] = CV_errors[i, k] + sum(GLASSO$wi !=
                                                  0) * log(n) / 2
      }

      # save estimate if path = TRUE
      if (path) {
        Path[, , i] = GLASSO$wi
      }

      # update progress bar
      if (trace == "progress") {
        utils::setTxtProgressBar(progress, i + (k - 1) * length(lam))

        # if not quiet, then print progress lambda
      } else if (trace == "print") {
        cat("\nFinished lam = ", paste(lam_, sep = ""))
      }
    }

    # if not quiet, then print progress kfold
    if (trace == "print") {
      cat("\nFinished fold ", paste(k, sep = ""))
    }
  }

  # determine optimal tuning parameters
  AVG = apply(CV_errors, 1, mean)
  best_lam = lam[which.min(AVG)]
  error = min(AVG)


  # return best lam and alpha values
  return(
    list(
      lam = best_lam,
      path = Path,
      min.error = error,
      avg.error = AVG,
      cv.error = CV_errors
    )
  )

}


CVP_mod = function(X = NULL,
                   lam = exp(seq(log(1 / 100), log(1), length.out = 100)),
                   seed = 1234,
                   tol = 1e-04,
                   maxit = 10000,
                   adjmaxit = NULL,
                   K = 5,
                   crit.cv = c("loglik", "AIC", "BIC"),
                   start = c("warm", "cold"),
                   cores = 1,
                   trace = c("progress", "print", "none")) {
  set.seed(seed)

  # match values
  crit.cv = match.arg(crit.cv)
  start = match.arg(start)
  trace = match.arg(trace)
  lam = sort(lam)

  # make cluster and register cluster
  num_cores = parallel:::detectCores()
  if (cores > num_cores) {
    cat("\nOnly detected", paste(num_cores, "cores...", sep = " "))
  }
  if (cores > K) {
    cat("\nNumber of cores exceeds K... setting cores = K")
    cores = K
  }

  cluster = parallel::makeCluster(cores, type = "FORK")
  doParallel::registerDoParallel(cluster)

  # use cluster for each fold in CV
  n = nrow(X)
  ind = sample(n)
  k = NULL
  CV = foreach::foreach(
    k = 1:K,
    .packages = "glasso",
    .combine = "cbind",
    .inorder = FALSE
  ) %dopar% {
    # set progress bar
    if (trace == "progress") {
      progress = utils::txtProgressBar(max = length(lam), style = 3)
    }

    # training set
    leave.out = ind[(1 + floor((k - 1) * n / K)):floor(k * n / K)]
    X.train = X[-leave.out, , drop = FALSE]
    X_bar = apply(X.train, 2, mean)
    X.train = scale(X.train, center = X_bar, scale = FALSE)
    n.train = nrow(X.train)

    # validation set
    X.valid = X[leave.out, , drop = FALSE]
    X.valid = scale(X.valid, center = X_bar, scale = FALSE)
    n.valid = nrow(X.valid)

    C.train = stats::cor(X.train)
    C.valid = stats::cor(X.valid)

    # initial estimates
    init = C.train
    initOmega = diag(ncol(C.train))
    CV_error = array(0, length(lam))


    # provide estimate that is pd and dual feasible
    Cminus = C.train
    diag(Cminus) = 0
    alpha = min(c(lam[1] / max(abs(Cminus)), 1))
    init = (1 - alpha) * C.train
    diag(init) = diag(C.train)

    # loop over all tuning parameters
    for (i in 1:length(lam)) {
      # set temporary tuning parameter
      lam_ = lam[i]

      # compute the penalized likelihood precision matrix estimator
      GLASSO = glasso::glasso(
        s = C.train,
        rho = lam_,
        thr = tol,
        maxit = maxit,
        penalize.diagonal = FALSE,
        start = "warm",
        w.init = init,
        wi.init = initOmega,
        trace = FALSE
      )

      if (start == "warm") {
        # option to save initial values for warm starts
        init = GLASSO$w
        initOmega = GLASSO$wi
        maxit = adjmaxit

      }

      # compute the observed negative validation loglikelihood
      CV_error[i] = sum(GLASSO$wi * C.valid) - determinant(GLASSO$wi, logarithm = TRUE)$modulus[1]

      # update for crit.cv, if necessary
      if (crit.cv == "AIC") {
        CV_error[i] = CV_error[i] + sum(GLASSO$wi != 0)
      }
      if (crit.cv == "BIC") {
        CV_error[i] = CV_error[i] + sum(GLASSO$wi != 0) *
          log(n) / 2
      }

      # update progress bar
      if (trace == "progress") {
        utils::setTxtProgressBar(progress, i)

        # if not quiet, then print progress lambda
      } else if (trace == "print") {
        cat("\nFinished lam = ", paste(lam_, sep = ""))
      }
    }

    # return CV errors
    return(CV_error)

  }

  # determine optimal tuning parameters
  AVG = apply(CV, 1, mean)
  best_lam = lam[which.min(AVG)]
  error = min(AVG)

  # stop cluster
  parallel::stopCluster(cluster)

  # return best lam and alpha values
  return(list(
    lam = best_lam,
    min.error = error,
    avg.error = AVG,
    cv.error = CV
  ))

}
