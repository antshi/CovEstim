#' Condition-Number Regularized Covariance Estimation
#'
#' Computes the condition-number regularized (CONDREG) estimator of the covariance matrix.
#'
#' @param data an nxp data matrix.
#' @param k a double, indicating the regularization parameter (or the maximum condition number for the estimated covariance matrix).
#' Default value is NULL and the optimal k is found with a cross-validation (CV), optimizing the negative likelihood.
#' @param k_seq a vector of doubles, specifying the grid of k values to search over within the CV.
#' Default value is NULL and the sequence is generated in dependence of the sample covariance matrix, the user-supplied length and the minimum deviation ratio of the values.
#' @param k_seq_len an integer, indicating the length of k_seq. Default value is 50.
#' @param nfolds an integer, specifying the number of folds for the CV. Default value is 5.
#' @param zeromean_log a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE). Default value is FALSE.
#'
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the regularization parameter k.
#' }
#'
#' @details The CONDREG estimator is elaborated in detail in \insertCite{won2013condition;textual}{CovEstim}. More information on the functionality can be found in
#' \insertCite{condregpackage;textual}{CovEstim}.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' results_condreg <- sigma_estim_condreg(sp_rets)
#' sigma_condreg <- results_condreg[[1]]
#' param_condreg <- results_condreg[[2]]
#' sigma_condreg <- sigma_estim_condreg(sp_rets, k=100)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#'\insertAllCited
#'
#' @export sigma_estim_condreg
#'
sigma_estim_condreg <-
  function(data,
           k = NULL,
           k_seq = NULL,
           k_seq_len = 50,
           nfolds = 5,
           zeromean_log = FALSE) {
    data <- as.matrix(data)
    if (!zeromean_log) {
      centered <- apply(data, 2, function(x)
        x - mean(x))
      n <- dim(data)[1] - 1
    } else{
      centered <- data
      n <- dim(data)[1]
    }

    sigma_sample <- t(centered) %*% centered / n

    svd_sigma_sample <- svd(sigma_sample)
    eigenval <- svd_sigma_sample$d
    u <- svd_sigma_sample$u

    if (is.null(k)) {
      k <- kmax_cv_condreg(data, k_seq, k_seq_len, nfolds, zeromean_log)
    }

    sol <- lbar_solver(eigenval, k)
    eigenval_condreg <- as.numeric(sol$lbar)
    sigma_mat <- u %*% diag(eigenval_condreg) %*% t(u)

    colnames(sigma_mat) <- colnames(data)
    rownames(sigma_mat) <- colnames(data)

    return(list(sigma_mat, k))

  }

kseq_calc <- function(mat,
                      k_seq_len = 50) {
  k_max <- kappa(mat, exact = TRUE)
  k_seq <-
    sort(10 ^ seq(log10(1), log10(k_max), length = k_seq_len), decreasing = TRUE)
  return(k_seq)
}


kseq_calc_orig <- function(k_max = 30,
                           k_seq_len = 50) {
  k_seq <- 1 / (seq(1, k_max, length.out = k_seq_len))
  k_seq <- k_seq - min(k_seq)
  k_seq <-  k_seq / max(k_seq)
  k_seq <- (k_seq * (k_max - 1)) + 1
  return(k_seq)
}

kmax_cv_condreg <-
  function(data,
           k_seq = NULL,
           k_seq_len = 50,
           nfolds = 5,
           zeromean_log = FALSE) {
    data <- as.matrix(data)
    if (!zeromean_log) {
      centered <- apply(data, 2, function(x)
        x - mean(x))
      n <- dim(data)[1] - 1
    } else{
      centered <- data
      n <- dim(data)[1]
    }

    p <- dim(data)[2]
    sigma_sample <- t(centered) %*% centered / n

    if (is.null(k_seq)) {
      k_seq <- kseq_calc(sigma_sample, k_seq_len)
    }
    g <- length(k_seq)
    sections <- cut(1:n, breaks = nfolds, labels = c(1:nfolds))
    condmax <- 1
    negloglikelihood <- matrix(0, nfolds, g)

    for (i in 1:nfolds) {
      tsindx <- which(sections == i)
      trindx <- which(sections != i)
      xtrain <- data[trindx, , drop = FALSE]
      xtest <- data[tsindx, , drop = FALSE]
      ntr <- nrow(xtrain)
      nts <- nrow(xtest)

      soln <- condreg_bulk(xtrain, k_seq, zeromean_log)
      ytest <- xtest %*% soln$Q
      a <- array(ytest ^ 2, dim = c(nts, p, g))
      a <- aperm(a, c(2, 1, 3))
      b <- array(1 / soln$lbar, dim = c(g, p, nts))
      b <- aperm(b, c(2, 3, 1))
      ztest <- a * b
      negloglikelihood[i, ] <-
        rowSums(aperm(ztest, dim = c(3, 1, 2))) / nts + rowSums(log(soln$lbar))

      L <- rep(0, p)
      L[1:min(ntr, p)] <- soln$L[1:min(ntr, p)]
      condmax <- max(condmax, L[1] / L[min(ntr, p)])

    }

    nL <- colSums(negloglikelihood)
    minind <- floor(stats::median(which(nL == min(nL))))
    kmaxopt <- min(k_seq[minind], condmax)

    return(kmaxopt)

  }

lbar_solver <- function(L, k, dir = "forward") {

  p <- length(L)
  g <- length(k)

  lbar <- matrix(0, g, p)
  uopt <- rep(0, g)
  intv <- rep(0, g)

  L[L < .Machine$double.eps] <- .Machine$double.eps

  degenindx <- (k > (L[1] / L[p]))

  if (sum(degenindx) > 0) {
    lbar[which(degenindx), ] <-
      matrix(L, sum(degenindx), p, byrow = TRUE)
    uopt[degenindx] <- pmax(1 / k[degenindx] / L[p], 1 / L[1])
    intv[degenindx] <- TRUE
  }

  if (any(!degenindx)) {
    kmax1 <- k[!degenindx]

    if (dir == "forward") {
      path <- path_forward(L)
    } else if (dir == "backward") {
      path <- path_backward(L)
    } else {
      stop("dir is either 'forward' or 'backward'\n")
    }

    tmp <-
      stats::approx(path$k, 1 / path$u, kmax1) #linear interpolation
    uopt <- 1 / tmp$y

    lambda <-
      pmin(matrix(rep(kmax1 * uopt, p), ncol = p), pmax(matrix(rep(uopt, p), ncol =
                                                                 p), matrix(
                                                                   rep(1 / L, length(kmax1)), ncol = p, byrow = TRUE
                                                                 )))

    lbar[!degenindx, ] <- 1 / lambda
    uopt[!degenindx] <- uopt
    intv[!degenindx] <- FALSE
  }

  return(list(
    lbar = lbar,
    uopt = uopt,
    intv = intv
  ))

}


path_forward <- function(L) {

  p <- length(L)
  idxzero <- L < .Machine$double.eps
  numzero <- sum(idxzero)
  L[idxzero] <- .Machine$double.eps

  u_cur <- 1 / mean(L)
  v_cur <- u_cur

  alpha <- 0
  while (u_cur > 1 / L[alpha + 1]) {
    alpha <- alpha + 1
  }
  beta <- alpha + 1
  slope_num <- sum(L[1:alpha])
  slope_denom <- sum(L[beta:p])

  u <- u_cur
  v <- v_cur
  kmax <- 1

  r <- p - numzero

  while (alpha >= 1 && beta <= r) {
    # intersection of the line passing (u_cur,v_cur) of slope 'slope' with
    # the rectangle [ 1/d(alpha), 1/d(alpha+1) ] x [ 1/d(beta-1), 1/d(beta) ]
    h_top <- 1 / L[beta]
    v_left <- 1 / L[alpha]

    # First, compute intersection with the horizontal line v=1/d(beta)
    v_new <- h_top
    u_new <- u_cur - slope_denom * (v_new - v_cur) / slope_num

    # If u_new is outside the rectangle,
    # compute intersection with the vertical line u=1/d(alpha+1).
    if (u_new < v_left) {
      u_new <- v_left
      v_new <- v_cur - slope_num * (u_new - u_cur) / slope_denom
    }

    # update
    if (abs(u_new - v_left) < .Machine$double.eps) {
      # keep this order!
      slope_num <- slope_num - L[alpha]   # first
      alpha <- alpha - 1                  # second
    }
    if (abs(v_new - h_top) < .Machine$double.eps) {
      # keep this order!
      slope_denom <- slope_denom - L[beta]  # first
      beta <- beta + 1                      # second
    }

    new_kmax <- v_new / u_new
    u <- c(u, u_new)
    v <- c(v, v_new)
    kmax <- c(kmax, new_kmax)

    u_cur <- u_new
    v_cur <- v_new
  }

  kmax <- c(kmax, Inf)
  u <- c(u, u_new)
  v <- c(v, Inf)

  return(list(k = kmax, u = u, v = v))
}


path_backward <- function(L) {
  p <- length(L)

  idxzero <- L < .Machine$double.eps
  numzero <- sum(idxzero)
  L[idxzero] <- .Machine$double.eps

  r <- p - numzero   # rank

  # ending point finding algorithm
  alpha <- 1
  slope_num <- sum(L[1:alpha])
  u_cur <- (alpha + p - r) / slope_num
  while (u_cur < 1 / L[alpha] || u_cur > 1 / L[alpha + 1]) {
    alpha <- alpha + 1
    slope_num <- slope_num + L[alpha]
    u_cur <- (alpha + p - r) / slope_num
  }

  v_cur <- 1 / L[r]

  beta <- r
  slope_denom <- sum(L[beta:p])

  # vertical half-infinite line segment
  u <- c(u_cur, u_cur)
  v <- c(v_cur, Inf)
  kmax <- c(v_cur / u_cur, Inf)

  is_done <- FALSE
  while (!is_done) {
    # intersection of the line passing (u_cur,v_cur) of slope 'slope' with
    # the rectangle [ 1/d(alpha), 1/d(alpha+1) ] x [ 1/d(beta-1), 1/d(beta) ]
    h_bottom <- 1 / L[beta - 1]
    v_right <- 1 / L[alpha + 1]

    # First, check the intersection with the diagonal line v=u.
    u_new <-
      (slope_num * u_cur + slope_denom * v_cur) / (slope_num + slope_denom)
    v_new <- u_new
    if (u_new < v_right && v_new > h_bottom) {
      is_done <- TRUE
      u <- c(u_new, u)
      v <- c(v_new, v)
      kmax <- c(1, kmax)
      break
    }

    # Compute intersection with the horizontal line v=1/d(beta-1)
    v_new <- h_bottom
    u_new <- u_cur - slope_denom * (v_new - v_cur) / slope_num

    # If u_new is outside the rectangle,
    # compute intersection with the vertical line u=1/d(alpha+1).
    if (u_new > v_right) {
      u_new <- v_right

      v_new <- v_cur - slope_num * (u_new - u_cur) / slope_denom
    }

    # update
    if (abs(u_new - v_right) < .Machine$double.eps) {
      # keep this order!
      alpha <- alpha + 1                   # first
      slope_num <- slope_num + L[alpha]    # second
    }
    if (abs(v_new - h_bottom) < .Machine$double.eps) {
      # keep this order!
      beta <- beta - 1                      # first
      slope_denom <- slope_denom + L[beta]  # second
    }

    new_kmax <- v_new / u_new

    u <- c(u_new, u)
    v <- c(v_new, v)
    kmax <- c(new_kmax, kmax)

    u_cur <- u_new
    v_cur <- v_new

  }

  return(list(k = kmax, u = u, v = v))
}


condreg_bulk <- function(data, k_seq, zeromean_log = FALSE) {
  data <- as.matrix(data)
  if (!zeromean_log) {
    centered <- apply(data, 2, function(x)
      x - mean(x))
    n <- dim(data)[1] - 1
  } else{
    centered <- data
    n <- dim(data)[1]
  }

  sigma_sample <- t(centered) %*% centered / n
  p <- dim(data)[2]

  g <- length(k_seq)
  svdS <- svd(sigma_sample)
  soln <- lbar_solver(svdS$d, k_seq)

  ## if n<p there are some 0 eigenvalues
  if (n < p) {
    soln$lbar <- cbind(soln$lbar, matrix(0, g, max(p - n, 0)))
  }

  list(Q = svdS$u,
       lbar = soln$lbar,
       L = svdS$d)
}
