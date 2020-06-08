#' t-Lasso (TLASSO) Covariance Estimation
#'
#' Computes the TLASSO estimator of the covariance matrix.
#'
#' @param data an nxp data matrix.
#' @param rho a double, the non-negative regularization parameter \eqn{\rho} for lasso. \eqn{\rho=0} means no regularization.
#' Can be a scalar (usual) or a symmetric p by p matrix, or a vector of length p. In the latter case, the penalty matrix has jkth element \eqn{\sqrt{\rho_j*\rho_k}}.
#' Default value is 0.
#' @param df an integer, indicating the degrees of freedom of the assumed t-distribution. Default value is 3.
#' @param pendiag_log a logical, indicating whether the diagonal of the sample covariance matrix is to be penalized (TRUE) or not (FALSE). Default value is FALSE.
#' @param tol a double, indicating the tolerance for the glasso algorithm. Default value is set to 1e-05.
#' @param maxit an integer, indicating the maximum number of iterations for the glasso algorithm. Default value is set to 10000.
#' @param symmetric_log a logical, indicating whether the output should be a symmetric matrix (TRUE) or not necessarily (FALSE). Default value is set to TRUE.
#' @param theta_init a pxp initial matrix for the inverse of the covariance matrix. Default value is NULL and the sample inverse for the t-distribution is used.
#'
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the lasso penalty.
#' }
#'
#' @details The TLASSO estimator is elaborated in detail in \insertCite{finegold2011robust;textual}{CovEstim}.
#' Originally developed by \insertCite{torri2019sparse;textual}{CovEstim}.
#'
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_tlasso <- sigma_estim_tlasso(sp_rets, rho=0.001)[[1]]
#'
#'
#' @import glasso
#' @importFrom Rdpack reprompt
#' @references
#'\insertAllCited
#'
#' @export sigma_estim_tlasso
#'
sigma_estim_tlasso <-
  function(data,
           rho = NULL,
           pendiag_log = FALSE,
           df = 3,
           tol = 1e-5,
           maxit = 10000,
           symmetric_log = TRUE,
           theta_init = NULL) {
    data <- as.matrix(data)
    n <- dim(data)[1]
    p <- dim(data)[2]
    mu <- apply(data, 2, mean)
    centered <- apply(data, 2, function(x)
      x - mean(x))
    sigma_sample <- t(centered) %*% centered / (n - 1)
    tau <- rep(0, n)

    if (is.null(theta_init)) {
      theta_mat <- solve(sigma_sample * ((df - 2) / df))
    } else {
      theta_mat <- theta_init
    }

    # loop E and M steps
    i <- 1
    conv_hist <- c()
    conv_hist[1] <- 1000

    while (conv_hist[i] > tol && i < maxit) {
      #E step
      tau <-
        apply(data, 1, function(x) {
          (df + p) / (t(x - mu) %*% theta_mat %*% (x - mu) + df)
        })

      #M step
      mu <- (t(data) %*% tau) / sum(tau)
      mu_large <- t(matrix(rep(mu, n), p, n))
      tau_large <- matrix(rep(tau, p), n, p)
      sigma_mat <-
        t(tau_large * (data - mu_large)) %*% (data - mu_large) / n
      theta_old <- theta_mat

      glasso_results <-
        glasso::glasso(sigma_mat, rho, penalize.diagonal = pendiag_log)

      # force symmetry (required sometimes for numerical issues)
      if (symmetric_log == TRUE) {
        theta_mat <- 0.5 * (glasso_results$wi + t(glasso_results$wi))
        sigma_mat <- solve(theta_mat)
      } else {
        theta_mat <- glasso_results$wi
        sigma_mat <- glasso_results$w
      }

      conv_hist[i + 1] <- max(abs(theta_mat - theta_old))

      i <- i + 1
    }

    return(list(sigma_mat, rho))
  }
