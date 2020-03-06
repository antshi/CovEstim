#' Graphical Lasso (GLASSO) Covariance Estimation
#'
#' Computes the GLASSO estimator of the covariance matrix.
#'
#' @param data an nxp data matrix.
#' @param rho a double, the non-negative regularization parameter \eqn{\rho} for lasso. \eqn{\rho=0} means no regularization.
#' Can be a scalar (usual) or a symmetric p by p matrix, or a vector of length p. In the latter case, the penalty matrix has jkth element \eqn{\sqrt{\rho_j*\rho_k}}. Default value is NULL
#' and an optimal regularization is computed with the cross-validation (CV) procedure as in \insertCite{cvglassopackage;textual}{CovEstim}.
#' @param rho_seq a user-sapplied sequence of possible values for the regularization parameter \eqn{\rho}. Default value is NULL and
#' a sequence is defined as in \insertCite{cvglassopackage;textual}{CovEstim}.
#' @param rho_seq_len an integer, indicating the length for setting the rho sequence.
#' @param rho_seq_min_ratio a double, the minimum ratio for deviation between the rho values in the sequence.
#' @param nfolds an integer, indicating the number of folds for the CV. Default value is 5.
#' @param crit a character, indicating which selection criterion within the CV. Possible values are "loglik", "AIC" and "BIC". Default is set to "BIC".
#' @param pendiag_log a logical, indicating whether the diagonal of the sample covariance matrix is to be penalized (TRUE) or not (FALSE). Default value is FALSE.
#' @param start a character, specifying the start type of the glasso algorithm. Possible values are "warm" or "cold". Default value is "cold".
#' @param tol a double, indicating the tolerance for the glasso algorithm. Default value is set to 1e-05.
#' @param maxit an integer, indicating the maximum number of iterations for the glasso algorithm. Default value is set to 10000.
#' @param cores an integer, indicating how many cores should be used for the CV. Default value is 1. cores cannot be higher than the maximum number of cores of the processor in use.
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
#' sigma_glasso <- sigma_estim_glasso(sp_rets, rho=0.0001)[[1]]
#'
#' @import glasso
#' @import CVglasso
#' @importFrom Rdpack reprompt
#' @references
#'\insertAllCited
#'
#' @export sigma_estim_glasso
#'
sigma_estim_glasso <-
  function(data,
           rho = NULL,
           rho_seq = NULL,
           rho_seq_len = 10,
           rho_seq_min_ratio = 0.001,
           nfolds = 5,
           crit = "BIC",
           pendiag_log = FALSE,
           start = "cold",
           tol = 1e-05,
           maxit = 10000,
           cores = 1) {
    data <- as.matrix(data)
    n <- dim(data)[1]
    centered <- apply(data, 2, function(x)
      x - mean(x))
    sigma_sample <- t(centered) %*% centered / (n - 1)

    if (is.null(rho)) {
      if (is.null(rho_seq)) {
        sigma_sample_diag0 <- sigma_sample
        diag(sigma_sample_diag0) <- 0
        rho_max <- max(abs(sigma_sample_diag0))
        rho_min <- rho_seq_min_ratio * rho_max
        rho_seq <-
          10 ^ seq(log10(rho_min), log10(rho_max), length = rho_seq_len)
      }
      rho <-
        CVglasso::CVglasso(
          X = data,
          lam = rho_seq,
          K = nfolds,
          crit.cv = crit,
          diagonal = pendiag_log,
          start = start,
          tol = tol,
          maxit = maxit,
          cores = cores,
          trace = "none"
        )$Tuning[2]
    }
    sigma_mat <-
      glasso::glasso(
        sigma_sample,
        rho = rho,
        penalize.diagonal = pendiag_log,
        start = start,
        thr = tol,
        maxit = maxit
      )$w

    rownames(sigma_mat) <- colnames(data)
    colnames(sigma_mat) <- colnames(data)

    return(list(sigma_mat, rho))

  }
