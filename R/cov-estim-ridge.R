#' Ridge-Penalized Covariance Estimation
#'
#' Computes the covariance matrix estimator with ridge-penalty.
#'
#' @param data an nxp data matrix.
#' @param rho a double, indicating the ridge penalty.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the user-supplied ridge penalty \eqn{\rho}.
#' }
#'
#' @details The ridge-penalized covariance matrix for a data matrix X is computed with the following formula:
#' \deqn{\hat{\Sigma}=V\left(\rho C+(1-\rho)I\right)V,}
#' where V is the sample volatility matrix, I is a pxp identity matrix,
#' C is the sample correlation matrix and \eqn{\rho} is the user-sapplied ridge penalty parameter.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_ridge <- sigma_estim_ridge(sp_rets, rho=0.01)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#'\insertRef{warton2008penalized}{CovEstim}
#'
#' @export sigma_estim_ridge
#'
sigma_estim_ridge <- function(data, rho) {
  data <- as.matrix(data)
  p <- dim(data)[2]
  corr_mat <- stats::cor(data)
  vola_mat <- diag(apply(data, 2, stats::sd, na.rm = TRUE))
  ident_mat <- diag(1, p)
  corr_mat_ridge <- rho * corr_mat + (1 - rho) * ident_mat

  sigma_mat <- as.matrix(vola_mat %*% corr_mat_ridge %*% vola_mat)
  rownames(sigma_mat) <- colnames(data)
  colnames(sigma_mat) <- colnames(data)

  return(list(sigma_mat, rho))
}
