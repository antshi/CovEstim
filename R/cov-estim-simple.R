#' Sample Covariance Estimation
#'
#' Computes the sample estimator of the covariance matrix.
#'
#' @param data an nxp data matrix.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here an NA.
#' }
#'
#' @details The sample estimator of the covariance matrix for a data matrix X is computed with the following formula:
#' \deqn{\hat{\Sigma}=\frac{1}{n-1} \left(X - \widehat{{\mu}} {1} \right)' \left({X} -  \widehat{{\mu}}{1}\right)}
#' where \eqn{\mu=\bar{x}_{j}=\frac{1}{n}\sum_{i=1}^{n}x_{ij}} (for \eqn{i=1,\ldots, n} and \eqn{j=1,\ldots,p}) is the sample mean vector and \eqn{1} is an 1xp vector of ones.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_sample <- sigma_estim_sample(sp_rets)[[1]]
#'
#' @export sigma_estim_sample
#'
sigma_estim_sample <- function(data) {
  data <- as.matrix(data)
  names_data <- colnames(data)
  n <- dim(data)[1]
  centered <- apply(data, 2, function(x)
    x - mean(x))
  sigma_mat <- t(centered) %*% centered / (n - 1)

  rownames(sigma_mat) <- names_data
  colnames(sigma_mat) <- names_data

  return(list(sigma_mat, NA))
}


#' Maximum-Likelihood Covariance Estimation
#'
#' Computes the Maximum-Likelihood estimator of the covariance matrix.
#'
#' @param data an nxp data matrix.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here an NA.
#' }
#'
#' @details The Maximum-Likelihood estimator of the covariance matrix for a data matrix X is computed with the following formula:
#' \deqn{\hat{\Sigma}=\frac{1}{n} \left(X - \widehat{{\mu}} {1} \right)' \left({X} -  \widehat{{\mu}}{1}\right)}
#' where \eqn{\mu=\bar{x}_{j}=\frac{1}{n}\sum_{i=1}^{n}x_{ij}} for (for \eqn{i=1,\ldots, n} and \eqn{j=1,\ldots,p}) is the sample mean vector and \eqn{1} is an 1xp vector of ones.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_ml <- sigma_estim_ml(sp_rets)[[1]]
#'
#' @export sigma_estim_ml
#'
sigma_estim_ml <- function(data) {
  data <- as.matrix(data)
  names_data <- colnames(data)
  n <- dim(data)[1]
  centered <- apply(data, 2, function(x)
    x - mean(x))
  sigma_mat <- t(centered) %*% centered / n

  rownames(sigma_mat) <- names_data
  colnames(sigma_mat) <- names_data

  return(list(sigma_mat, NA))
}


#' Bayes-Stein Covariance Estimation
#'
#' Computes the Bayes-Stein estimator of the covariance matrix.
#'
#' @param data an nxp data matrix.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here an NA.
#' }
#'
#' @details The Bayes-Stein estimator of the covariance matrix is computed according to \insertCite{jorion1986bayes;textual}{CovEstim}.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_bs <- sigma_estim_bs(sp_rets)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#'\insertAllCited
#'
#' @export sigma_estim_bs
#'
sigma_estim_bs <- function(data) {
  data <- as.matrix(data)
  names_data <- colnames(data)
  n <- dim(data)[1]
  p <- dim(data)[2]
  mu <- colMeans(data)
  centered <- apply(data, 2, function(x)
    x - mean(x))
  rm(data)
  gc()
  sigma_ml <- t(centered) %*% centered / n
  sigma_ml_inv <- solve(sigma_ml)
  sigma <- sigma_ml * (n / (n - p - 2))
  sigma_inv <- solve(sigma)

  ones <- rep.int(1, p)
  mug <-
    as.numeric(ones %*% sigma_ml_inv %*% mu / as.numeric(ones %*% sigma_ml_inv %*%
                                                           ones))
  lambda <-
    as.numeric((p + 2) / as.numeric(t(mu - mug * ones) %*% sigma_inv %*% (mu - mug *
                                                                            ones)))

  sigma_mat <-
    (1 + 1 / (n + lambda)) * sigma + (lambda / (n * (n + 1 + lambda))) * (ones %*%
                                                                            t(ones)) / as.numeric(ones %*% sigma_inv %*% ones)

  rownames(sigma_mat) <- names_data
  colnames(sigma_mat) <- names_data

  return(list(sigma_mat, NA))
}
