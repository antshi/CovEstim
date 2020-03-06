#' Stein-Haff Covariance Estimation
#'
#' Computes the Stein-Haff estimator of the covariance matrix.
#'
#' @param data an nxp data matrix.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here an NA.
#' }
#' @details The estimation procedure follows \insertCite{stein1975estimation;textual}{CovEstim} and \insertCite{haff1991variational;textual}{CovEstim}.
#' Originally found under \insertCite{stcovpackage;textual}{CovEstim}.
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_sh <- sigma_estim_sh(sp_rets)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#'\insertAllCited
#'
#' @export sigma_estim_sh
#'
sigma_estim_sh <- function(data) {
  data <- as.matrix(data)
  n <- dim(data)[1]
  p <- dim(data)[2]
  centered <- apply(data, 2, function(x)
    x - mean(x))
  sigma_sample <- t(centered) %*% centered / (n - 1)
  eigen_tmp <- eigen(sigma_sample)
  eigenval <- eigen_tmp$values
  eigenvec <- eigen_tmp$vectors
  k <- min(n, p)
  eigenval_subset <- eigenval[1:k]
  eigenval_sumdiffs <-
    rowSums(1 / (
      outer(eigenval_subset, eigenval_subset, "-") + diag(Inf, length(eigenval_subset))
    ))
  eigenval_stein <-
    eigenval_subset * n / (abs(n - p) + 1 + 2 * eigenval_subset * eigenval_sumdiffs)
  if (n < p) {
    eigenval_stein <- c(eigenval_stein, rep(0, p - n))
  }
  eigenval_haff <- 1 / stats::isoreg(1 / eigenval_stein[1:k])$yf
  if (n < p) {
    eigenval_haff <- c(eigenval_haff, rep(0, p - n))
  }

  sigma_mat <- eigenvec %*% diag(eigenval_haff) %*% t(eigenvec)

  rownames(sigma_mat) <- colnames(data)
  colnames(sigma_mat) <- colnames(data)

  return(list(sigma_mat, NA))
}

#' Eigenvalue Clipping Covariance Estimation (Marcenko-Pastur)
#'
#' Computes the eigenvalue clipping estimator of the covariance matrix with the Marcenko-Pastur edge.
#'
#' @param data an nxp data matrix
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here an NA.
#' }
#' @details  The eigenvalue clipping covariance matrix estimator is computed with the following formula:
#' \deqn{\hat{\Sigma}=\Delta\hat{\Lambda}\Delta',}
#' where \eqn{\Delta} is the matrix with the sample eigenvectors of the data matrix and \eqn{\hat{\Lambda}} is a diagonal matrix with the "clipped" sample eigenvalues.
#' The clipping procedure follows \insertCite{laloux1999;textual}{CovEstim}. In particular, when assuming i.i.d returns, the eigenvalues of the sample correlation matrix
#' are distributed according to a Marcenko-Pastur distribution \insertCite{marvcenko1967distribution}{CovEstim} with
#' \deqn{\lambda_{min, max}=(1\mp\sqrt{p/n})^2}
#' as the smallest and largest eigenvalues of a random correlation matrix. Therefore, only eigenvalues which lie outside this interval can bring useful information.
#' In this eigenvalue clipping procedure the sample eigenvalues bigger that \eqn{\lambda_{max}} are kept and the rest are substituted with their average as in \insertCite{bouchaudpotters2009;textual}{CovEstim}.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_evc_mp <- sigma_estim_evc_mp(sp_rets)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#'\insertAllCited
#'
#' @export sigma_estim_evc_mp
#'
sigma_estim_evc_mp <- function(data) {
  data <- as.matrix(data)
  n <- dim(data)[1]
  p <- dim(data)[2]
  vola_mat <- diag(apply(data, 2, stats::sd, na.rm = TRUE))
  corr_mat <- stats::cor(data)
  eigen_tmp <- eigen(corr_mat)
  eigenval <- eigen_tmp$values
  eigenvec <- eigen_tmp$vectors
  meanvar <- mean(eigenval)
  lambdamax <- meanvar * ((1 + sqrt(p / n)) ^ 2)
  index <- which(eigenval < lambdamax)
  eigenval_aver <- mean(eigenval[index])
  eigenval_rmt <- eigenval
  eigenval_rmt[index] <- eigenval_aver

  corr_mat_rmt <- eigenvec %*% diag(eigenval_rmt) %*% t(eigenvec)

  sigma_mat <- as.matrix(vola_mat %*% corr_mat_rmt %*% vola_mat)

  rownames(sigma_mat) <- colnames(data)
  colnames(sigma_mat) <- colnames(data)

  return(list(sigma_mat, NA))
}

#'  Eigenvalue Clipping Covariance Estimation (Bouchaud-Potters)
#'
#' Computes the eigenvalue clipping estimator of the covariance matrix with the Bouchaud-Potters technique.
#'
#' @param data an nxp data matrix
#' @param cut_edge a double, indicating the proportion for the applied eigenvalue clipping.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the user-supplied cut edge.
#' }
#' @details The eigenvalue clipping covariance matrix estimator is computed with the following formula:
#' \deqn{\hat{\Sigma}=\Delta\hat{\Lambda}\Delta',}
#' where \eqn{\Delta} is the matrix with the sample eigenvectors of the data matrix and \eqn{\hat{\Lambda}} is a diagonal matrix with the clipped sample eigenvalues.
#' The clipping procedure follows \insertCite{bouchaudpotters2009;textual}{CovEstim}.
#' In particular, the user-defined cutting edge \eqn{s} gives the proportion for the applied eigenvalue clipping so that the \eqn{(1-s)\times p} largest eigenvalues are kept
#' and the remaining \eqn{s\times p} eigenvalues are substituted by their average.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_evc_bp <- sigma_estim_evc_bp(sp_rets, cut_edge=0.3)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#'\insertAllCited
#'
#' @export sigma_estim_evc_bp
#'
sigma_estim_evc_bp <- function(data, cut_edge) {
  data <- as.matrix(data)
  p <- dim(data)[2]
  vola_mat <- diag(apply(data, 2, stats::sd, na.rm = TRUE))
  corr_mat <- stats::cor(data)
  eigen_tmp <- eigen(corr_mat)
  eigenval <-  eigen_tmp$values
  eigenvec <- eigen_tmp$vectors

  keep_edge <- (1 - cut_edge) * p
  if (keep_edge == 0) {
    eigenval_rmt <- rep.int(1, p)
  } else if (keep_edge == p) {
    eigenval_rmt <- eigenval
  } else{
    eigenval_rmt <- eigenval
    eigenval_rmt[(keep_edge + 1):p] <-
      mean(eigenval[(keep_edge + 1):p])
  }
  corr_mat_rmt <- eigenvec %*% diag(eigenval_rmt) %*% t(eigenvec)
  sigma_mat <- as.matrix(vola_mat %*% corr_mat_rmt %*% vola_mat)

  rownames(sigma_mat) <- colnames(data)
  colnames(sigma_mat) <- colnames(data)

  return(list(sigma_mat, cut_edge))
}

#'  PCA Covariance Estimation
#'
#' Computes a PCA estimator of the covariance matrix.
#'
#' @param data an nxp data matrix
#' @param number_pc an integer, indicating the number of principal components. Default value is NULL and the
#' number of principal components is set according to the Marcenko-Pastur edge as in \insertCite{marvcenko1967distribution;textual}{CovEstim}.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the number of principal components.
#' }
#' @examples
#' data(sp200)
#' sp_rets <- sp200[1:100,-1]
#' sigma_pca2 <- sigma_estim_pca(sp_rets, number_pc=2)[[1]] # user-defined number of factors
#' results_pca_mp <- sigma_estim_pca(sp_rets) # number of factors, defined with MP cut-edge
#' sigma_pca_mp <- results_pca_mp[[1]]
#' number_pc <- results_pca_mp[[2]]
#'
#'@importFrom Rdpack reprompt
#' @references
#'\insertRef{johnson2002applied}{CovEstim}
#'
#'\insertRef{fan2016overview}{CovEstim}
#'
#' @export sigma_estim_pca
#'
sigma_estim_pca <- function(data, number_pc = NULL) {
  data <- as.matrix(data)
  centered <- apply(data, 2, function(x)
    x - mean(x))
  centered_t <- t(centered)
  p <- dim(centered_t)[1]
  n <- dim(centered_t)[2]

  if (is.null(number_pc)) {
    corr_mat <- stats::cor(centered)
    cor_eigen_tmp <- eigen(corr_mat)
    cor_eigenval <- cor_eigen_tmp$values
    lambdamax <- ((1 + sqrt(p / n)) ^ 2)
    index <- which(cor_eigenval  >= lambdamax)
    if (length(index) != 0) {
      number_pc <- index[length(index)]
    } else{
      number_pc <- 0
    }
  }

  if (number_pc != 0) {
    eigen_tmp <- eigen(t(centered_t) %*% centered_t)
    eigenvec <- eigen_tmp$vectors
    fac_pca <- sqrt(n) * eigenvec[, 1:number_pc]
    lam_pca <- centered_t %*% fac_pca / n
    uhat <-  centered_t - lam_pca %*% t(fac_pca)
    low_rank <- lam_pca %*% t(lam_pca)
    su_pca <- uhat %*% t(uhat) / n
    su_diag <- diag(diag(su_pca))

    sigma_mat <- low_rank + su_diag
  } else{
    sigma_mat <- diag(1, p)
  }
  rownames(sigma_mat) <- colnames(data)
  colnames(sigma_mat) <- colnames(data)

  return(list(sigma_mat, number_pc))

}
