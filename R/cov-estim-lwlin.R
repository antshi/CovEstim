#' Ledoit-Wolf Linear Shrinkage Covariance Estimation I
#'
#' Computes the Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards the identity matrix.
#'
#' @param data an nxp data matrix.
#' @param shrink_int a double, indicating the shrinkage intensity. Default is the optimal shrinkage intensity as in \insertCite{ledoit2003identity;textual}{CovEstim}.
#' @param zeromean_log a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE). Default value is FALSE.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the shrinkage intensity.
#' }
#'
#' @details The Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards the identity matrix is calculated with the following formula:
#' \deqn{\hat{\Sigma}= s\Sigma_{T} + (1-s)\Sigma,}
#' where \eqn{\Sigma} is the sample covariance matrix, s is the user-supplied or optimal shrinkage intensity and \eqn{\Sigma_{T}} is a pxp identity matrix.
#' This covariance estimator assumes a zero correlation and variances of one as the underlying covariance structure of the data.
#' A corresponding MATLAB code for the estimator can be accessed under \url{https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html}.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_lwident <- sigma_estim_lwident(sp_rets)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#'\insertAllCited
#'
#' @export sigma_estim_lwident
#'
sigma_estim_lwident <-
  function(data,
           shrink_int = NULL,
           zeromean_log = FALSE) {
    data <- as.matrix(data)
    p <- dim(data)[2]
    if (!zeromean_log) {
      centered <- apply(data, 2, function(x)
        x - mean(x))
      n <- dim(data)[1] - 1
    } else{
      centered <- data
      n <- dim(data)[1]
    }

    sigma_sample <- t(centered) %*% centered / n

    sigma_target <- diag(1, p)

    if (is.null(shrink_int)) {
      asyvar_mat <-
        t(centered ^ 2) %*% (centered ^ 2) / n - 2 * t(centered) %*% (centered) *
        sigma_sample / n + sigma_sample ^ 2
      asyvar <- sum(asyvar_mat)
      gamma <- sum((sigma_target - sigma_sample) ^ 2)
      kappa <- asyvar / gamma
      shrink_int <- max(0, min(1, kappa / n))
    }

    sigma_mat <-
      shrink_int * sigma_target + (1 - shrink_int) * sigma_sample

    rownames(sigma_mat) <- colnames(data)
    colnames(sigma_mat) <- colnames(data)

    return(list(sigma_mat, shrink_int))
  }


#' Ledoit-Wolf Linear Shrinkage Covariance Estimation II
#'
#' Computes the Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards the one-parameter matrix.
#'
#' @param data an nxp data matrix.
#' @param shrink_int a double, indicating the shrinkage intensity. Default is the optimal shrinkage intensity as in \insertCite{ledoit2004oneparam;textual}{CovEstim}.
#' @param zeromean_log a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE). Default value is FALSE.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the shrinkage intensity.
#' }
#'
#' @details The Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards the diagonal matrix of equal variances is calculated with the following formula:
#' \deqn{\hat{\Sigma}= s\Sigma_{T} + (1-s)\Sigma,}
#' where \eqn{\Sigma} is the sample covariance matrix, s is the user-supplied or optimal shrinkage intensity and
#' \eqn{\Sigma_{T}} is a diagonal matrix with the average sample variance \eqn{\bar{\sigma}^2} on the diagonal.
#' This covariance estimator assumes a zero correlation and equal variances as the underlying covariance structure of the data.
#' A corresponding MATLAB code for the estimator can be accessed under \url{https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html}.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_lwone <- sigma_estim_lwone(sp_rets)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#'\insertAllCited
#'
#' @export sigma_estim_lwone
#'
sigma_estim_lwone <-
  function(data,
           shrink_int = NULL,
           zeromean_log = FALSE) {
    data <- as.matrix(data)
    p <- dim(data)[2]
    if (!zeromean_log) {
      centered <- apply(data, 2, function(x)
        x - mean(x))
      n <- dim(data)[1] - 1
    } else{
      centered <- data
      n <- dim(data)[1]
    }

    sigma_sample <- t(centered) %*% centered / n

    aver_var <- mean(diag(sigma_sample))
    sigma_target <- aver_var * diag(1, p)

    if (is.null(shrink_int)) {
      asyvar_mat <-
        t(centered ^ 2) %*% (centered ^ 2) / n - 2 * t(centered) %*% (centered) *
        sigma_sample / n + sigma_sample ^ 2
      asyvar <- sum(asyvar_mat)
      gamma <- sum((sigma_target - sigma_sample) ^ 2)
      kappa <- asyvar / gamma
      shrink_int <- max(0, min(1, kappa / n))
    }

    sigma_mat <-
      shrink_int * sigma_target + (1 - shrink_int) * sigma_sample

    rownames(sigma_mat) <- colnames(data)
    colnames(sigma_mat) <- colnames(data)

    return(list(sigma_mat, shrink_int))
  }

#' Ledoit-Wolf Linear Shrinkage Covariance Estimation III
#'
#' Computes the Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards the constant correlation covariance matrix.
#'
#' @param data an nxp data matrix.
#' @param shrink_int a double, indicating the shrinkage intensity. Default is the optimal shrinkage intensity as in \insertCite{ledoit2004cc;textual}{CovEstim}.
#' @param zeromean_log a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE). Default value is FALSE.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the shrinkage intensity.
#' }
#'
#' @details The Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards the constant correlation covariance matrix is calculated with the following formula:
#' \deqn{\hat{\Sigma}= s\Sigma_{T} + (1-s)\Sigma,}
#' where \eqn{\Sigma} is the sample covariance matrix, s is the user-supplied or optimal shrinkage intensity and
#' \eqn{\Sigma_{T}} is the constant correlation covariance matrix.
#' This covariance estimator assumes a constant correlation and individual variances as the underlying covariance structure of the data.
#' A corresponding MATLAB code for the estimator can be accessed under \url{https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html}.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_lwcc <- sigma_estim_lwcc(sp_rets)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#'\insertAllCited
#'
#' @export sigma_estim_lwcc
#'
sigma_estim_lwcc <-
  function(data,
           shrink_int = NULL,
           zeromean_log = FALSE) {
    data <- as.matrix(data)
    p <- dim(data)[2]
    if (!zeromean_log) {
      centered <- apply(data, 2, function(x)
        x - mean(x))
      n <- dim(data)[1] - 1
    } else{
      centered <- data
      n <- dim(data)[1]
    }

    sigma_sample <- t(centered) %*% centered / n
    corr_mat <- stats::cor(data)
    rbar <- (sum(colSums(corr_mat)) - p) / (p * (p - 1))
    sigma_target <- rbar * sigma_sample / corr_mat
    diag(sigma_target) <- diag(sigma_sample)

    if (is.null(shrink_int)) {
      asyvar_mat <-
        t(centered ^ 2) %*% (centered ^ 2) / n - 2 * t(centered) %*% (centered) *
        sigma_sample / n + sigma_sample ^ 2
      asyvar <- sum(asyvar_mat)
      term1 <- t((centered) ^ 3) %*% centered
      term2 <- diag(sigma_sample) * (t(centered) %*% centered)
      term3 <-
        sigma_sample * (t(centered ^ 2) %*% matrix(rep(1, p * nrow(centered)), ncol =
                                                     p))
      term4 <- (diag(sigma_sample) %o% rep(1, p)) * sigma_sample
      term_all <- (term1 - term2 - term3 + term4) / n
      ratios <-
        (diag(sigma_sample) %o% diag(sigma_sample) ^ -1) ^ 0.5
      rhos <-
        0.5 * rbar * (ratios * term_all + t(ratios) * t(term_all))
      rho <-
        sum(diag(asyvar_mat), na.rm = TRUE) + sum(rhos[lower.tri(rhos)], na.rm =
                                                    TRUE) + sum(rhos[upper.tri(rhos)], na.rm = TRUE)
      gamma <- sum((sigma_target - sigma_sample) ^ 2)
      kappa <- (asyvar - rho) / gamma
      shrink_int <- max(0, min(1, kappa / n))
    }

    sigma_mat <-
      shrink_int * sigma_target + (1 - shrink_int) * sigma_sample

    rownames(sigma_mat) <- colnames(data)
    colnames(sigma_mat) <- colnames(data)

    return(list(sigma_mat, shrink_int))
  }
