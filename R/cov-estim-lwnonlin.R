#' Ledoit-Wolf Covariance Estimation (Nonlinear Shrinkage)
#'
#' Computes the analytical Ledoit-Wolf nonlinear shrinkage estimator of the covariance matrix.
#'
#' @param data an nxp data matrix.
#' @param bandwidth_speed a double, indicating the speed at which the bandwidth vanishes in the number of variables p.
#' Default value is -0.35.
#' @param zeromean_log a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE). Default value is FALSE.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the bandwidth speed.
#' }
#'
#' @details The Ledoit-Wolf nonlinear shrinkage estimator of the covariance matrix is computed according to \insertCite{ledoit2018analytical;textual}{CovEstim}
#' with the following formula:
#' \deqn{\hat{\Sigma}=\Delta\hat{\Lambda}\Delta',}
#' where \eqn{\Delta} is the matrix with the sample eigenvectors of the data matrix and \eqn{\hat{\Lambda}} is a diagonal matrix with the sample eigenvalues, shrunk in a nonlinear way.
#' The optimal solution is achieved using a nonparametric variable bandwidth kernel estimation of the limiting spectral density of the sample eigenvalues and its Hilbert transform.
#' The speed at which the bandwidth vanishes in the number of assets is set to -0.35.
#' A corresponding MATLAB code for the estimator can be accessed under \url{https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html}.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_lwnl <- sigma_estim_lwnl(sp_rets)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#'\insertAllCited
#'
#' @export sigma_estim_lwnl
#'
sigma_estim_lwnl <-
  function(data,
           bandwidth_speed = NULL,
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

    sigma_ml <- t(centered) %*% centered / n

    eigen_tmp <- eigen(sigma_ml)
    eigenval <- eigen_tmp$values
    eigenvec <- eigen_tmp$vectors
    eigenval_sort <- sort(eigenval, index.return = TRUE)$x
    eigenvec_sort <-
      eigenvec[, sort(eigenval, index.return = TRUE)$ix]
    lambda <- as.matrix(eigenval_sort[max(1, p - n + 1):p])
    ident_mat <-
      matrix(rep.int(1, min(p, n) * 1), nrow = 1, ncol = min(p, n))
    lambda_mat <- kronecker(ident_mat, lambda)
    if(is.null(bandwidth_speed)){
        bandwidth_speed <- -0.35
    }
    h <- n ^ (bandwidth_speed)
    h_mat <- h * t(lambda_mat)
    x <- (lambda_mat - t(lambda_mat)) / h_mat
    f_tilde <-
      (3 / 4 / sqrt(5)) * rowMeans(pmax(1 - (x ^ 2) / 5, 0) / h_mat)
    hf_temp <-
      (-3 / 10 / pi) * x + (3 / 4 / sqrt(5) / pi) * (1 - (x ^ 2) / 5) * log(abs((sqrt(5) - x) / (sqrt(5) + x)))
    hf_temp[abs(x) == sqrt(5)] <-
      (-3 / 10 / pi) * x[abs(x) == sqrt(5)]
    hf_tilde <- rowMeans(hf_temp / h_mat)

    if (p <= n) {
      dtilde <-
        lambda / (((pi * (p / n) * lambda * f_tilde) ^ 2) + ((1 - (p / n) - pi *
                                                                (p / n) * lambda * hf_tilde) ^ 2))
    } else{
      hf_tilde0 <-
        (1 / pi) * (3 / 10 / (h ^ 2) + 3 / 4 / sqrt(5) / h * (1 - 1 / 5 / (h ^ 2)) * log((1 + sqrt(5) * h) / (1 - sqrt(5) * h))) * mean(1 / lambda)
      dtilde0 <- 1 / (pi * ((p - n) / n) * hf_tilde0)
      dtilde1 <-
        lambda / ((pi ^ 2) * (lambda ^ 2) * (((f_tilde) ^ 2) + (hf_tilde ^ 2)))
      dtilde <- c(rep(dtilde0, p - n), dtilde1)
    }

    sigma_mat <-
      eigenvec_sort %*% diag(as.numeric(dtilde)) %*% t(eigenvec_sort)

    rownames(sigma_mat) <- colnames(data)
    colnames(sigma_mat) <- colnames(data)

    return(list(sigma_mat, bandwidth_speed))
  }
