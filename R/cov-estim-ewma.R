#' EWMA Covariance Estimation
#'
#' Computes the Exponentially Moving Average (EWMA) estimator of the covariance matrix.
#'
#' @param data an nxp data matrix.
#' @param lambda a double for the decay parameter \eqn{\lambda}.
#' Default is 0.97 - the standard for monthly returns according to \insertCite{riskmetrics1996;textual}{CovEstim}.
#' If the data consists of daily returns, lambda should be set to 0.94.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the decay parameter \eqn{\lambda}.
#' }
#' @details The EWMA estimator of the covariance matrix for an nxp data matrix X is computed with the following formula:
#' \deqn{\hat{\Sigma}_{t}= (1-\lambda)R_{t}R'_{t} + \lambda\hat{\Sigma}_{t-1},}
#' where \eqn{R_{t}} is the matrix of demeaned returns for the time period \eqn{t} and \eqn{\hat{\Sigma}_{t-1}} is the EWMA covariance estimator for the period \eqn{t-1}.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_ewma <- sigma_estim_ewma(sp_rets)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#'\insertAllCited
#'
#' @export sigma_estim_ewma
#'
sigma_estim_ewma <- function(data, lambda = 0.97) {
    data <- as.matrix(data)
    centered <- apply(data, 2, function(x)
        x - mean(x))
    n <- dim(centered)[1]

    for (i.e in seq_len(n)) {
        if (i.e == 1) {
            sigma_mat <- (centered[i.e, ] %*% t(centered[i.e, ]))
        } else{
            sigma_mat <-
                (1 - lambda) / (1 - lambda ^ n) * (centered[i.e, ] %*% t(centered[i.e, ])) + lambda * (sigma_mat)
        }
    }

    rownames(sigma_mat) <- colnames(data)
    colnames(sigma_mat) <- colnames(data)

    return(list(sigma_mat, lambda))
}
