#' Simulated Covariance Matrix
#'
#' Simulates a covariance matrix with specific properties
#'
#' @param p an integer, indicating the number of variables.
#' @param corr_min a double, indicating the minimum possible correlation across the p variables.
#' @param corr_max a double, indicating the maximum possible correlation across the p variables (except the correlation with self of 1).
#' @param vola_min a double, indicating the minimum possible volatility of the p variables.
#' @param vola_max a double, indicating the maximum possible volatility of the p variables.
#' @return a pxp covariance matrix.
#' @details The correlation and volatility values are drown from a uniform distribution.
#' A condition for positive-definiteness of the resulting covariance matrix is not implemented.
#'
#' @examples
#' sim_sigma <- sigma_sim(100, corr_min=0.01, corr_max=0.40, vola_min=0.05, vola_max=0.30)
#'
#' @export sigma_sim
#'
sigma_sim <- function(p, corr_min, corr_max, vola_min, vola_max) {
    corr_mat <- diag(0.5, p)
    corr_mat[lower.tri(corr_mat)] <-
        stats::runif((p ^ 2 - p) / 2, min = corr_min, max = corr_max)
    corr_mat <- corr_mat + t(corr_mat)
    vola_mat <- diag(stats::runif(p, min = vola_min, max = vola_max))

    sigma_mat <- as.matrix(vola_mat %*% corr_mat %*% vola_mat)

    return(sigma_mat)

}
