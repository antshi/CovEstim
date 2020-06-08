#' Nonparametric eigenvalue-regularized Precision or Covariance Matrix Estimation (NERCOME)
#'
#' Computes the NERCOME estimator of the covariance matrix.
#'
#' @param data an nxp data matrix.
#' @param split an integer, indicating the position of the split.
#' Default value is NULL and an optimal split is calculated with permutation along a sequence as in \insertCite{lam2016nonparametric;textual}{CovEstim}.
#'
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the split.
#' }
#'
#'@details More information on the NERCOME estimator is found in \insertCite{lam2016nonparametric;textual}{CovEstim}.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_nercome <- sigma_estim_nercome(sp_rets)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#'\insertAllCited
#'
#' @export sigma_estim_nercome
#'
sigma_estim_nercome <- function(data, split = NULL) {
  data <- as.matrix(data)
  n <- dim(data)[1]
  centered <- apply(data, 2, function(x)
    x - mean(x))

  if (is.null(split)) {
    split_seq <-
      round(c(
        2 * sqrt(n),
        0.2 * n,
        0.4 * n,
        0.6 * n,
        0.8 * n,
        n - 2.5 * sqrt(n),
        n - 1.5 * sqrt(n)
      ),
      0)
    critval_perm <- matrix(, nrow = 500, ncol = length(split_seq))

    for (i in 1:500) {
      centered_perm <- centered[sample(1:n, n, replace = FALSE), ]

      for (s in seq_along(split_seq)) {
        split_ind <- seq_len(split_seq[s])
        y1 <- centered_perm[split_ind, ]
        nm <- dim(y1)[1]
        y2 <- centered_perm[-split_ind, ]
        nk <- dim(y2)[1]
        sigma1 <- (1 / nm) * t(y1) %*% y1
        eigenvec_m1 <- eigen(sigma1)$vectors
        sigma2 <- (1 / nk) * t(y2) %*% y2
        sigma_mat <-
          eigenvec_m1 %*% diag(diag(t(eigenvec_m1) %*% sigma2 %*% eigenvec_m1)) %*%
          t(eigenvec_m1)
        critval_perm[i, s] <- sum((sigma_mat - sigma2) ^ 2)
      }
    }

    critval <- colMeans(critval_perm, na.rm = TRUE)
    split <- split_seq[which.min(critval)]

  }
  split_ind <- 1:split
  y1 <- centered[split_ind, ]
  nm <- dim(y1)[1]
  y2 <- centered[-split_ind, ]
  nk <- dim(y2)[1]
  sigma1 <- (1 / nm) * t(y1) %*% y1
  eigenvec_m1 <- eigen(sigma1)$vectors
  sigma2 <- (1 / nk) * t(y2) %*% y2
  sigma_mat <-
    eigenvec_m1 %*% diag(diag(t(eigenvec_m1) %*% sigma2 %*% eigenvec_m1)) %*%
    t(eigenvec_m1)

  rownames(sigma_mat) <- colnames(data)
  colnames(sigma_mat) <- colnames(data)

  return(list(sigma_mat, split))
}
