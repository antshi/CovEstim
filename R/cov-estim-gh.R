#' Graphical Horseshoe Covariance Estimation
#'
#' Computes the Graphical Horseshoe estimator of the covariance matrix.
#'
#' @param data an nxp data matrix.
#' @param burnin an integer, indicating the number of MCMC burnins.
#' @param nmc an integer, indicating the number of samples.
#'
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here NA.
#' }
#'
#' @details The Graphical Horseshoe estimator is elaborated in detail in \insertCite{li2019graphical;textual}{CovEstim}.
#' More information on the functionality can be found in the Matlab package \url{https://github.com/liyf1988/GHS}.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_ghs <- sigma_estim_ghs(sp_rets)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#'\insertAllCited
#'
#' @export sigma_estim_glasso
#'

sigma_estim_ghs <- function(data, burnin = 100, nmc = 5000) {
  data <- as.matrix(data)
  n <- dim(data)[1]
  p <- dim(data)[2]
  centered <- apply(data, 2, function(x)
    x - mean(x))
  S <- t(centered) %*% centered

  omega_save <- array(0, c(p, p, nmc))
  lambda_sq_save <- array(0, c(p * (p - 1) / 2, nmc))
  tau_sq_save <- matrix(0, nmc, nmc)

  Omega <- diag(1, p)
  Sigma <- diag(1, p)
  Lambda_sq <- matrix(1, p, p)
  Nu <- matrix(1, p, p)
  tau_sq <- 1
  xi <- 1

  for (iter in 1:(burnin + nmc)) {
    for (i in 1:p) {
      if (i == 1) {
        ind <- 2:p
      } else if (i == p) {
        ind <- 1:(p - 1)
      } else{
        ind <- c(1:(i - 1), (i + 1):p)
      }

      Sigma_11 <- Sigma[ind, ind]
      sigma_12 <- Sigma[ind, i]
      sigma_22 <- Sigma[i, i]
      s_21 <- S[ind, i]
      s_22 <- S[i, i]
      lambda_sq_12 <- Lambda_sq[ind, i]
      nu_12 <- Nu[ind, i]
      gamma <- stats::rgamma(1, shape = n / 2 + 1, rate = s_22 / 2)
      inv_Omega_11 <-
        Sigma_11 - (sigma_12 %*% t(sigma_12)) / sigma_22
      inv_C <-
        s_22 * inv_Omega_11 + diag(1 / (lambda_sq_12 * tau_sq))

      inv_C_chol <- chol(inv_C)
      mu_i <- solve(-inv_C, s_21)
      beta <- mu_i + solve(inv_C_chol, stats::rnorm(p - 1, 1))
      omega_12 <- beta
      omega_22 <- gamma + t(beta) %*% inv_Omega_11 %*% beta
      rate <- (omega_12 ^ 2) / (2 * tau_sq) + 1 / nu_12
      lambda_sq_12 <-
        1 / (sapply(rate, function(x)
          stats::rgamma(
            1, shape = 1, rate = x
          )))
      nu_12 <-
        1 / (sapply(lambda_sq_12, function(x)
          stats::rgamma(1, shape = 1, rate = 1 + 1 / x)))
      Omega[i, ind] <- omega_12
      Omega[ind, i] <- omega_12
      Omega[i, i] <- omega_22
      temp <- inv_Omega_11 %*% beta
      Sigma_11 <- inv_Omega_11 + temp %*% t(temp) / gamma
      sigma_12 <- -temp / gamma
      sigma_22 <- 1 / gamma
      Sigma[ind, ind] <- Sigma_11
      Sigma[i, i] <- sigma_22
      Sigma[i, ind] <- sigma_12
      Sigma[ind, i] <- sigma_12
      Lambda_sq[i, ind] <- lambda_sq_12
      Lambda_sq[ind, i] <- lambda_sq_12
      Nu[i, ind] <- nu_12
      Nu[ind, i] <- nu_12

    }

    omega_vector <-
      Omega[lower.tri(matrix(TRUE, p, p), diag = FALSE)]
    lambda_sq_vector <-
      Lambda_sq[lower.tri(matrix(TRUE, p, p), diag = FALSE)]
    rate <-
      1 / xi + sum((omega_vector ^ 2) / (2 * lambda_sq_vector))
    tau_sq <-
      1 / stats::rgamma(1, shape = (p * (p - 1) / 2 + 1) / 2, rate = rate)
    xi <- 1 / stats::rgamma(1, shape = 1, rate = (1 + 1 / tau_sq))

    if (iter > burnin) {
      omega_save[, , iter - burnin] <- Omega
      lambda_sq_save[, iter - burnin] <- lambda_sq_vector
      tau_sq_save[iter - burnin] <- tau_sq
    }

  }
  sigma_mat <- as.matrix(apply(omega_save, 3, mean))
  rownames(sigma_mat) <- colnames(data)
  colnames(sigma_mat) <- colnames(data)
  return(list(sigma_mat, NA))

}
