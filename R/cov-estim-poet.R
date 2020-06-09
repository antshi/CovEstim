#' Principal Orthogonal ComplEment Thresholding (POET) Covariance Estimation
#'
#' Computes the POET estimator of the covariance matrix.
#'
#' @param data an nxp data matrix.
#' @param K an integer, the number of principal components (factors).
#' Default values is NULL and the optimal K is calculated as in \insertCite{poetpackage;textual}{CovEstim}.
#' K=0 corresponds to threshoding the sample covariance directly.
#' @param C a double, the positive constant for thresholding. Default value is NULL and the optimal K is calculated as in \insertCite{poetpackage;textual}{CovEstim}.
#' @param thres a character, indicating the choice of thresholding. Possible values are "soft" for soft-thresholding, "hard" for hard thresholding and
#' "scad" for scad thresholding. Default is "soft".
#' @param matrix a character, indicating the option of thresholding either correlation or covairance matrix.
#' Possible values are "cor" for thresholding the error correlation matrix then transform back to covariance matrix and "vad" for thresholding the error covariance matrix directly.
#' Default is "cor".
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the number of principal components K.
#' }
#'
#' @details The POET estimator of the covariance matrix is computed according to \insertCite{fan2013poet;textual}{CovEstim}.
#' The POET estimation is originally found under \insertCite{poetpackage;textual}{CovEstim}.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_poet <- sigma_estim_poet(sp_rets, K=2)[[1]] # user-defined K
#' sigma_poet <- sigma_estim_poet(sp_rets)[[1]] # optimal K
#'
#'
#' @importFrom Rdpack reprompt
#' @references
#'\insertAllCited
#'
#' @export sigma_estim_poet
#'
sigma_estim_poet <-
  function(data,
           K = NULL,
           C = NULL,
           thres = "soft",
           matrix = "vad") {
    if (is.null(K)) {
      result <- poet_estim(data,
                           C = C,
                           thres = thres,
                           matrix = matrix)
      sigma_mat <- result$SigmaY
      K <- result$K
    } else{
      sigma_mat <-
        poet_estim(
          data,
          K = K,
          C = C,
          thres = thres,
          matrix = matrix
        )$SigmaY
    }

    return(list(sigma_mat, K))
  }


##
poet_estim <-
  function(data,
           K = NULL,
           C = NULL,
           thres = "soft",
           matrix = "vad") {
    data <- as.matrix(data)
    centered <- apply(data, 2, function(x)
      x - mean(x))
    Y <- t(centered)
    p <- dim(Y)[1]
    n <- dim(Y)[2]

    if (is.null(K)) {
      K1 <- mean(unlist(poet_kopt(data)))
      K <- floor(K1) + 1
    }

    if (K > 0) {
      eigen_tmp <- eigen(t(Y) %*% Y)
      V <- eigen_tmp$vectors
      F <- sqrt(n) * V[, 1:K]
      LamPCA <- Y %*% F / n
      uhat <- Y - LamPCA %*% t(F)
      Lowrank <- LamPCA %*% t(LamPCA)
      rate <- 1 / sqrt(p) + sqrt((log(p)) / n)
    } else{
      uhat <- Y
      rate <- sqrt((log(p)) / n)
      Lowrank <- matrix(0, p, p)
    }

    SuPCA <- uhat %*% t(uhat) / n
    SuDiag <- diag(diag(SuPCA))

    if (matrix == "cor") {
      SuDiagSqrt <- SuDiag ^ (1 / 2)
      SuDiagSqrt_inv <- solve(SuDiagSqrt)
      R <- SuDiagSqrt_inv %*% SuPCA %*% SuDiagSqrt_inv
    }
    if (matrix  == "vad") {
      R <- SuPCA
    }

    if (is.null(C)) {
      C1 <- poet_copt(data, K, thres, matrix)
      C <- C1 + 0.1
    }

    uu <- array(0, dim = c(p, p, n))
    roottheta <- array(0, dim = c(p, p))
    lambda <- array(0, dim = c(p, p))

    for (i in 1:p) {
      for (j in 1:i) {
        uu[i, j,] <- uhat[i,] * uhat[j,]
        roottheta[i, j] <- stats::sd(uu[i, j,])
        lambda[i, j] <- roottheta[i, j] * rate * C
        lambda[j, i] <- lambda[i, j]
      }
    }

    Rthresh <- matrix(0, p, p)

    if (thres == "soft") {
      for (i in 1:p) {
        for (j in 1:i) {
          if (abs(R[i, j]) < lambda[i, j] && j < i) {
            Rthresh[i, j] <- 0
          } else{
            if (j == i) {
              Rthresh[i, j] <- R[i, j]
            } else {
              Rthresh[i, j] <- sign(R[i, j]) * (abs(R[i, j]) - lambda[i, j])

            }
          }
          Rthresh[j, i] <- Rthresh[i, j]

        }
      }
    }


    if (thres == "hard") {
      for (i in 1:p) {
        for (j in 1:i) {
          if (abs(R[i, j]) < lambda[i, j] && j < i) {
            Rthresh[i, j] <- 0
          } else{
            Rthresh[i, j] <- R[i, j]
          }

          Rthresh[j, i] <- Rthresh[i, j]
        }
      }
    }



    if (thres == "scad") {
      for (i in 1:p) {
        for (j in 1:i) {
          if (j == i) {
            Rthresh[i, j] <- R[i, j]
          } else{
            if (abs(R[i, j]) < lambda[i, j]) {
              Rthresh[i, j] <- 0
            } else{
              if (abs(R[i, j]) < 2 * lambda[i, j]) {
                Rthresh[i, j] <- sign(R[i, j]) * (abs(R[i, j]) - lambda[i, j])
              } else{
                if (abs(R[i, j]) < 3.7 * lambda[i, j]) {
                  Rthresh[i, j] <-
                    ((3.7 - 1) * R[i, j] - sign(R[i, j]) * 3.7 * lambda[i, j]) / (3.7 - 2)
                } else{
                  Rthresh[i, j] <- R[i, j]
                }
              }
            }
          }
          Rthresh[j, i] <- Rthresh[i, j]
        }
      }
    }

    SigmaU <- matrix(0, p, p)

    if (matrix == "cor") {
      SigmaU <- SuDiagSqrt %*% Rthresh * SuDiagSqrt
    }
    if (matrix == "vad") {
      SigmaU <- Rthresh
    }


    SigmaY <- Lowrank + SigmaU


    if (K > 0) {
      result <-
        list(
          SigmaU = SigmaU,
          SigmaY = SigmaY,
          factors = t(F),
          loadings = LamPCA,
          K = K,
          C = C
        )

    } else{
      result <- list(
        SigmaU = SigmaU,
        SigmaY = SigmaY,
        K = K,
        C = C
      )
    }

    return(result)
  }


##
poet_copt <- function(data,
                      K,
                      thres = "soft",
                      matrix = "vad") {
  mineig <- function(x) {
    sigma_u <-
      poet_estim(
        data = data,
        K = K,
        x,
        thres = thres,
        matrix = matrix
      )$SigmaU
    f <- min(eigen(sigma_u)$values)
    return(c(f))
  }

  if (mineig(50) * mineig(-50) < 0) {
    return(max(0, stats::uniroot(mineig, c(-50, 50), tol = 0.001)$root))
  } else{
    return(0)
  }
}

##
poet_kopt <- function(data, reps = 20, k_max = 10) {
  data <- as.matrix(data)
  centered <- apply(data, 2, function(x)
    x - mean(x))
  Y <- t(centered)
  p <- dim(Y)[1]
  n <- dim(Y)[2]

  # Hallin and Liska method
  penalty_seq <-
    seq(0.05, 5, length = 100) # or seq(0.01, 3, by=0.01) #as in Hallin, Liska (2007)
  penalty_len <- length(penalty_seq)

  IC <- array(, c(2, reps, k_max, penalty_len))
  gT1HL <- array(, c(reps))
  gT2HL <- array(, c(reps))
  pi <- array(, c(reps))
  ni <- array(, c(reps))

  # generate the subsets (reps)
  for (i in 1:reps) {
    pi <- min(i * floor(p / reps) + min(p, 5), p)
    ni <- min(i * floor(n / reps) + min(n, 5), n)
    if (i == reps) {
      pi <- p
      ni <- n
    }

    Yi <- Y[1:pi, 1:ni]
    eigen_tmpi <- eigen(t(Yi) %*% Yi)
    Vi <- as.matrix(eigen_tmpi$vectors)
    gT1HL[i] <- log((pi * ni) / (pi + ni)) * (pi + ni) / (pi * ni)
    gT2HL[i] <- log(min(pi, ni)) * (pi + ni) / (pi * ni)

    frob <- array(, c(k_max))
    minval <- min(pi, ni, k_max)

    for (k in 1:minval) {
      F <- Vi[, 1:k]
      LamPCA <- Yi %*% F / ni
      uhat <- Yi - LamPCA %*% t(F)
      frob[k] <- sum(diag(uhat %*% t(uhat))) / (pi * ni)

      for (l in 1:penalty_len) {
        # only fills in the ICs up to k, which may be <k_max
        IC[1, i, k, l] <-
          log(frob[k]) + penalty_seq[l] * k * gT1HL[i]
        IC[2, i, k, l] <-
          log(frob[k]) + penalty_seq[l] * k * gT2HL[i]
      }
    }
  }

  rhat <- array(, c(2, reps, penalty_len))
  for (i in 1:reps) {
    for (l in 1:penalty_len) {
      rhat[1, i, l] <- which.min(IC[1, i, 1:minval, l])
      rhat[2, i, l] <- which.min(IC[2, i, 1:minval, l])
    }
  }

  Sc1 <- array(, c(penalty_len))
  Sc2 <- array(, c(penalty_len))

  for (l in 1:penalty_len) {
    Sc1[l] <- stats::sd(rhat[1, , l])
    Sc2[l] <- stats::sd(rhat[2, , l])
  }

  K1HL <- rhat[1, 1, which(Sc1 == 0)[1]]
  K2HL <- rhat[2, 1, which(Sc2 == 0)[1]]

  # Bai and Ng method - penalty corresponds to c=1

  gT1BN <- log((p * n) / (p + n)) * (p + n) / (p * n)
  gT2BN <- log(min(p, n)) * (p + n) / (p * n)

  ICBN <- matrix(, 2, k_max)
  logfrob <- c()

  eigen_tmp <- eigen(t(Y) %*% Y)
  V <- as.matrix(eigen_tmp$vectors)

  for (k in 1:k_max) {
    F <- V[, 1:k]
    LamPCA <- Y %*% F / n
    uhat <- Y - LamPCA %*% t(F)
    logfrob[k] <- log(sum(diag(uhat %*% t(uhat))) / (p * n))

    ICBN[1, k] <- logfrob[k] + k * gT1BN
    ICBN[2, k] <- logfrob[k] + k * gT2BN
  }

  K1BN <- which.min(ICBN[1, ])
  K2BN <- which.min(ICBN[2, ])

  return(list(
    K1HL = K1HL,
    K2HL = K2HL,
    K1BN = K1BN,
    K2BN = K2BN
  ))

}
