#' Matrix Sparsity Check
#'
#' Performs a sparsity check on a data matrix.
#'
#' @param mat an nxp data matrix.
#' @return a logical, indicating whether the matrix is sparse (TRUE) or not (FALSE).
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sparsity_logic <- is_sparse(sp_rets)
#'
#' @export is_sparse
#'
is_sparse <- function(mat) {
  mat <- as.matrix(mat)
  n <- dim(mat)[1]
  p <- dim(mat)[2]
  check <- length(which(mat == 0)) > (n * p) / 2
  return(check)
}


#' Matrix Positive-Definiteness Check
#'
#' Performs a tolerance check for positive-definiteness of a symmetric matrix.
#'
#' @param mat a pxp matrix.
#' @param tol a double with tolerance for the eigenvalues. Default value is 1e-8.
#' @return a logical. TRUE, if the matrix mat is positive-definite.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' ml_sigma <- sigma_estim(sp_rets, "ML")
#' is_posdef(ml_sigma)
#'
#' ml_sigma <- sigma_estim(sp_rets[1:100,], "ML")
#' is_posdef(ml_sigma)
#'
#' @export is_posdef
#'
is_posdef <- function(mat, tol = 1e-8) {
  eigenval <- eigen(mat)$values
  n <- dim(mat)[1]
  for (i in 1:n) {
    if (abs(eigenval[i]) < tol) {
      eigenval[i] <- 0
    }
  }
  check <- !(any(eigenval <= 0))
  return(check)
}

#' Matrix Positive-Definiteness
#'
#' Produces a positive-definite matrix (to a certain tolerance). Originally found under \insertCite{corpcorpackage;textual}{CovEstim}.
#'
#' @param mat a pxp matrix.
#' @param tol a double, the tolerance for the relative positiveness of eigenvalues compared to the largest eigenvalue. Default value is tol=1e-8. If NULL, a machine tolerance approximation is applied.
#' @return a positive-definite (to the tolerance of tol) matrix.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[1:100,-1]
#' ml_sigma <- sigma_estim(sp_rets, "ML")
#' is_posdef(ml_sigma)
#' ml_sigma_posdef <- make_posdef(ml_sigma)
#' is_posdef(ml_sigma_posdef)
#'
#' @importFrom Rdpack reprompt
#' @references
#'\insertRef{corpcorpackage}{CovEstim}
#'
#' @export make_posdef
#'
make_posdef <- function(mat, tol = 1e-8) {
  p <- dim(mat)[1]
  eigen_tmp <- eigen(mat, symmetric = TRUE)
  eigenval <- eigen_tmp$values
  eigenvec <- eigen_tmp$vectors
  if (is.null(tol)) {
    tol <- p * max(abs(eigenval)) * .Machine$double.eps
  }
  delta <- 2 * tol
  tau <- pmax(0, delta - eigenval)
  dm <- eigenvec %*% diag(tau, p) %*% t(eigenvec)
  mat_posdef <- mat + dm

  return(mat_posdef)
}


#' Nearest Positive-Definite Matrix
#'
#' Implements the algorithm of \insertCite{higham2002computing;textual}{CovEstim} to compute the nearest positive-definite matrix to an approximate one, typically a correlation or variance-covariance matrix.
#'
#'
#' @param mat a pxp matrix.
#' @param corr a logical, indicating if the result should be a correlation matrix. Default value is FALSE.
#' @param keep_diag a logical, indicating if the result should have the same diagonal as the original matrix mat.
#' @param do_dykstra a logical, indicating if Dykstra's correlation is to be used. Default value is TRUE.
#' @param eig_tol a double, defining the relative positiveness of eigenvalues compared to the largest eigenvalue. Default value is 1e-6.
#' @param conv_tol a double, defining the convergence tolerance for the Higham algorithm. Default value is 1e-7.
#' @param posd_tol a double, defining the tolerance for enforcing positive definiteness. Default value is 1e-8.
#' @param maxit maximum number of iterations. Default value is 100.
#' @param trace a logical, indicating whether iterations are to be traced (printed out). Default value is FALSE.
#' @return a positive-definite matrix.
#' @details Note that setting corr = TRUE just sets diag(.) <- 1 within the algorithm. The near_posdef() is originally found under \insertCite{matrixpackage;textual}{CovEstim}.
#' @examples
#'
#' data(sp200)
#' sp_rets <- sp200[1:100,-1]
#' ml_sigma <- sigma_estim(sp_rets, "ML")
#' ml_sigma_near_posdef <- near_posdef(ml_sigma, eig_tol=1e-8)
#'
#' @importFrom Rdpack reprompt
#' @references
#'\insertAllCited
#'
#'
#' @export near_posdef
#'
near_posdef <-
  function(mat,
           corr = FALSE,
           keep_diag = FALSE,
           do_dykstra = TRUE,
           eig_tol = 1e-6,
           conv_tol = 1e-7,
           posd_tol = 1e-8,
           maxit = 100,
           trace = FALSE) {
    n <- dim(mat)[2]
    if (keep_diag) {
      diag_mat0 <- diag(mat)
    }
    if (do_dykstra) {
      ds_mat <- matrix(0, n, n)
    }

    iter <- 0
    converged <- FALSE
    conv <- Inf

    while (iter < maxit && !converged) {
      mat_y <- mat
      if (do_dykstra) {
        r_mat <- mat_y - ds_mat
      }

      ## project onto PSD matrices  X_k  =  P_S (R_k)
      eigen_tmp <-
        eigen(if (do_dykstra)
          r_mat
          else
            mat_y, symmetric = TRUE)
      eigenvec <- eigen_tmp$vectors
      eigenval <- eigen_tmp$values

      ## create mask from relative positive eigenvalues
      eigenval_log <-  eigenval > eig_tol * eigenval[1]
      if (!any(eigenval_log))
        stop("Matrix seems negative semi-definite")

      ## use p mask to only compute 'positive' part
      eigenvec <- eigenvec[, eigenval_log, drop = FALSE]
      mat <-
        tcrossprod(eigenvec * rep(eigenval[eigenval_log], each = nrow(eigenvec)), eigenvec)

      if (do_dykstra) {
        ## update Dykstra's correction D_S = \Delta S_k
        ds_mat <- mat - r_mat
      }

      if (corr) {
        diag(mat) <- 1
      } else if (keep_diag) {
        diag(mat) <- diag_mat0
      }

      conv <- norm(mat_y - mat, "I") / norm(mat_y, "I")
      iter <- iter + 1

      if (trace) {
        cat(sprintf(
          "iter %3d : #{p}=%d, ||Y-X|| / ||Y||= %11g\n",
          iter,
          sum(eigenval),
          conv
        ))
      }

      converged <- (conv <= conv_tol)
    }

    if (!converged) {
      warning(gettextf("'nearPosDef()' did not converge in %d iterations", iter),
              domain = NA)
    }

    eigen_tmp <- eigen(mat, symmetric = TRUE) #e
    eigenval <- eigen_tmp$values
    eps <- posd_tol * abs(eigenval[1])
    if (eigenval[n] < eps) {
      eigenval[eigenval < eps] <- eps
      eigenvec <- eigen_tmp$vectors
      o_diag <- diag(mat)
      mat <- eigenvec %*% (eigenval * t(eigenvec))
      d_mat <- sqrt(pmax(eps, o_diag) / diag(mat))
      mat[] <- d_mat * mat * rep(d_mat, each = n)
    }

    if (corr) {
      diag(mat) <- 1
    } else if (keep_diag) {
      diag(mat) <- diag_mat0
    }
    sigma_mat <- as.matrix(mat)

    return(sigma_mat)
  }


#' Squared Root of a Matrix
#'
#' Calculates the squared root of a matrix.
#'
#' @param mat a positive-definite pxp matrix.
#' @return a pxp matrix, the squared root of mat.
#'
#' @details The squared root of a positive-definite symmetric matrix X is calculated with the following formula:
#' \deqn{\sqrt{X}=\Delta \sqrt{\Lambda}\Delta^{-1},}
#' where \eqn{\Delta} is the matrix of eigenvectors and \eqn{\sqrt{\Lambda}} is a diagonal matrix with the squared roots of the eigenvalues of X.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_ml <- sigma_estim(sp_rets, "ML")
#' sigma_ml_sqrt_root <- sqrt_root_calc(sigma_ml)
#'
#' @export sqrt_root_calc
#'
sqrt_root_calc <- function(mat) {
  eigen_tmp <- eigen(mat)
  eigenval <- eigen_tmp$values
  eigenvec <- eigen_tmp$vectors
  sqrt_roor_mat <-
    eigenvec %*% diag(sqrt(eigenval)) %*% solve(eigenvec)

  return(sqrt_roor_mat)
}
