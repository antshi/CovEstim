#' #' Ledoit-Wolf Linear Shrinkage Covariance Estimation IV
#'
#' Computes the Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards the one-factor covariance matrix.
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
#' @details The Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards the one-factor covariance matrix is calculated with the following formula:
#' \deqn{\hat{\Sigma}= s\Sigma_{T} + (1-s)\Sigma,}
#' where \eqn{\Sigma} is the sample covariance matrix, s is the user-supplied or optimal shrinkage intensity and \eqn{\Sigma_{T}} is the covariance matrix estimator, given by a one-factor model,
#' where the factor is equal to the cross-sectional average of all the variables.
#' This covariance estimator assumes a zero correlation and variances of one as the underlying covariance structure of the data.
#' A corresponding MATLAB code for the estimator can be accessed under \url{https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html}.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_lwcc_sf <- sigma_estim_lwcc_sf(sp_rets)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#'\insertAllCited
#'
#' @export sigma_estim_lwcc_sf
#'
sigma_estim_lwcc_sf <-
  function(data,
           shrink_int = NULL,
           zeromean_log = FALSE) {
    data <- as.matrix(data)
    names_data <- colnames(data)
    p <- dim(data)[2]

    if (!zeromean_log) {
      centered <- apply(data, 2, function(x)
        x - mean(x))
    } else{
      centered <- data
    }
    factors <- rowMeans(centered)

    data_all <- cbind(data, factors)
    rm(data)
    gc()
    if (!zeromean_log) {
      centered_all <- apply(data_all, 2, function(x)
        x - mean(x))
      centered <- centered_all[, 1:p]
      n <- dim(data_all)[1] - 1
    } else{
      centered_all <- apply(data_all, 2, function(x)
        x - mean(x))
      centered <- centered_all[, 1:p]
      n <- dim(data_all)[1]
    }

    sigma_sample_all <- t(centered_all) %*% centered_all / n
    sigma_factors <- sigma_sample_all[1:p, p + 1]
    var_factors <- as.numeric(sigma_sample_all[p + 1, p + 1])
    sigma_sample <- sigma_sample_all[-(p + 1), -(p + 1)]
    sigma_target <- (sigma_factors %o% sigma_factors) / var_factors
    diag(sigma_target) <- diag(sigma_sample)

    if (is.null(shrink_int)) {
      asyvar_mat <-
        (t(centered ^ 2) %*% (centered ^ 2)) / n - sigma_sample ^
        2
      asyvar <- sum(asyvar_mat)
      rhomat_diag <- sum(diag(asyvar_mat))
      term1 <-
        1 / n * t(centered ^ 2) %*% (centered * factors) - sigma_factors * sigma_sample
      term2 <-
        1 / n * t(centered * factors) %*% (centered * factors) - var_factors *
        sigma_sample
      rhomat1 <-
        sum(term1 * t(matrix(
          rep(sigma_factors, p), ncol = p, byrow = FALSE
        )) / var_factors) - sum(diag(term1) * sigma_factors / var_factors)
      rhomat2 <-
        sum(term2 * (sigma_factors %o% sigma_factors) / var_factors ^ 2) - sum(diag(term2) *
                                                                                 (sigma_factors ^ 2) / var_factors ^ 2)
      rho <- 2 * rhomat1 - rhomat2 + rhomat_diag
      gamma <- norm(sigma_target - sigma_sample, "F") ^ 2
      kappa <- (asyvar - rho) / gamma
      shrink_int <- max(0, min(1, kappa / n))
    }
    sigma_mat <-
      shrink_int * sigma_target + (1 - shrink_int) * sigma_sample

    rownames(sigma_mat) <- names_data
    colnames(sigma_mat) <- names_data

    return(list(sigma_mat, shrink_int))

  }

#' Exact Factor Model (EFM) Covariance Estimation
#'
#' Computes the EFM estimator of the covariance matrix.
#'
#' @param data an nxp data matrix.
#' @param factors a nxf matrix with factors. Default value is NULL and the factor is equal to the cross-sectional average of all the variables in the data.
#' @param zeromean_log a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE). Default value is FALSE.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here NA.
#' }
#'
#' @details The EFM covariance estimator is calculated with the following formula:
#' \deqn{\hat{\Sigma} = B\hat{\Sigma}_F B' + \hat{\Sigma_u},}
#' where \eqn{\hat{\Sigma}_F} is the sample covariance matrix of the common factors and \eqn{\hat{\Sigma}_u} is the covariaance matrix of residuals, assumed to have zero correlation.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_efm <- sigma_estim_efm(sp_rets)[[1]]
#'
#'
#' @export sigma_estim_efm
#'
sigma_estim_efm <- function(data,
                            factors = NULL) {

  data <- as.matrix(data)
  names_data <- colnames(data)
  if (is.null(factors)) {
      factors <- as.matrix(rowMeans(data))
  }
  sigma_factors <-
    stats::var(factors, use = "pairwise", na.rm = TRUE)
  fm <- stats::lm(data ~ factors)
  rm(data, factors)
  gc()
  factor_betas <- fm$coefficients[2,]
  sigma_res <- diag(diag(stats::var(fm$residuals)))
  sigma_fm <-
    as.matrix(factor_betas %*% sigma_factors %*% t(factor_betas))
  sigma_mat <-  as.matrix(sigma_fm + sigma_res)

  colnames(sigma_mat) <- names_data
  rownames(sigma_mat) <- names_data

  return(list(sigma_mat, NA))
}

#' Approximate Factor Model (AFM) Covariance Estimation
#'
#' Computes the AFM estimator of the covariance matrix.
#'
#' @param data an nxp data matrix.
#' @param factors a nxf matrix with factors. Default value is NULL and the factor is equal to the cross-sectional average of all the variables in the data.
#' @param zeromean_log a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE). Default value is FALSE.
#' @param resid_est_func a covariance estimation function, applied to the residuals covariance matrix.
#' @param ... further arguments to be parsed to resid_est_func
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, depending on the estimation function type.
#' }
#'
#' @details The AFM covariance estimator is calculated with the following formula:
#' \deqn{\hat{\Sigma} = B\hat{\Sigma}_F B' + \hat{\Sigma_u},}
#' where \eqn{\hat{\Sigma}_F} is the sample covariance matrix of the common factors and \eqn{\hat{\Sigma}_u} is the residuals covariance matrix, estimated with the user-sapplied estim_func.
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_afm <- sigma_estim_afm(sp_rets, resid_est_func=sigma_estim_lwnl)[[1]]
#' results_afm <- sigma_estim_afm(sp_rets, resid_est_func=sigma_estim_lwone, shrink_int=0.1)
#' sigma_afm <- results_afm[[1]]
#' param_afm <- results_afm[[2]]
#'
#' @export sigma_estim_afm
#'
sigma_estim_afm <-
  function(data,
           factors = NULL,
           zeromean_log = FALSE,
           resid_est_func,
           ...) {
    data <- as.matrix(data)
    names_data <- colnames(data)
    if (is.null(factors)) {
      if (!zeromean_log) {
        factors <- as.matrix(rowMeans(apply(data, 2, function(x)
          x - mean(x))))
      } else{
        factors <- as.matrix(rowMeans(data))
      }
    }

    sigma_factors <-
      stats::var(factors, use = "pairwise", na.rm = TRUE)
    fm <- stats::lm(data ~ factors - 1)
    rm(data, factors)
    gc()
    sigma_fm <-
      t(fm$coefficients) %*% sigma_factors %*% t(t(fm$coefficients))
    res_estim <-
      sigma_estim_wrapper(fm$residuals, res_all = TRUE, resid_est_func, ...)
    sigma_res <- res_estim[[1]]
    param_res <- res_estim[[2]]

    sigma_mat <- as.matrix(sigma_fm + sigma_res)

    rownames(sigma_mat) <- names_data
    colnames(sigma_mat) <- names_data

    return(list(sigma_mat, param_res))

  }

#' Covariance Estimation - Preconditioned
#'
#' Computes a specific estimator of the covariance matrix after preconditioning
#' the data with a single factor model - the implicit market portfolio.
#'
#' @param data an nxp data matrix.
#' @param zeromean_log a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE). Default value is FALSE.
#' @param precond_est_func a function for estimating the precondtioned covariance matrix.
#' @param ... further arguments to be parsed to precond_est_func.
#'
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific tuning parameter, here the bandwidth speed.
#' }
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_lwnl_sf <- sigma_estim_precond(sp_rets,
#' precond_est_func=sigma_estim_lwnl, bandwidth_speed=NULL)[[1]]
#' sigma_glasso_sf <- sigma_estim_precond(sp_rets,
#' precond_est_func=sigma_estim_glasso, rho=0.01)[[1]]
#'
#' @importFrom Rdpack reprompt
#' @references
#'\insertAllCited
#'
#' @export sigma_estim_precond
#'
sigma_estim_precond <-
  function(data,
           zeromean_log = FALSE,
           precond_est_func = NULL,
           ...) {
    data <- as.matrix(data)
    names_data <- colnames(data)

    sqrt_sigma_mat_fm <-
      sqrt_root_calc(sigma_estim_efm(data, zeromean_log = zeromean_log)[[1]])
    data_precond <- data %*% solve(sqrt_sigma_mat_fm)
    rm(data)
    gc()
    results_precond <-
      sigma_estim_wrapper(data = data_precond,
                          est_func = precond_est_func,
                          res_all = TRUE,
                          ...)
    sigma_mat_precond <- results_precond[[1]]
    param_precond <- results_precond[[2]]

    sigma_mat <-
      sqrt_sigma_mat_fm %*% sigma_mat_precond %*% sqrt_sigma_mat_fm

    rownames(sigma_mat) <- names_data
    colnames(sigma_mat) <- names_data

    return(list(sigma_mat, param_precond))
  }
