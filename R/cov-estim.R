#' Wrapper Function for Covariance Estimation I
#'
#' Estimates the covariance matrix of a dataset
#' according to the user-defined function.
#'
#' @param data an nxp data matrix.
#' @param res_all a logical, defining the return object.
#' If FALSE, only the estimated covariance matrix is returned.
#' If TRUE, a list with two entries is returned. The first entry is the estimated covariance matrix.
#' The second entry is the estimation specific (tuning) parameter. Default value is FALSE.
#' @param estim_func an estimation function.
#' @param ... additional parameters, parsed to estim_func.
#' @return a list with the following entries
#' \itemize{
#' \item a pxp estimated covariance matrix.
#' \item an estimation specific (tuning) parameter.
#' }
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_ml <- sigma_estim_wrapper(sp_rets, estim_func=sigma_estim_ml)
#'
#'
#' @export sigma_estim_wrapper
#'
sigma_estim_wrapper <-
  function(data, estim_func, res_all = FALSE, ...) {
    result <- estim_func(data, ...)

    if (res_all) {
      return(result)
    } else {
      return(result[[1]])
    }

  }

#' Wrapper Function for Covariance Estimation II
#'
#' Estimates the covariance matrix of a dataset
#' according to the user-defined character string.
#'
#' @param data an nxp data matrix.
#' @param est_type a character string, defining the estimation method.
#' @param param a double, setting an estimation specific (tuning) parameter.
#' @param factors a nxm data matrix of factors.
#' @param res_all a logical, defining the return object.
#' If FALSE, only the estimated covariance matrix is returned.
#' If TRUE, a list with two entries is returned. The first entry is the estimated covariance matrix.
#' The second entry is the estimation specific (tuning) parameter. Default value is FALSE.
#'
#' @return a list with following entries
#' \itemize{
#' \item a pxp estimated covariance matrix
#' \item an estimation specific tuning parameter
#' }
#'
#' @examples
#' data(sp200)
#' sp_rets <- sp200[,-1]
#' sigma_ml <- sigma_estim(sp_rets, "ML")
#' sigma_lwcc <- sigma_estim(sp_rets, "LW-CC", param=0.3, res_all=TRUE)[[1]]
#'
#' @export sigma_estim
#'
sigma_estim <-
  function(data,
           est_type,
           param = NULL,
           factors = NULL,
           res_all = FALSE) {
    if (est_type == "SAMPLE") {
      result <- sigma_estim_sample(data)

    } else if (est_type == "ML") {
      result <- sigma_estim_ml(data)

    } else if (est_type == "BS") {
      result <- sigma_estim_bs(data)

    } else if (est_type == "STEIN-HAFF") {
      result <- sigma_estim_sh(data)

    } else if (est_type == "EVC-MP") {
      result <- sigma_estim_evc_mp(data)

    } else if (est_type == "EVC-BP") {
      result <- sigma_estim_evc_bp(data, param)

    } else if (est_type == "PCA") {
      result <- sigma_estim_pca(data, param)

    } else if (est_type == "RIDGE") {
      result <- sigma_estim_ridge(data, param)

    } else if (est_type == "LW-IDENT") {
      result <- sigma_estim_lwident(data, param)

    } else if (est_type == "LW-ONE") {
      result <- sigma_estim_lwone(data, param)

    } else if (est_type == "LW-CC") {
      result <- sigma_estim_lwcc(data, param)

    } else if (est_type == "LW-NONLIN") {
      result <- sigma_estim_lwnl(data, param)

    } else if (est_type == "POET") {
      result <- sigma_estim_poet(data, param)

    } else if (est_type == "GLASSO") {
      result <- sigma_estim_glasso(data, param)

    } else if (est_type == "EWMA") {
      result <- sigma_estim_ewma(data, param)

    } else if (est_type == "EFM") {
      result <- sigma_estim_efm(data, factors)

    } else if (est_type == "LW-CC-SF") {
      result <- sigma_estim_lwcc_sf(data, param)

    } else if (est_type == "AFM+RIDGE") {
      result <- sigma_estim_afm(data, factors, sigma_estim_ridge, param)

    } else if (est_type == "AFM+LWNL") {
      result <-
        sigma_estim_afm(data, factors, sigma_estim_lwnl, param)

    }

    if (res_all) {
      return(result)
    } else {
      return(result[[1]])
    }
  }
