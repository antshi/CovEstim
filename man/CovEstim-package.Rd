\name{CovEstim-package}
\alias{CovEstim}
\docType{package}
\title{
\packageTitle{CovEstim}
}
\description{
\packageDescription{CovEstim}
}

\author{
\packageAuthor{CovEstim}
}

\examples{
data(sp200)
sp_rets <- sp200[,-1]
sigma_ml <- sigma_estim_wrapper(sp_rets, estim_func=sigma_estim_ml) # ML estimator
sigma_lwnl_cpp <- sigma_estim_lwnl_cpp(as.matrix(sp_rets)) #LW-Nonlinear Shrinkage with Rcpp
}
