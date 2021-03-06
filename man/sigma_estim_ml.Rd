% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cov-estim-simple.R
\name{sigma_estim_ml}
\alias{sigma_estim_ml}
\title{Maximum-Likelihood Covariance Estimation}
\usage{
sigma_estim_ml(data)
}
\arguments{
\item{data}{an nxp data matrix.}
}
\value{
a list with the following entries
\itemize{
\item a pxp estimated covariance matrix.
\item an estimation specific tuning parameter, here an NA.
}
}
\description{
Computes the Maximum-Likelihood estimator of the covariance matrix.
}
\details{
The Maximum-Likelihood estimator of the covariance matrix for a data matrix X is computed with the following formula:
\deqn{\hat{\Sigma}=\frac{1}{n} \left(X - \widehat{{\mu}} {1} \right)' \left({X} -  \widehat{{\mu}}{1}\right)}
where \eqn{\mu=\bar{x}_{j}=\frac{1}{n}\sum_{i=1}^{n}x_{ij}} for (for \eqn{i=1,\ldots, n} and \eqn{j=1,\ldots,p}) is the sample mean vector and \eqn{1} is an 1xp vector of ones.
}
\examples{
data(sp200)
sp_rets <- sp200[,-1]
sigma_ml <- sigma_estim_ml(sp_rets)[[1]]

}
