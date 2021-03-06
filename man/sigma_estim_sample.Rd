% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cov-estim-simple.R
\name{sigma_estim_sample}
\alias{sigma_estim_sample}
\title{Sample Covariance Estimation}
\usage{
sigma_estim_sample(data)
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
Computes the sample estimator of the covariance matrix.
}
\details{
The sample estimator of the covariance matrix for a data matrix X is computed with the following formula:
\deqn{\hat{\Sigma}=\frac{1}{n-1} \left(X - \widehat{{\mu}} {1} \right)' \left({X} -  \widehat{{\mu}}{1}\right)}
where \eqn{\mu=\bar{x}_{j}=\frac{1}{n}\sum_{i=1}^{n}x_{ij}} (for \eqn{i=1,\ldots, n} and \eqn{j=1,\ldots,p}) is the sample mean vector and \eqn{1} is an 1xp vector of ones.
}
\examples{
data(sp200)
sp_rets <- sp200[,-1]
sigma_sample <- sigma_estim_sample(sp_rets)[[1]]

}
