% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cov-estim-nercome.R
\name{sigma_estim_nercome}
\alias{sigma_estim_nercome}
\title{Nonparametric eigenvalue-regularized Precision or Covariance Matrix Estimation (NERCOME)}
\usage{
sigma_estim_nercome(data, split = NULL)
}
\arguments{
\item{data}{an nxp data matrix.}

\item{split}{an integer, indicating the position of the split.
Default value is NULL and an optimal split is calculated with permutation along a sequence as in \insertCite{lam2016nonparametric;textual}{CovEstim}.}
}
\value{
a list with the following entries
\itemize{
\item a pxp estimated covariance matrix.
\item an estimation specific tuning parameter, here the split.
}
}
\description{
Computes the NERCOME estimator of the covariance matrix.
}
\details{
More information on the NERCOME estimator is found in \insertCite{lam2016nonparametric;textual}{CovEstim}.
}
\examples{
data(sp200)
sp_rets <- sp200[,-1]
sigma_nercome <- sigma_estim_nercome(sp_rets)[[1]]

}
\references{
\insertAllCited
}
