% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cov-estim-lwnonlin.R
\name{sigma_estim_lwnl}
\alias{sigma_estim_lwnl}
\title{Ledoit-Wolf Covariance Estimation (Nonlinear Shrinkage)}
\usage{
sigma_estim_lwnl(data, bandwidth_speed = NULL, zeromean_log = FALSE)
}
\arguments{
\item{data}{an nxp data matrix.}

\item{bandwidth_speed}{a double, indicating the speed at which the bandwidth vanishes in the number of variables p.
Default value is -1/3.}

\item{zeromean_log}{a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE). Default value is FALSE.}
}
\value{
a list with the following entries
\itemize{
\item a pxp estimated covariance matrix.
\item an estimation specific tuning parameter, here the bandwidth speed.
}
}
\description{
Computes the analytical Ledoit-Wolf nonlinear shrinkage estimator of the covariance matrix.
}
\details{
The Ledoit-Wolf nonlinear shrinkage estimator of the covariance matrix is computed according to \insertCite{ledoit2018analytical;textual}{CovEstim}
with the following formula:
\deqn{\hat{\Sigma}=\Delta\hat{\Lambda}\Delta',}
where \eqn{\Delta} is the matrix with the sample eigenvectors of the data matrix and \eqn{\hat{\Lambda}} is a diagonal matrix with the sample eigenvalues, shrunk in a nonlinear way.
The optimal solution is achieved using a nonparametric variable bandwidth kernel estimation of the limiting spectral density of the sample eigenvalues and its Hilbert transform.
The speed at which the bandwidth vanishes in the number of assets is set to -1/3.
A corresponding MATLAB code for the estimator can be accessed under \url{https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html}.
}
\examples{
data(sp200)
sp_rets <- sp200[,-1]
sigma_lwnl <- sigma_estim_lwnl(sp_rets)[[1]]

}
\references{
\insertAllCited
}
