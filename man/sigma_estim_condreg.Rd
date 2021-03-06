% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cov-estim-condreg.R
\name{sigma_estim_condreg}
\alias{sigma_estim_condreg}
\title{Condition-Number Regularized Covariance Estimation}
\usage{
sigma_estim_condreg(
  data,
  k = NULL,
  k_seq = NULL,
  k_seq_len = 50,
  nfolds = 5,
  zeromean_log = FALSE
)
}
\arguments{
\item{data}{an nxp data matrix.}

\item{k}{a double, indicating the regularization parameter (or the maximum condition number for the estimated covariance matrix).
Default value is NULL and the optimal k is found with a cross-validation (CV), optimizing the negative likelihood.}

\item{k_seq}{a vector of doubles, specifying the grid of k values to search over within the CV.
Default value is NULL and the sequence is generated in dependence of the sample covariance matrix, the user-supplied length and the minimum deviation ratio of the values.}

\item{k_seq_len}{an integer, indicating the length of k_seq. Default value is 50.}

\item{nfolds}{an integer, specifying the number of folds for the CV. Default value is 5.}

\item{zeromean_log}{a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE). Default value is FALSE.}
}
\value{
a list with the following entries
\itemize{
\item a pxp estimated covariance matrix.
\item an estimation specific tuning parameter, here the regularization parameter k.
}
}
\description{
Computes the condition-number regularized (CONDREG) estimator of the covariance matrix.
}
\details{
The CONDREG estimator is elaborated in detail in \insertCite{won2013condition;textual}{CovEstim}. More information on the functionality can be found in
\insertCite{condregpackage;textual}{CovEstim}.
}
\examples{
data(sp200)
sp_rets <- sp200[,-1]
results_condreg <- sigma_estim_condreg(sp_rets)
sigma_condreg <- results_condreg[[1]]
param_condreg <- results_condreg[[2]]
sigma_condreg <- sigma_estim_condreg(sp_rets, k=100)[[1]]

}
\references{
\insertAllCited
}
