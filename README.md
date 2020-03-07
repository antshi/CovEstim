
# CovEstim

<!-- badges: start -->
<!-- badges: end -->

The package CovEstim includes various functions for estimating the covariance matrix of a dataset, especially in a high-dimensional setting. 
A dataset of 200 monthly returns of stocks from the S&P 500 index is provided for testing purposes.

Here is a short overview of the estimation methods:

* Sample
* Maximum-Likelihood
* Bayes-Stein
* Stein-Haff
* Eigenvalue Clipping
  * Marcenko-Pastur Edge
  * Bouchaud-Potters
* Linear Shrinkage
  * target with constant variance equal to 1 and correlation equal to 0
  * target with constant variance equal to the average of variances and correlation equal to 0
  * target with constant correlation equal to the average of correlations
* Nonlinear Shrinkage
  * Ledoit-Wolf
  * NERCOME
* Factor Models
  * Exact Factor Model
  * Single Market Factor Model (for stock returns)
  * Approximate Factor Model
* PCA
* POET 
* Ridge Regularization
* Condition Number Regularization
* Graphical Lasso Regularization
* t-distributed Lasso Regularization
* Graphical Horseshoe 
* EWMA


## Installation

You can install this version of the package from Github with:

``` r
install.packages("devtools")
library(devtools)
install_github("antshi/CovEstim")
library(CovEstim)
```

## 

These are basic examples which show you how to use the wrapper function sigma_estim_wrapper. First, let's prepare with

``` r
# include the package
library(covEstim)

# data
data(sp200)
str(sp200)

# delete the date column, since the estimation functions expect a matrix.
sp_rets <- as.matrix(sp200[,-1])
```

### Example 1

Maximum-Likelihood Estimator
```r
sigma_ml <- sigma_estim_ml(sp_rets)[[1]]

sigma_ml <- sigma_estim_wrapper(sp_rets, estim_func=sigma_estim_ml)

sigma_ml <- sigma_estim(sp_rets, est_type="ML")

```

### Example 2

Linear Shrinkage towards the Constant Correlation Covariance Matrix
```r
sigma_lwcc_results <- sigma_estim_lwcc(sp_rets)
sigma_lwcc <- sigma_lwcc_results[[1]]
sigma_lwcc_param <- sigma_lwcc_results[[2]]

# The argument res_all=FALSE (default) in the wrapper functions sigma_estim and sigma_estim_wrapper returns only the covariance estimator.
sigma_lwcc <- sigma_estim_wrapper(sp_rets, sigma_estim_lwcc, res_all=FALSE)
```
You can also parse the user-defined tuning parameters such as shrinking intensity to the estimation function. 

```r
sigma_lwcc_results <- sigma_estim_wrapper(sp_rets, sigma_estim_lwcc, res_all=TRUE, shrink_int=0.3)
str(sigma_lwcc_results)
```

### Example 3

The linear and nonlinear shrinkage estimators by Ledoit and Wolf are implemented additionally in Rcpp for higher efficiency.
An important note here is that while the R-function can deal with data frames and matrices, 
the Rcpp-function requires the data to be preformated as a matrix.

```r
sigma_lwcc_results <- sigma_estim_wrapper(sp_rets, sigma_estim_lwcc, res_all=TRUE)
str(sigma_lwcc_results)

sigma_lwcc_results <- sigma_estim_wrapper(sp_rets, sigma_estim_lwcc_cpp, res_all=TRUE)
str(sigma_lwcc_results)
```
### Example 4

Some more complicated application would be the Approximate Factor Model with a nonlinear shrinkage estimation on the covariance matrix of residuals

```r
sigma_afm_lwnl_results <- sigma_estim_wrapper(sp_rets, estim_func = sigma_estim_afm, resid_estim_func=sigma_estim_lwnl)
str(sigma_afm_lwnl_results)
```
### Example 5

While the function sigma_estim_wrapper can take up any estimation function from the package, the function sigma_estim is interesting for practical applications.
For example, if you want to calculate different covariance estimators on the same dataset for comparison purposes, sigma_estim is easier to work with.

```r
# Maximum-Likelihood, Exact Factor Model, Ledoit-Wolf linear shrinkage and Approximate Factor Model with nonlinear shrinkage on the residuals:
est_types <- c("ML", "EFM", "LW-IDENT", "AFM-LWNL") 
sigma_list <- lapply(est_types, sigma_estim, data=sp_rets)
```



