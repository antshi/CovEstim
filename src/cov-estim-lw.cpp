
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Ledoit-Wolf Linear Shrinkage Covariance Estimation I (CPP)
//'
//' Computes the Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards the identity matrix.
//'
//' @param data an nxp data matrix.
//' @param shrink_int a double, indicating the shrinkage intensity. Default is the optimal shrinkage intensity as in \insertCite{ledoit2003identity;textual}{CovEstim}.
//' @param zeromean_log a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE). Default value is FALSE.
//' @return a pxp estimated covariance matrix.
//'
//' @details The Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards the identity matrix is calculated with the following formula:
//' \deqn{\hat{\Sigma}= s\Sigma_{T} + (1-s)\Sigma,}
//' where \eqn{\Sigma} is the sample covariance matrix, s is the user-supplied or optimal shrinkage intensity and \eqn{\Sigma_{T}} is a pxp identity matrix.
//' This covariance estimator assumes a zero correlation and variances of one as the underlying covariance structure of the data.
//' A corresponding MATLAB code for the estimator can be accessed under \url{https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html}.
//'
//' @examples
//' data(sp200)
//' sp_rets <- sp200[,-1]
//' sigma_lwident <- sigma_estim_lwident_cpp(as.matrix(sp_rets))
//'
//' @importFrom Rcpp evalCpp
//' @useDynLib CovEstim
//' @importFrom Rdpack reprompt
//' @references
//'\insertAllCited
//'
//' @export sigma_estim_lwident_cpp
//'
// [[Rcpp::export]]
arma::mat sigma_estim_lwident_cpp(arma::mat data, double shrink_int = -1, bool zeromean_log = false){

	int p = data.n_cols;
	int n = data.n_rows;

	if(zeromean_log==false){

		arma::mat mean_vec = mean(data,0);

		for(int i=0;i<p;i=i+1){
		data.col(i) += - mean_vec(i);
		}
		n = n - 1;

	}

	arma::mat sigma_sample = trans(data)*data/n;
	arma::mat sigma_target;
    sigma_target.eye(p,p);

	if(shrink_int == -1){

		arma::mat asyvarmat = trans(square(data))*(square(data))/n- 2*trans(data)*(data)%sigma_sample/n + square(sigma_sample);

		double asyvar = accu(asyvarmat);

		double gamma = accu(square(sigma_target-sigma_sample));

		double kappa = (asyvar/gamma)/n;

		arma::vec minf = {1,kappa};
		double minval = minf.min();

		arma::vec shrinkf = {0,minval};
		shrink_int = shrinkf.max();
	}

	arma::mat sigma_estim = shrink_int * sigma_target + (1 - shrink_int) * sigma_sample;


	return sigma_estim;
}

//' Ledoit-Wolf Linear Shrinkage Covariance Estimation II (CPP)
//'
//' Computes the Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards the one-parameter matrix.
//'
//' @param data an nxp data matrix.
//' @param shrink_int a double, indicating the shrinkage intensity. Default is the optimal shrinkage intensity as in \insertCite{ledoit2004oneparam;textual}{CovEstim}.
//' @param zeromean_log a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE). Default value is FALSE.
//' @return a pxp estimated covariance matrix.
//'
//' @details The Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards the diagonal matrix of equal variances is calculated with the following formula:
//' \deqn{\hat{\Sigma}= s\Sigma_{T} + (1-s)\Sigma,}
//' where \eqn{\Sigma} is the sample covariance matrix, s is the user-supplied or optimal shrinkage intensity and
//' \eqn{\Sigma_{T}} is a diagonal matrix with the average sample variance \eqn{\bar{\sigma}^2} on the diagonal.
//' This covariance estimator assumes a zero correlation and equal variances as the underlying covariance structure of the data.
//' A corresponding MATLAB code for the estimator can be accessed under \url{https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html}.
//'
//' @examples
//' data(sp200)
//' sp_rets <- sp200[,-1]
//' sigma_lwone <- sigma_estim_lwone_cpp(as.matrix(sp_rets))
//'
//' @importFrom Rcpp evalCpp
//' @useDynLib CovEstim
//' @importFrom Rdpack reprompt
//' @references
//'\insertAllCited
//'
//' @export sigma_estim_lwone_cpp
//'
// [[Rcpp::export]]
arma::mat sigma_estim_lwone_cpp(arma::mat data, double shrink_int = -1, bool zeromean_log = false){

  int p = data.n_cols;
  int n = data.n_rows;

  if(zeromean_log==false){

    arma::mat mean_vec = mean(data,0);

    for(int i=0;i<p;i=i+1){
      data.col(i) += - mean_vec(i);
    }
    n = n - 1;

  }

  arma::mat sigma_sample = trans(data)*data/n;
  arma::mat sigma_target;

  double aver_var = mean(sigma_sample.diag());
  sigma_target.eye(p,p);
  sigma_target = aver_var*sigma_target;

  if(shrink_int == -1){

    double aver_var = mean(sigma_sample.diag());
    sigma_target.eye(p,p);
    sigma_target = aver_var*sigma_target;

    arma::mat asyvarmat = trans(square(data))*(square(data))/n- 2*trans(data)*(data)%sigma_sample/n + square(sigma_sample);
    double asyvar = accu(asyvarmat);

    double gamma = accu(square(sigma_target-sigma_sample));

    double kappa = (asyvar/gamma)/n;

    arma::vec minf = {1,kappa};
    double minval = minf.min();

    arma::vec shrinkf = {0,minval};
    shrink_int = shrinkf.max();

  }

  arma::mat sigma_estim = shrink_int * sigma_target + (1 - shrink_int) * sigma_sample;


  return sigma_estim;
}

//' Ledoit-Wolf Linear Shrinkage Covariance Estimation III (CPP)
//'
//' Computes the Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards the constant correlation covariance matrix.
//'
//' @param data an nxp data matrix.
//' @param shrink_int a double, indicating the shrinkage intensity. Default is the optimal shrinkage intensity as in \insertCite{ledoit2004cc;textual}{CovEstim}.
//' @param zeromean_log a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE). Default value is FALSE.
//' @return a pxp estimated covariance matrix.
//'
//' @details The Ledoit-Wolf linear shrinkage estimator of the covariance matrix towards the constant correlation covariance matrix is calculated with the following formula:
//' \deqn{\hat{\Sigma}= s\Sigma_{T} + (1-s)\Sigma,}
//' where \eqn{\Sigma} is the sample covariance matrix, s is the user-supplied or optimal shrinkage intensity and
//' \eqn{\Sigma_{T}} is the constant correlation covariance matrix.
//' This covariance estimator assumes a constant correlation and individual variances as the underlying covariance structure of the data.
//' A corresponding MATLAB code for the estimator can be accessed under \url{https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html}.
//'
//' @examples
//' data(sp200)
//' sp_rets <- sp200[,-1]
//' sigma_lwcc <- sigma_estim_lwcc_cpp(as.matrix(sp_rets))
//'
//' @importFrom Rcpp evalCpp
//' @useDynLib CovEstim
//' @importFrom Rdpack reprompt
//' @references
//'\insertAllCited
//'
//' @export sigma_estim_lwcc_cpp
//'
// [[Rcpp::export]]
arma::mat sigma_estim_lwcc_cpp(arma::mat data, double shrink_int = -1, bool zeromean_log = false){

  int p = data.n_cols;
  int n = data.n_rows;
  int obs = data.n_rows;

  if(zeromean_log==false){

    arma::mat mean_vec = mean(data,0);

    for(int i=0;i<p;i=i+1){
      data.col(i) += - mean_vec(i);
    }
    n = n - 1;

  }

  arma::mat sigma_sample = trans(data)*data/n;
  arma::mat sigma_target;

  arma::mat Corr = cor(data, data);
  Corr.diag().zeros();
  double rbar = mean(mean(Corr))*p/(p-1);

  sigma_target = rbar * sigma_sample/Corr;
  sigma_target.diag() = diagvec(sigma_sample);

  if(shrink_int == -1){

    arma::mat asyvarmat = trans(square(data))*(square(data))/n- 2*trans(data)*(data)%sigma_sample/n + square(sigma_sample);
    double asyvar = accu(asyvarmat);

    arma::mat term1 = trans(pow(data, 3))*data;

    arma::mat term2 = trans(data)*data;
    arma::vec vars = diagvec(sigma_sample);

    for(int i=0;i<p;i=i+1){
      term2.row(i) *= vars(i);
    }

    arma::mat onesmat = ones(obs, p);
    arma::mat term3 = sigma_sample%(trans(square(data))*onesmat);

    arma::mat term4 = diagvec(sigma_sample)%(ones(p));
    term4 = repmat(term4, 1, p)%sigma_sample;

    arma::mat terms = (term1 - term2 - term3 + term4)/n;

    arma::mat ratios = pow(diagvec(sigma_sample)*trans(pow(diagvec(sigma_sample), -1)), 0.5);

    arma::mat rhos = 0.5*rbar*(ratios%terms + trans(ratios)%trans(terms));
    rhos.diag().zeros();

    double rho = accu(asyvarmat.diag()) + accu(rhos);

    double gamma = accu(square(sigma_target-sigma_sample));

    double kappa = ((asyvar-rho)/gamma)/n;

    arma::vec minf = {1,kappa};
    double minval = minf.min();

    arma::vec shrinkf = {0,minval};
    shrink_int = shrinkf.max();

  }

  arma::mat sigma_estim = shrink_int * sigma_target + (1 - shrink_int) * sigma_sample;


  return sigma_estim;
}

arma::mat pmax(arma::mat K, double comp_num){ /*replaces pmax*/
  int comp_p = K.n_cols;
  int comp_n = K.n_rows;
  arma::vec compvals;
  for(int i=0;i<comp_n;i++){
    for(int j=0;j<comp_p;j++){
      compvals = {comp_num,K(i,j)};
      K(i,j) = compvals.max();
    }
  }
  return K;
}

//' Ledoit-Wolf Covariance Estimation (Nonlinear Shrinkage) (CPP)
//'
//' Computes the analytical Ledoit-Wolf nonlinear shrinkage estimator of the covariance matrix.
//'
//' @param data an nxp data matrix.
//' @param bandwidth_speed a double, indicating the speed at which the bandwidth vanishes in the number of variables p.
//' Default value is -0.35.
//' @param zeromean_log a logical, indicating whether the data matrix has zero means (TRUE) or not (FALSE). Default value is FALSE.
//' @return a pxp estimated covariance matrix.
//'
//' @details The Ledoit-Wolf nonlinear shrinkage estimator of the covariance matrix is computed according to \insertCite{ledoit2018analytical;textual}{CovEstim}
//' with the following formula:
//' \deqn{\hat{\Sigma}=\Delta\hat{\Lambda}\Delta',}
//' where \eqn{\Delta} is the matrix with the sample eigenvectors of the data matrix and \eqn{\hat{\Lambda}} is a diagonal matrix with the sample eigenvalues, shrunk in a nonlinear way.
//' The optimal solution is achieved using a nonparametric variable bandwidth kernel estimation of the limiting spectral density of the sample eigenvalues and its Hilbert transform.
//' The speed at which the bandwidth vanishes in the number of assets is set to -0.35.
//' A corresponding MATLAB code for the estimator can be accessed under \url{https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html}.
//'
//' @examples
//' data(sp200)
//' sp_rets <- sp200[,-1]
//' sigma_lwnl <- sigma_estim_lwnl_cpp(as.matrix(sp_rets))
//'
//' @import Rcpp
//' @importFrom Rcpp evalCpp
//' @useDynLib CovEstim
//' @importFrom Rdpack reprompt
//' @references
//'\insertAllCited
//'
//' @export sigma_estim_lwnl_cpp
//'
// [[Rcpp::export]]
arma::mat sigma_estim_lwnl_cpp(arma::mat data, double bandwidth_speed=-1, bool zeromean_log = false){

  double p = data.n_cols;
  double n = data.n_rows;

  if(zeromean_log==false){

    arma::mat mean_vec = mean(data,0);

    for(int i=0;i<p;i=i+1){
      data.col(i) += - mean_vec(i);
    }
    n = n - 1;

  }

  arma::mat Sigma = trans(data)*data/n;
  arma::vec eigval;
  arma::mat eigvec;

  eig_sym(eigval, eigvec, Sigma);

  double comp_val = p-n;
  arma::vec maxf = {0,comp_val};
  double start_ind = maxf.max();

  arma::vec lambda = eigval(span(start_ind,p-1));

  arma::vec p_n = {p,n};
  double minpn = p_n.min();

  arma::mat Ident(1,minpn);
  Ident.fill(1);
  arma::mat L = kron(Ident,lambda);

  if(bandwidth_speed==-1){

    bandwidth_speed = -0.35;
  }
  double h = pow(n,bandwidth_speed);

  arma::mat H = h*trans(L);
  arma::mat x = (L - trans(L))/H;

  double sqrtfive = pow(5, 0.5);
  double pi_val = datum::pi;

  double express_five = double(3)/double(4)/sqrtfive;
  double express_pi = double(-3)/double(10)/pi_val;

  arma::mat help_compare = 1 - square(x)/5;
  arma::mat compare = pmax(help_compare, 0);
  arma::mat ftildemat = compare/H;
  arma::vec ftilde = express_five*mean(ftildemat,1);

  arma::mat Hftemp = express_pi*x + (express_five/pi_val)*help_compare%(log(abs((sqrtfive - x)/(sqrtfive + x))));

  arma::vec equalelem = x.elem(find(abs(x)==sqrtfive));
  Hftemp.elem(find(abs(x)==sqrtfive)) = express_pi*equalelem;
  arma::mat Hfmat = Hftemp/H;
  arma::vec Hftilde = mean(Hfmat, 1);

  arma::vec dtilde(p);
  arma::vec dtilde1;

  double Hftilde0;
  double dtilde0;

  if (p <= n) {
    dtilde = lambda/((square(pi_val*(p/n)*lambda%ftilde))+(square(1-(p/n)-pi_val*(p/n)*lambda%Hftilde)));
  } else {
    arma::vec revlambda = 1/lambda;
    Hftilde0 = (double(1)/pi_val)*((double(3)/double(10)/pow(h,2)) + (express_five/h)*(1-double(1)/double(5)/pow(h,2))*log((1+sqrtfive*h)/(1-sqrtfive*h)))*mean(revlambda);
    dtilde0 = 1/(pi_val*((p-n)/n)*Hftilde0);
    dtilde1 = lambda/(pow(pi_val,2)*square(lambda)%((square(ftilde))+square(Hftilde)));
    dtilde.fill(dtilde0);
    dtilde(span(p-n,dtilde.n_elem-1)) = dtilde1;
  }

  arma::vec eigenval_lwnl = dtilde;
  arma::mat sigma_estim = eigvec*diagmat(eigenval_lwnl)*trans(eigvec);

  return sigma_estim;
}

/*
 vec iso_reg(vec y){

 int n = y.n_elem, i, ip, known, n_ip;
 double tmp, slope;
 vec yc(n+1), yf(n);

 yc(0) = 0;
 tmp = 0;
 for (i = 0; i < n; i++) {
 tmp += y(i);
 yc(i + 1) = tmp;
 }
 known = 0; ip = 0, n_ip = 0;
 do {
 slope = 10000000;
 for (i = known + 1; i <= n; i++) {
 tmp = (yc(i) - yc(known)) / (i - known);
 if (tmp < slope) {
 slope = tmp;
 ip = i;
 }
 }

 for (i = known; i < ip; i++)
 yf(i) = (yc(ip) - yc(known)) / (ip - known);
 } while ((known = ip) < n);

 return yf;

 }
 */
