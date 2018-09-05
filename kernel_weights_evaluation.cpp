#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double phi(double x){
  double kernel = 1;
  if (abs(x) < 1){
    kernel = exp( - 1/(1 - x * x));
  }
  else{
    kernel = 0;
  }
  return(kernel);
}

// [[Rcpp::export]]
double log_lik_1(double y, double x, NumericVector mu_star,
               NumericVector mu, NumericMatrix omega, 
               double h_K, double sigma){
  int K = omega(0, _).size();
  int m = omega(_, 0).size();
  NumericVector w(K);
  double phi(double x);
  for (int k = 0; k < K; k++){
    w(k) = phi((x - mu[k])/h_K);
    // Rcout<<w(k)<<std::endl;
  }
  double D = sum(w);
  if (D > 0){ w = w/D; }
  NumericVector regression(K);
  for (int k = 0; k < K; k++){
    regression[k] = 0;
    for (int alpha = 0; alpha < m; alpha ++){
      regression[k] = regression[k] + w[k] * omega(alpha, k) * pow(x - mu_star[k], alpha);
    }
    // Rcout<<regression(k)<<std::endl;
  }
  // Rcout<<(y - sum(regression))<<std::endl;
  // Rcout<<(y - sum(regression)) * (y - sum(regression))<<std::endl;
  // Rcout<<-1/(2 * sigma * sigma) * (y - sum(regression)) * (y - sum(regression))<<std::endl;
  return(-1/(2 * sigma * sigma) * (y - sum(regression)) * (y - sum(regression)));
}

// [[Rcpp::export]]
double log_lik_joint(NumericVector y, NumericVector x, NumericVector mu_star,
                 NumericVector mu, NumericMatrix omega, 
                 double h_K, double sigma){
  int n = y.size();
  double log_lik_1(double y, double x, NumericVector mu_star,
                   NumericVector mu, NumericMatrix omega, 
                   double h_K, double sigma);
  NumericVector log_lik_tmp(n);
  for (int i = 0; i < n; i++){
    log_lik_tmp(i) = log_lik_1(y[i], x[i], mu_star, mu, omega, h_K, sigma);
  }
  return(sum(log_lik_tmp));
}

// [[Rcpp::export]]
NumericMatrix W(NumericVector x, NumericVector mu, double h_K){
  int K = mu.size();
  int n = x.size();
  NumericMatrix W_mat(n, K);
  for (int i = 0; i < n; i++){
    NumericVector w(K);
    for (int k = 0; k < K; k++){
      w(k) = phi((x[i] - mu[k])/(h_K));
    }
    double D = sum(w);
    if (D > 0){ w = w/D; }
    W_mat(i, _) = w;
  }
  return(W_mat);
}

// [[Rcpp::export]]
NumericMatrix Design_mat_alpha(NumericVector x, NumericVector mu_star, NumericMatrix W, int alpha){
  int K = mu_star.size();
  int n = W(_, 0).size();
  NumericMatrix W_new(n, K);
  for (int i = 0; i < n; i++){
    for (int k = 0; k < K; k++){
      W_new(i, k) = W(i, k) * pow(x[i] - mu_star[k], alpha - 1);
    }
  }
  return(W_new);
}

// // [[Rcpp::export]]
// double Gibbs_type_prior(NumericMatrix mu, double g0){
//   int K = mu(0, _).size();
//   double log_Psi;
//   log_Psi = 0;
//   for (int k1 = 0; k1 < K - 1; k1++){
//     for (int k2 = k1 + 1; k2 < K; k2++){
//       double d;
//       d = sqrt(sum((mu(_, k1) - mu(_, k2)) * (mu(_, k1) - mu(_, k2))));
//       log_Psi = log_Psi + log(d/(g0 + d))/K;
//     }
//   }
//   return(log_Psi);
// }

// [[Rcpp::export]]
NumericVector f_mean_p1(NumericVector x, NumericVector mu_star, 
                        NumericVector mu, NumericMatrix omega, double h_K){
  int n = x.size();
  int K = mu.size();
  int m = omega(_, 0).size();
  NumericVector f(n);
  for (int i = 0; i < n; i ++){
    f(i) = 0;
    double phi(double x);
    NumericVector w(K);
    for (int k = 0; k < K; k ++){
      w(k) = phi((x[i] - mu[k])/h_K);
    }
    double D = sum(w);
    if (D > 0){ w = w/D; }
    NumericVector regression(K);
    for (int k = 0; k < K; k++){
      regression[k] = 0;
      for (int alpha = 0; alpha < m; alpha ++){
        regression[k] = regression[k] + w[k] * omega(alpha, k) * pow(x[i] - mu_star[k], alpha);
      }
    }
    f(i) = sum(regression);
  }
  return(f);
}

// // [[Rcpp::export]]
// NumericVector volterra_function(NumericVector x, int N){
//   int n = x.size();
//   NumericVector f0(n);
//   for (int i = 0; i < n; i ++){
//     NumericVector fx(N);
//     for (int j = 0; j < N; j++){
//       fx[j] = 1/sqrt((j + 1) * (j + 1) * (j + 1)) * sin(j + 1) * cos((j + 1/2) * 3.14159265358979323846 * x[i]);
//     }
//     f0[i] = sum(fx);
//   }
//   return(sqrt(2) * f0);
// }

// // [[Rcpp::export]]
// double power_function(double x, double alpha){
//   return(pow(x, alpha));
// }
// 
// /*** R
// power_function(3, 3/2)
// */

