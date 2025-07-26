#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_MATRIX(WX);
  DATA_SPARSE_MATRIX(W);
  DATA_VECTOR(W_eigs);
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(theta);
  PARAMETER(log_sigma);
  PARAMETER(atanh_rho);
  Type sigma = exp(log_sigma);
  Type rho = tanh(atanh_rho); 
  vector<Type> mu = X * beta + WX * theta;
  int n = y.size();
  vector<Type> y_transformed = y - rho * (W * y);
  vector<Type> res = y_transformed - mu;
  Type nll = -dnorm(res, Type(0.0), sigma, true).sum();
  Type log_det = (log(Type(1.0) - rho * W_eigs)).sum();
  Type total_nll = nll - log_det;
  ADREPORT(rho);
  ADREPORT(sigma);
  
  return total_nll;
}