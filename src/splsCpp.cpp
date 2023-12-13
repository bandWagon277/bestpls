// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::export]]
vec Ust(vec b, double eta) {
  int n = b.n_elem;
  vec b_ust(n);
  vec sign_b = sign(b);
  if (eta < 1) {
    vec valb = abs(b) - eta * max(abs(b));
    for (int i = 0; i < n; i++) {
      if (valb(i) >= 0) {
        b_ust(i) = valb(i) * sign_b(i);
      }
    }
  }
  return b_ust;
}

// [[Rcpp::export]]
mat case_1(mat M,vec c, double eps,int maxstep,double eta){
  double dis = 10;
  int i = 1;
  vec c_new = c;

  while (dis > eps && i <= maxstep) {
    // Optimize a for fixed c
    mat U;
    mat V;
    vec s;
    svd(U, s, V, M * c);
    mat a = U.col(0) * V.t();


    // Soft thresholding for optimizing c (assuming lambda2 -> Inf)
    mat ma = M*a;
    vec ma_v = ma.col(0);
    c_new = Ust(ma_v, eta);

    // Calculate discrepancy between a & c
    dis = max(abs(c-c_new));
    c = c_new;
    i++;
  }
  return c;
}
