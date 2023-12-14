#' @title an improved implementation of constrained simpls
#' @name simSPLS
#' @param X description a n*q matrix
#' @param Y description a n*p matrix
#' @return a list of pls components
#' @examples
#' X=matrix(rnorm(200),10,20)
#' Y=matrix(rnorm(100),5,20)
#' fit = sparsePLS(X,Y,2)
#' @importFrom pls plsr
#' @import Rcpp
#' @import RcppArmadillo
#' @export

#sourceCpp("src/splsCpp.cpp")

sparsePLS_1<-function (x, y, K=3, eta=0.5, kappa = 0.5, eps = 1e-04, maxstep = 100){
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  ip <- c(1:p)
  y <- as.matrix(y)
  q <- ncol(y)
  one <- matrix(1, 1, n)
  mu <- one %*% y/n
  y <- scale(y, drop(mu), FALSE)
  meanx <- drop(one %*% x)/n
  x <- scale(x, meanx, FALSE)
  normx <- sqrt(drop(one %*% (x^2))/(n - 1))
  if (any(normx < .Machine$double.eps)) {
    stop("Some of the columns of the predictor matrix have zero variance.")
  }
  x <- scale(x, FALSE, normx)
  normy <- rep(1,q)
  y <- scale(y, FALSE, normy)

  betahat <- matrix(0, p, q)
  betamat <- list()
  x1 <- x
  y1 <- y

  if (is.null(colnames(x))) {
    xnames <- c(1:p)
  }
  else {
    xnames <- colnames(x)
  }
  new2As <- list()

  for (k in 1:K) {
    Z <- t(x1) %*% y1
    what <- Weight_vec_cal(Z, eta, kappa, eps, maxstep)
    A <- unique(ip[what != 0 | betahat[, 1] != 0])
    new2A <- ip[what != 0 & betahat[, 1] == 0]
    xA <- x[, A, drop = FALSE]
    plsfit <- pls::plsr(y ~ xA, ncomp = min(k, length(A)),
                        method = "simpls", scale = FALSE)
    betahat <- matrix(0, p, q)
    betahat[A, ] <- matrix(coef(plsfit), length(A), q)
    betamat[[k]] <- betahat
    pj <- plsfit$projection
    y1 <- y - x %*% betahat
    new2As[[k]] <- new2A
  }
  return(list("model" = plsfit,"index" = A))
}

QCQR <- function(A,b,alpha){
  svdA <- svd(A)
  c <- t(svdA$u) %*% b
  sigma <- svdA$d
  f<-function(miu){
    nume <- (sigma^2 + miu)^2
    denom <- (sigma^2)*(c^2)
    return(sum(nume/(denom+1e-4))-1)
  }
  miu <- uniroot(f,c(1e-4,1e30))$root
  inv <- diag(sigma/(sigma^2+miu))
  return(svdA$v%*%(inv%*%c))
}

Weight_vec_cal <-function( Z, eta, kappa, eps, maxstep ){
  # initialization
  p <- nrow(Z)
  q <- ncol(Z)
  Znorm1 <- median(abs(Z))
  Z <- Z/Znorm1
  
  # main iterations
  if ( q==1 ){
    # if univariate response
    # soft thresholding(assuming lambda2 -> Inf)
    c <- Ust(Z,eta)
  }

  if ( q>1 ){
    # if multivariate response
    M <- Z %*% t(Z)
    dis <- 10
    i <- 1

    # main iteration: optimize c and a iteratively
    # use svd solution if kappa==0.5
    if (kappa==0.5){
      # initial value for c (outside the unit circle)
      c <- matrix(10,p,1)
      c <- case_1(M,c,eps,maxstep,eta)
      }

    # solve equation if 0<kappa<0.5
    if (kappa>0 & kappa<0.5){
      kappa2 <- (1-kappa)/(1-2*kappa )
      # initial value for c (outside the unit circle)
      c <- matrix(255,p,1)
      c.old <- c

      #solve by QCQR method
      while (dis>eps & i<=maxstep){
        # optimize a for fixed c
        a <- QCQR(t(Z),kappa2*(t(Z)%*%c),1)

        # optimize c for fixed a
        # soft thresholding (assuming lambda2 -> Inf)
        c <- Ust(M%*%a,eta)

        # calculate discrepancy between a & c
        dis <- max(abs(c-c.old))
        c.old <- c
        i <- i+1
      }
    }
  }
  return(c)
}

cppFunction('arma::vec Ust(arma::vec b, double eta) {
   int n = b.n_elem;
   arma::vec b_ust(n);
   arma::vec sign_b = sign(b);
   if (eta < 1) {
     arma::vec valb = abs(b) - eta * max(abs(b));
     for (int i = 0; i < n; i++) {
       if (valb(i) >= 0) {
         b_ust(i) = valb(i) * sign_b(i);
       }
     }
   }
   return b_ust;
 }',depends = "RcppArmadillo")

cppFunction('arma::mat case_1(arma::mat M,arma::vec c, double eps,int maxstep,double eta){
  double dis = 10;
  int i = 1;
  arma::vec c_new = c;

  while (dis > eps && i <= maxstep) {
    // Optimize a for fixed c
    arma::mat U;
    arma::mat V;
    arma::vec s;
    svd(U, s, V, M * c);
    arma::mat a = U.col(0) * V.t();


    // Soft thresholding for optimizing c (assuming lambda2 -> Inf)
    arma::mat ma = M*a;
    arma::vec ma_v = ma.col(0);
    c_new = Ust(ma_v, eta);

    // Calculate discrepancy between a & c
    dis = max(abs(c-c_new));
    c = c_new;
    i++;
  }
  return c;
}',depends = "RcppArmadillo")


