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
#' @export

library(spls)
sparsePLS_1<-function (x, y, K=3, eta=0.5, kappa = 0.5, eps = 1e-04, maxstep = 100)
{
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
    what <- spls.dv(Z, eta, kappa, eps, maxstep)
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
  if (!is.null(colnames(x))) {
    rownames(betahat) <- colnames(x)
  }
  if (q > 1 & !is.null(colnames(y))) {
    colnames(betahat) <- colnames(y)
  }
  object <- list(x = x, y = y, betahat = betahat, A = A, betamat = betamat,
                 new2As = new2As, mu = mu, meanx = meanx, normx = normx,
                 normy = normy, eta = eta, K = K, kappa = kappa, projection = pj)
  class(object) <- "bspls"
  object
}

QCQR <- function(A,b,alpha){
  svdA <- svd(A)
  c <- t(svdA$u) %*% b
  sigma <- svdA$d
  f<-function(miu){
    nume <- (sigma^2 + miu)^2
    denom <- (sigma^2) * (c^2)
    return(sum(nume/(denom+1e-4))-1)
  }
  miu <- uniroot(f,c(1e-4,1e30))$root
  inv <- diag(sigma/(sigma^2+miu))
  return(svdA$v%*%(inv%*%c))
}

spls.dv <-function( Z, eta, kappa, eps, maxstep )
{
  # initialization

  p <- nrow(Z)
  q <- ncol(Z)
  Znorm1 <- median( abs(Z) )
  Z <- Z / Znorm1

  # main iterations

  if ( q==1 )
  {
    # if univariate response, then just soft thresholding

    c <- ust( Z, eta )
  }

  if ( q>1 )
  {
    # if multivariate response

    M <- Z %*% t(Z)
    dis <- 10
    i <- 1

    # main iteration: optimize c and a iteratively

    # use svd solution if kappa==0.5

    if ( kappa==0.5 )
    {
      # initial value for a & c (outside the unit circle)

      c <- matrix( 10, p, 1 )
      c.old <- c

      while ( dis>eps & i<=maxstep )
      {
        # optimize a for fixed c

        mcsvd <- svd( M%*%c )
        a <- mcsvd$u %*% t(mcsvd$v)

        # optimize c for fixed a
        # soft thresholding ( assuming lambda2 -> Inf )

        c <- ust( M%*%a, eta )

        # calculate discrepancy between a & c

        dis <- max( abs( c - c.old ) )
        c.old <- c
        i <- i + 1
      }
    }

    # solve equation if 0<kappa<0.5

    if ( kappa>0 & kappa<0.5 )
    {
      kappa2 <- ( 1 - kappa ) / ( 1 - 2*kappa )

      # initial value for c (outside the unit circle)

      c <- matrix( 255, p, 1 )
      c.old <- c

      # define function for Lagrange part

      while ( dis>eps & i<=maxstep )
      {

        # optimize a for fixed c

        a<-QCQR(t(Z),kappa2*(t(Z)%*%c),1)

        # optimize c for fixed a
        # soft thresholding ( assuming lambda2 -> Inf )

        c <- ust( M%*%a, eta )

        # calculate discrepancy between a & c

        dis <- max( abs( c - c.old ) )
        c.old <- c
        i <- i + 1
      }
    }
  }

  return(c)
}

ust <-function( b, eta )
{
  b.ust <- matrix( 0, length(b), 1 )
  if ( eta < 1 )
  {
    valb <- abs(b) - eta * max( abs(b) )
    b.ust[ valb>=0 ] <- valb[ valb>=0 ] * (sign(b))[ valb>=0 ]
  }
  return(b.ust)
}
