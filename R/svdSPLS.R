#' A basic implementation of constrained PLSSVD
#' @name svdSPLS
#' @param X description a n*q matrix
#' @param Y description a n*p matrix
#' @return a list of pls components
#' @examples
#' X=matrix(rnorm(200),10,20)
#' Y=matrix(rnorm(100),5,20)
#' fit = sparsePLS(X,Y,2)
#' @importFrom pls plsr
#' @importFrom irlba irlba
#' @export

sparsePLS_0<-function(X,Y,h=2,fullrank=TRUE,maxiter=100,lambda1=0.5,lambda2=0.5){
  #get dimension
  p <- ncol(X)
  q <- ncol(Y)
  #centralize
  X <- scale(X, center=TRUE, scale=TRUE)
  Y <- scale(Y, center=TRUE, scale=TRUE)
  #initialization
  loadingX <- NULL
  loadingY <- NULL
  Xscore<-NULL
  Yscore<-NULL
  #ScoreX<-matrix(0, nrow=p, ncol=h)
  #ScoreY<-matrix(0, nrow=q, ncol=h)

  constrain<-function(lambda,V){
    dis = abs(V)-lambda
    dis[dis<0]=0
    return(sign(V)*dis)
  }

  #iteration
  for(i in c(1:h)){
    M<-t(X) %*% Y
    if(fullrank){
      M_svd<-svd(M)
    }else{
      M_svd<-irlba(M, 1)
    }
    u_old = M_svd$u[,1]
    v_old = M_svd$v[,1]

    tol = 1e-4
    iter = 1
    u_new = rep(0,p)
    v_new = rep(0,q)
    #update direction vector until convergence
    while(iter<maxiter){
      #add lasso and normalize
      u_temp = constrain(lambda1,M%*%v_old)
      u_new = u_temp/sqrt(sum(u_temp^2))
      v_temp = constrain(lambda2,t(M)%*%u_old)
      v_new = v_temp/sqrt(sum(v_temp^2))
      if(max(abs(u_old-u_new)) < tol){break}
      u_old = u_new
      v_old = v_new
      iter = iter + 1
    }

    #calculate PLS components of X and Y
    norm_u = sum(u_new^2)
    norm_v = sum(v_new^2)
    e = (X %*% u_new)/norm_u
    w = (Y %*% v_new)/norm_v

    #save weight vector

    Xscore=cbind(Xscore,e)
    Yscore=cbind(Yscore,w)


    #calculate loading of X and Y
    norm_e = sum(e^2)
    norm_w = sum(w^2)
    inner_ew = sum(e*w)
    c = (t(X) %*% e)/norm_e
    d = (t(Y) %*% e)/inner_ew
    #f = (t(Y) %*% w)/norm_w
    loadingX = cbind(c,loadingX)
    loadingY = cbind(d,loadingY)

    #deflation
    X = X - e%*%t(c)
    Y = Y - e%*%t(d)
  }
  return(list(loadingX=loadingX,loadingY=loadingY,Xscore=Xscore,Yscore=Yscore))
}
