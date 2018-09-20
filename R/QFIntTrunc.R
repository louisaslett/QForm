QFIntTrunc<-function(obs, M, mu, sigma, k, neg.log.p = TRUE){
  
  
  # All input must be on the original scale.  obs= y^T M y
  
  if(class(M) != "dsyMatrix") {
    stop("Matrix M must be symmetric of class dsyMatrix (from Matrix package).")
  }
  if(!all(diag(M) == 0)) {
    stop("M must be a hollow matrix (zero on diagonal).")
  }
  if(!is.vector(obs)) {
    stop("obs must be a vector.")
  }
  # obs <- suppressMessages(as.numeric(crossprod(crossprod(M, y), y)))
  N <- nrow(M)
  if(!is.vector(sigma) || length(sigma) != N) {
    stop("sigma must be a vector matching the number of rows of M.")
  }
  if(!is.vector(mu) || length(mu) != N) {
    stop("mu must be a vector matching the number of rows of M.")
  }
  
  # Eigen-decompose and rescale M by N in order to keep magnitude of observations in a more maneagable range
  q <- obs/N
  mu.tilde <- mu/sigma
  sigma <- sigma/sqrt(N) # Note: observation vector is rescaled just above here
  M.tilde <- sweep(M * sigma, 2, sigma, "*")
  e <- eigs(M.tilde, k, which = "LM")
  evec.tilde <- e$vectors[, order(abs(e$values), decreasing=TRUE)]
  eval.tilde <- e$values[order(abs(e$values), decreasing=TRUE)]
  ncps <- c(crossprod(evec.tilde, mu.tilde))^2

  
  # Run FFT  
  t.cdf <- Tcdf_new(eval.tilde,ncps)
  
  if(q<=t.cdf$x[t.cdf$index.l]){
    res <- t.cdf$a.l+t.cdf$b.l*q
  }
  
  if(t.cdf$x[t.cdf$index.l] < q){
    res <- -log(c(t.cdf$y[t.cdf$index.l:t.cdf$index.r],1)[which.min(q > c(t.cdf$x[t.cdf$index.l:t.cdf$index.r],Inf))])
  }
  
  if(neg.log.p==T){  
    return(res)
  }
  else{
    return(exp(-res))
  }
  
}

# 
# QFIntTrunc(2*18408000,abo_min_mat,mu=rep(0.5,nrow(abo_min_mat)),sigma = rep(0.5,nrow(abo_min_mat)),10)
# 
# ## Why is Exp returning that the area under the curve is so massive to the left?????  
# # Bug in ExpInt?  Seems to be throwing off everything
# 
# xx <- seq(18408000,18408000*2,len=10)
# yy <- xx
# require(Matrix)
# for(i in 1:length(xx)){
#   yy[i] <- QFIntTrunc(xx[i],abo_min_mat,mu=rep(0.5,nrow(abo_min_mat)),sigma = rep(0.5,nrow(abo_min_mat)),10)
# }
# plot(xx,yy)
  