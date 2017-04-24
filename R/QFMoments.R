QFMoments<-function(M,mu,sigma){

  if(class(M) != "dsyMatrix") {
    stop("Matrix M must be symmetric of class dsyMatrix (from Matrix package).")
  }
  if(!all(diag(M) == 0)) {
    stop("M must be a hollow matrix (zero on diagonal).")
  }
  sigma_2<-sigma^2

  Mmu<-mat%*%mu
  E_yBy<-sum(mu*Mmu)
  Var_yBy<-sum(4*(sigma_2)*Mmu^2+2*sigma_2*((M^2)%*%sigma_2))

  return(list("mean"=E_yBy,"var"=Var_yBy))
}
