ETruncNorm <- function(upper,mu,sigma){
  # Return log(E(Z)) where Z \sim TruncNormal(mu,sigma,0,upper)

  if(upper==0){return(-Inf)}

  a <- -mu/sigma
  b <- (upper-mu)/sigma

  if(mu==0 & abs(a)==abs(b)){return(-Inf)}

  Pd<-Phidiff(a,b)

  lmu <- log(abs(mu))

  if(mu!=0 & abs(a)==abs(b)){
    # Then mu must be positive since we know overall integral must be positive.  (and the only way abs(a) and abs(b) could be equal is if mu = upper/2 >0)
    return( lmu + Pd )
  }

  apd <- absphidiff(a,b)


  if(mu==0 & abs(a)!=abs(b)){
    # Then a=0 so phi(a) > phi(b)
    return( log(sigma) + apd - Pd )
  }

  # If none of the above cases, we must have mu!=0 & abs(a)!=abs(b).

  core <- expm1(apd - 2*Pd + log(sigma)-lmu)

  if(mu < 0){
    # If mu < 0 , then a is some positive number smaller than b, so phi(a) > phi(b)
    return( log(core) )
  }

  if(abs(a) < abs(b) & mu > 0){
    return( log(2) + log1p(core/2) )
  }

  if(abs(a) > abs(b) & mu > 0){
    return( log(-core) )
  }

  # Since we're integrating from 0 to some upper, if mu is negative, then we must have that phi(a)

}


Phidiff <- function(a, b){
  # We expect a <= b
  # Calc log( \Phi(b) - \Phi(a) )
  pnorm(a,log.p = T) + log(expm1( pnorm(b,log.p = T) - pnorm(a,log.p = T)  ))
}

absphidiff <- function(a,b){
  # Calc log( |phi(a)-phi(b)| )
  a. <- dnorm(a,log=T)
  b. <- dnorm(b,log=T)

  a. + log(abs(expm1(b.-a.)))
}

ETruncNormVec<-Vectorize(FUN=ETruncNorm,vectorize.args = "upper",SIMPLIFY = T)
xx<-seq(0,100,len=1000)
plot(xx,exp(ETruncNormVec(xx,mu=0,sigma=1)))
