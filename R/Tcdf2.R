

calc.R2<-function(yy,xx,interval.length){

  x.bar <- roll_mean(xx, interval.length)
  y.bar <- roll_mean(yy, interval.length)
  x.sd <- roll_sd(xx, interval.length)
  y.sd <- roll_sd(yy, interval.length)
  xy.sum <- interval.length*roll_mean(xx*yy, interval.length)
  corr <- (xy.sum - interval.length*x.bar*y.bar)/(interval.length*x.sd*y.sd)

  list("corr"=corr,"x.bar"=x.bar,"y.bar"=y.bar,"x.sd"=x.sd,"y.sd"=y.sd)
}

fft.croppper3<-function(r2,interval.length,num.baseline.intervals){
  # this function returns the index of the r2 that is the most reliable extrapolation point minus 1
  # so that the result can be directly added to the index that began r2

  #   win.select<-seq(1,20*interval.length,by=interval.length)
  #   win.mean<-mean(r2[win.select])
  #   win.se<-sd(r2[win.select])/sqrt(20)
  #
  #   i<-20*interval.length+1


  win.select<-seq(1,num.baseline.intervals*floor(interval.length/2),by=floor(interval.length/2))
  win.mean<-mean(r2[win.select])
  win.se<-sd(r2[win.select])/sqrt(num.baseline.intervals)

  i<-max(win.select)+1

  while(i<=(length(r2)-interval.length)){

    if(r2[i] < win.mean-2*win.se){
      break
    }

    if(r2[i] > win.mean){
      win.select[which.min(r2[win.select])]<-i
      win.mean<-mean(r2[win.select])
      win.se<-sd(r2[win.select])/sqrt(num.baseline.intervals)
      i<-i+interval.length
    }else{
      i<-i+1
    }
  }
  # For added robustness, don't take the linear interval that is furthest out, but the one that is second to the furthest out
  sort(win.select,decreasing = T)[2]-1

}

extrapolate_tail<-function(log.cdf,xx.standardized,start,best,num.windows,right.side=T){

  side<-2*as.integer(right.side)-1

  # Choose the interval length so that the windows will be half overlapping
  # And such that the intervals half over lap and extend roughly from 1e-2 to 1e-5
  interval.length<-ifelse(right.side,floor(2*(which.max(log.cdf>-log(1e-5))-start)/(num.windows+1)),
                          -floor(2*(which.max(-log.cdf>log(1e-5))-1-start)/(num.windows+1)))

  if(interval.length<10){warning(paste("Too few evaluation points in the", side,"tail to make reliable extrapolation: result returned.  Try increasing the number of CDF domain grid points n (optional argument)",sep=""))}

  if(interval.length<4){
    warning(paste("Too few evaluation points in the", side,"tail to extrapolate: NAs returned for extrapolation.  Try increasing the number of CDF domain grid points n (optional argument)",sep=""))
    return(list("b" = NA,"best" = NA,"successful"=F))
    }


  r2 <- calc.R2(log.cdf[start:best], xx.standardized[start:best], interval.length)
  res <- fft.croppper3(side*r2$corr,interval.length,num.windows)
  best <- start+side*res
  b <- r2$corr[res+1]*r2$y.sd[res+1]/r2$x.sd[res+1]
  #a <- log.cdf[best]-xx.standardized[best]*b
  #browser()

  list("b" = b,"best" = best,"successful"=T)
}


# This finds cdf for Y.s - mu.s where Y = \sum_i evals.s_i * (Z_i + \delta_i)^2
log.rho.Q.easy.centered<-function(t,evals.s,ncps,a,mu.s){
  complex(imaginary = -t*(a+mu.s))+sum( complex(imaginary=ncps*t*evals.s)/complex(real=1,imaginary=-2*t*evals.s) - 0.5*log(complex(real=1,imaginary=-2*t*evals.s))  )
}



Tcdf2 <- function(evals, ncps=rep(0,length(evals)), n = 2^16-1, qfft.apply = sapply) {
  # This function estimates the CDF of the truncated distribution with the identity function
  # as the function of interest

  # Sanity Checks
  ####################

  if(length(evals)!=length(ncps)){stop("length of evals does not equal length of ncps")}

  # Set Up
  #########################

  # Standardize by the component with the largest variance
  mu <- sum((ncps+1)*evals)
  v <- 2*(1+2*ncps)*(evals^2)
  sigma <- sqrt(sum(v))
  s <- sqrt(max(v))

  evals.s <- evals / s
  mu.s <- mu / s


  # Determine on how fine a grid we need to FFT (will determine the interval on which we estimate the density)
  ####################################

  esrange<-range(evals.s) # these eigenvalues are not centered

  b.standardized <- uniroot(function(z,nu,resid.op.norm.bd,bound) {
    # This is our concentration inequality
    ifelse(z <= nu / (4*resid.op.norm.bd),
           0.5*(z^2) / nu,
           z / (4*resid.op.norm.bd) - 0.5*nu / ((4*resid.op.norm.bd)^2) ) / log(10) - bound
  }, lower = 0, upper = 1e3*sigma/s,tol = .Machine$double.eps,
  nu = 8*sum(ncps*(evals.s^2))+4*sum(evals.s^2), resid.op.norm.bd = esrange[2], bound = 17,
  extendInt = "upX")$root

  if(esrange[2]<=0){ # All eigenvalues are zero or negative
    b.standardized <- min(-mu.s,b.standardized) # If you know the CDF must stop at 0, then don't worry about evaluating past that
  }


  a.standardized <- -uniroot(function(z,nu,resid.op.norm.bd,bound) {
    # This function is negative log of H
    ifelse(z <= nu / (4*resid.op.norm.bd),
           0.5*(z^2) / nu,
           z / (4*resid.op.norm.bd) - 0.5*nu / ((4*resid.op.norm.bd)^2) ) / log(10) - bound
  }, lower = 0, upper = 1e3*sigma/s,tol = .Machine$double.eps,
  nu = 8*sum(ncps*(evals.s^2))+4*sum(evals.s^2), resid.op.norm.bd = abs(esrange[1]), bound = 17,
  extendInt = "upX")$root

  if(esrange[1]>=0){ # All eigenvalues are zero or positive
    a.standardized <- max(-mu.s,a.standardized) # If you know the CDF must stop at 0 (on this scale, mu.s), then don't worry about evaluating past that
  }

  # Evaluate Density
  ######################

  rho <- sapply(1:n, function(k, a.standardized, b.standardized, n, evals.s, ncps,mu.s) {
    exp(log.rho.Q.easy.centered((pi*n/(b.standardized-a.standardized))*(2*((k-1)/n)-1), evals.s,ncps, a.standardized,mu.s)) / (pi*(2*((k-1)/n)-1))
  }, a.standardized = a.standardized, b.standardized = b.standardized, n = n, evals.s = evals.s, ncps = ncps, mu.s = mu.s)


  # Perform FFT
  fft_cdf <- 0.5 - (Im(fft(rho))*(-1)^(0:(n-1)))/n


  xx.standardized <- seq(a.standardized, b.standardized-(b.standardized-a.standardized)/n, length.out = n)

  a<-(mu.s+a.standardized)*s
  b<-(mu.s+b.standardized)*s
  xx<-seq(a,b-(b-a)/n,length.out = n)


  # Stop if FFT returns something nonsensical
  if(any(is.nan(fft_cdf))){
    stop("fft returned NaN values")
  }
  if(any(is.na(fft_cdf))){
    stop("fft returned NA values")
  }
  if(any(is.infinite(fft_cdf))){
    stop("fft returned infinite values")
  }


  # Find linear extrapolation points
  ####################################

  ctstart.r <- which.max(fft_cdf>0.99)
  ctstart.l <- which.max(fft_cdf>0.01)-1

  # First: Theoretically, the CDF should never hit 0 or 1...well depends on whether all the eigenvalues have the
  # same sign.  In that case, just force the CDF to hit zero there and kill all points leading up to 0 that hit or are below it.
  # best.l and best.r will be defined to be one point in from where the numerical instability starts

  if(esrange[2]<=0){ # All eigenvalues are zero or negative
    fft_cdf[n] <- 1
    r0 <- which( fft_cdf[ctstart.r:(n-1)]>=1 )
  }else{
    r0 <- which( fft_cdf[ctstart.r:n]>=1 )
  }

  if(esrange[1]>=0){ # All eigenvalues are zero or positive
    fft_cdf[1] <- 0
    l0 <- which( fft_cdf[2:ctstart.l]<=0 )+1
  }else{
    l0 <- which( fft_cdf[1:ctstart.l]<=0 )
  }

  best.r <- ifelse(length(r0)==0, n, ctstart.r-2+min(r0))
  best.l <- ifelse(length(l0)==0, 1, max(l0)+1)

  # Take -log density
  log.cdf.l <- suppressWarnings(-log(fft_cdf))
  log.cdf.r <- suppressWarnings(-log1p(-fft_cdf))

  # Second: CDF should be monotonic
  rblips <- which(diff(log.cdf.r[ctstart.r:best.r]) < -sqrt(.Machine$double.eps))
  best.r <- ifelse(length(rblips)==0,best.r,ctstart.r-2+min(rblips))
  lblips <- which(diff(log.cdf.l[best.l:ctstart.l]) > sqrt(.Machine$double.eps))
  best.l <- ifelse(length(lblips)==0,best.l,best.l+max(lblips))


  # Refine extrapolation points by cropping off the density when it starts becoming curvy due to numerical instability

  if(esrange[1]>=0){ # All eigenvalues are zero or positive
    limit.l <- 0
    limit.r <- Inf
    right.tail<-extrapolate_tail(log.cdf.r,xx.standardized,ctstart.r,best.r,num.windows=20,right.side=T)

    if(right.tail$successful){
      best.r<-right.tail$best
      b.r<-right.tail$b/s
      a.r <- log.cdf.r[best.r]-xx[best.r]*b.r
      c.r <- 0
    }else{
      c.r <- b.r <- a.r <- NA
    }


    if( best.l != 1 ){
      left.tail<-extrapolate_tail(log.cdf.l,xx.standardized,ctstart.l,best.l,num.windows=20,right.side=F)

      if(left.tail$successful){
        best.l<-left.tail$best
        b.l<-left.tail$b/s
        a.l <- log.cdf.l[best.l]-xx[best.l]*b.l
        c.l <- -exp(-(a.l+b.l*(-mu.s)))
      }else{
        c.l <- b.l <- a.l <- NA
      }

    }else{
      b.l <- NULL
      a.l <- NULL
      c.l <- NULL
    }
  }

  if(esrange[2]<=0){ # All eigenvalues are zero or negative
    limit.l <- Inf
    limit.r <- 0
    left.tail<-extrapolate_tail(log.cdf.l,xx.standardized,ctstart.l,best.l,num.windows=20,right.side=F)

    if(left.tail$successful){
      best.l<-left.tail$best
      b.l<-left.tail$b/s
      a.l <- log.cdf.l[best.l]-xx[best.l]*b.l
      c.l <- 0
    }else{
     c.l <- b.l <- a.l <- NA
    }

    if( best.r != n ){
      right.tail<-extrapolate_tail(log.cdf.r,xx.standardized,ctstart.r,best.r,num.windows=20,right.side=T)

      if(right.tail$successful){
        best.r<-right.tail$best
        b.r<-right.tail$b/s
        a.r <- log.cdf.r[best.r]-xx[best.r]*b.r
        c.r <- -exp(-(a.r+b.r*(-mu.s)))

      }else{
        c.r <- b.r <- a.r <- NA
      }

    }else{
      b.r <- NULL
      a.r <- NULL
      c.r <- NULL
    }
  }

  if(esrange[1]<0 & esrange[2]>0){ # We have both positive and negative eigenvalues

    limit.r <- Inf
    limit.l <- -Inf

    left.tail<-extrapolate_tail(log.cdf.l,xx.standardized,ctstart.l,best.l,num.windows=20,right.side=F)
    right.tail<-extrapolate_tail(log.cdf.r,xx.standardized,ctstart.r,best.r,num.windows=20,right.side=T)

    if(left.tail$successful){
      best.l<-left.tail$best
      b.l<-left.tail$b/s
      a.l <- log.cdf.l[best.l]-xx[best.l]*b.l
      c.l <- 0
    }else{
      c.l <- b.l <- a.l <- NA
    }

    if(right.tail$successful){
      best.r<-right.tail$best
      b.r<-right.tail$b/s
      a.r <- log.cdf.r[best.r]-xx[best.r]*b.r
      c.r <- 0
    }else{
      c.r <- b.r <- a.r <- NA
    }

  }



  list("x" = xx,
       "y" = fft_cdf,
       "interval.width"=(b-a)/n,
       "limit.l" = limit.l,
       "index.l" = best.l,
       "a.l" = a.l,
       "b.l" = b.l,
       "c.l" = c.l,
       "limit.r" = limit.r,
       "index.r" = best.r,
       "a.r" = a.r,
       "b.r" = b.r,
       "c.r" = c.r)


  # The function returns the full output of the FFT
  # if limit.l
}

system.time(test<-Tcdf2(evals=l.e.vals,n=2^14-1))

plot(test$x,test$y,type="l",col="blue")

#plot(log(test$x[1:1000]),-log(test$y)[1:1000]/log(10),type="l",col="blue",ylim=c(0,50))
abline(test$a.l,test$b.l)

xx<-seq(min(test$x[1:1000],max(test$x[1:1000])),len=10000)
lines(xx,-log(exp(-(test$a.l+xx*test$b.l))+test$c.l))
abline(v=test$x[test$index.l])

plot(test$x,-log1p(-test$y),type="l",col="blue")
abline(test$a.r,test$b.r)
abline(v=test$x[test$index.r])

plot(test$x,-log(test$y),type="l",col="blue")
abline(test$a.l,test$b.l)
abline(v=test$x[test$index.l])




xxx<-seq(0,1e-1,len=1000)
plot(log(xxx),log(-pchisq(xxx,df=2,lower.tail = T,log.p = T,ncp=100)),type="l")

# Gotta come back and make sure that former best.l are not discarded and also implement
# the interpolation back to 0.











#
#
# log.rho.gamma<-function(t,alpha,beta,a){
#   complex(imaginary = -t*a)-alpha*log(complex(real=1,imaginary = -t/beta))
# }
#
# # THIS WORKS!
# log.rho.Q.easy<-function(t,evals,ncps,a){
#   complex(imaginary = -t*a)+sum( complex(imaginary=ncps*t*evals)/complex(real=1,imaginary=-2*t*evals) - 0.5*log(complex(real=1,imaginary=-2*t*evals))  )
# }
#
# a.standardized <- -1e6
# b.standardized <- 1e6
# n<-2^16-1
#
# evals<-l.e.vals
# ncps<-rep(0,length(evals))
#
# rho <- sapply(1:n, function(k, a.standardized, b.standardized, n, evals, ncps) {
#   exp(log.rho.Q.easy((pi*n/(b.standardized-a.standardized))*(2*((k-1)/n)-1), evals,ncps, a.standardized)) / (pi*(2*((k-1)/n)-1))
# }, a.standardized = a.standardized, b.standardized = b.standardized, n = n, evals = evals, ncps = ncps)
#
# # Perform FFT
# fft_cdf <- 0.5 - (Im(fft(rho))*(-1)^(0:(n-1)))/n
# xx.standardized <- seq(a.standardized, b.standardized-(b.standardized-a.standardized)/n, length.out = n)
# #plot(xx.standardized,pgamma(xx.standardized,shape = alpha,rate=beta),type="l")
# lines(xx.standardized,fft_cdf,col="blue")
#
# ##########################
#
# # Trial below with centering and scaling
#
# s<-sqrt(max(2*(1+2*ncps)*(evals^2)))
# evals.s <- evals/s
# mu.s <- sum((ncps+1)*evals)/s
#
# # This finds cdf for Y.s - mu.s (both have been scaled by s)
# log.rho.Q.easy.centered<-function(t,evals.s,ncps,a,mu.s){
#   complex(imaginary = -t*(a+mu.s))+sum( complex(imaginary=ncps*t*evals.s)/complex(real=1,imaginary=-2*t*evals.s) - 0.5*log(complex(real=1,imaginary=-2*t*evals.s))  )
# }
#
# a.standardized <- -1e3
# b.standardized <- 1e3
# n<-2^16-1
#
# rho <- sapply(1:n, function(k, a.standardized, b.standardized, n, evals.s, ncps,mu.s) {
#   exp(log.rho.Q.easy.centered((pi*n/(b.standardized-a.standardized))*(2*((k-1)/n)-1), evals.s,ncps, a.standardized,mu.s)) / (pi*(2*((k-1)/n)-1))
# }, a.standardized = a.standardized, b.standardized = b.standardized, n = n, evals.s = evals.s, ncps = ncps, mu.s = mu.s)
#
# # Perform FFT
# fft_cdf.centered <- 0.5 - (Im(fft(rho))*(-1)^(0:(n-1)))/n
# xx.centered <- seq(a.standardized, b.standardized-(b.standardized-a.standardized)/n, length.out = n)*s + mu.s*s
# #plot(xx.standardized,pgamma(xx.standardized,shape = alpha,rate=beta),type="l")
# lines(xx.centered,fft_cdf.centered,col="blue")
#
#
# #
#
#
#
#
#
#
#
#
#
# log.rho.func<-function(t,phi,ncps,a){
#   lp1<-log1p((t*phi)^2)
#   Omega<-exp(log(ncps)+2*log(abs(t*phi))-lp1)
#   complex(real=-0.5*sum(Omega+lp1),imaginary = -t*a-0.5*sum((Omega+1)*t*phi-atan(t*phi)))
# }
#
# rho.old <- sapply(1:n, function(k, a.standardized, b.standardized, n, evals, ncps) {
#   exp(log.rho.func((pi*n/(b.standardized-a.standardized))*(2*((k-1)/n)-1), 2*evals,ncps, a.standardized)) / (pi*(2*((k-1)/n)-1))
# }, a.standardized = a.standardized, b.standardized = b.standardized, n = n, evals = evals, ncps = ncps)
#
# # Perform FFT
# fft_cdf.old <- 0.5 - (Im(fft(rho.old))*(-1)^(0:(n-1)))/n
# #plot(xx.standardized,pgamma(xx.standardized,shape = alpha,rate=beta),type="l")
# lines(xx.standardized+E,fft_cdf.old,col="red")
#
#
# # Wowza these agree!!!
#



