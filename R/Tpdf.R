Tpdf <- function(evals, ncps, n = 2^17, neg.log = FALSE, num.baseline.intervals = 20, min.interval.length = 20, qfft.apply = sapply) {
  # This function estimates the CDF of the truncated distribution with the identity function
  # as the function of interest


  # Sanity Checks
  ####################

  if(length(evals)!=length(ncps)){stop("length of evals does not equal length of ncps")}


  # Set Up
  #########################

  # Standardize by the component with the largest variance
  E <- sum((ncps+1)*evals)
  SD.max <- sqrt(max(2*(1+2*ncps)*(evals^2)))

  E.tilde <- E/SD.max
  evals.tilde <- 2*evals/SD.max

  # Determine on how fine a grid we need to FFT (will determine the interval on which we estimate the density)
  ####################################

  b.standardized <- uniroot(function(z,nu,resid.op.norm.bd,bound) {
    # This function is negative log of H
    ifelse(z <= nu / (4*resid.op.norm.bd),
           0.5*(z^2) / nu,
           z / (4*resid.op.norm.bd) - 0.5*nu / ((4*resid.op.norm.bd)^2) ) / log(10) - bound
  }, lower = 0, upper = 1e6,tol = .Machine$double.eps,
  nu = 2*sum(ncps*(evals.tilde^2))+sum(evals.tilde^2), resid.op.norm.bd = max(abs(evals.tilde))/2, bound = 20,
  extendInt = "upX")$root

  a.standardized <- - b.standardized

  # Evaluate Density
  ######################3
  rho <- qfft.apply(1:n, function(k, a.standardized, b.standardized, n, evals.tilde, ncps) {
    exp(log.rho.func((pi*n/(b.standardized-a.standardized))*(2*((k-1)/n)-1), evals.tilde, ncps, a.standardized))
  }, a.standardized = a.standardized, b.standardized = b.standardized, n = n, evals.tilde = evals.tilde, ncps = ncps)

  # Perform FFT
  ans <- Re(fft(rho))*(-1)^(0:(n-1))/(b.standardized-a.standardized)
  xx.standardized <- seq(a.standardized, b.standardized-(b.standardized-a.standardized)/n, length.out = n)

  # Stop if FFT returns something nonsensical
  if(any(is.nan(ans))){
    stop("fft returned NaN values")
  }
  if(any(is.na(ans))){
    stop("fft returned NA values")
  }
  if(any(is.infinite(ans))){
    stop("fft returned infinite values")
  }

  # Note: this is not a numerically safe cumsum because we have not sorted the values.
  # This is just to help up obtain the 0.1 and 0.01 quantiles
  # and verify if something went grossly wrong with fft

  fft_cdf <- cumsum(ans)*(b.standardized-a.standardized)/n

  if( (fft_cdf[n]-1)>1e-13 || (1-fft_cdf[n])>1e-13 ) {
    warning("CDF implied by fft did not normalize to 1.  Results may be numerically unstable.")
  }

  # Find linear extrapolation points
  ####################################

  ctstart.r <- min(which(fft_cdf>0.9))
  ctstart.l <- max(which(fft_cdf<0.1))

  tstart.r <- min(which(fft_cdf>0.99))
  tstart.l <- max(which(fft_cdf<0.01))

  # First: Screen for 0s in the density
  # best.l and best.r will be defined to be one point in from where the numerical instability starts
  r0 <- which( ans[tstart.r:n]<=0 )
  best.r <- ifelse(length(r0)==0, n, tstart.r-2+min(r0))
  l0 <- which( ans[1:tstart.l]<=0 )
  best.l <- ifelse(length(l0)==0, 1, max(l0)+1)

  # Take -log density
  ans2 <- suppressWarnings(-log(ans))

  # Second: screen for increases in density
  rblips <- which(diff(ans2[tstart.r:best.r]) < -sqrt(.Machine$double.eps))
  best.r <- ifelse(length(rblips)==0,best.r,tstart.r-2+min(rblips))
  lblips <- which(diff(ans2[best.l:tstart.l]) > sqrt(.Machine$double.eps))
  best.l <- ifelse(length(lblips)==0,best.l,best.l+max(lblips))

  interval.length.r<-floor((tstart.r-ctstart.r)/num.baseline.intervals)
  interval.length.l<-floor((ctstart.l-tstart.l)/num.baseline.intervals)

  # This is testing that the number of points between cstart.r and best.r are sufficient
  if(interval.length.r < min.interval.length){stop("Need to increase n")}
  if(floor((best.r-ctstart.r) / interval.length.r)<40){stop("Need to increase n")}

  if(interval.length.l < min.interval.length){stop("Need to increase n")}
  if(floor((ctstart.l-best.l) / interval.length.l)<40){stop("Need to increase n")}



  # Refine extrapolation points by cropping off the density when it starts becoming curvy due to numerical instability

  # Input to calc.R2 should be truncated and be roughly ordered to increase from left to right
  # (until region of numerical instability)

  r2.r <- calc.R2(ans2[ctstart.r:best.r], xx.standardized[ctstart.r:best.r], interval.length.r)
  res <- fft.croppper2(r2.r$corr,interval.length.r,num.baseline.intervals)
  best.r <- ctstart.r+res
  b.r <- r2.r$corr[res+1]*r2.r$y.sd[res+1]/r2.r$x.sd[res+1]
  a.r <- r2.r$y.bar[res+1]-b.r*r2.r$x.bar[res+1]
  limit.r <- Inf

  r2.l <- calc.R2(ans2[ctstart.l:best.l], xx.standardized[ctstart.l:best.l], interval.length.l)
  res <- fft.croppper2(-r2.l$corr, interval.length.l, num.baseline.intervals)
  best.l <- ctstart.l-res
  b.l <- r2.l$corr[res+1]*r2.l$y.sd[res+1]/r2.l$x.sd[res+1]
  a.l <- r2.l$y.bar[res+1]-b.l*r2.l$x.bar[res+1]
  limit.l <- -Inf


  # Scale and shift solutions back to the original real line

  a<-E+SD.max*a.standardized
  b<-E+SD.max*b.standardized
  xx<-seq(a,b-(b-a)/n,length.out = n)

  b.r <- b.r/SD.max
  b.l <- b.l/SD.max

  a.r <- ans2[best.r]-xx[best.r]*b.r+log(SD.max)
  a.l <- ans2[best.l]-xx[best.l]*b.l+log(SD.max)


  # Change results if the support of the distribution is the positive or negative real line


  if(all(evals.tilde>=0)){
    # Support is only over the positive real line
    limit.l <- 0
    if(xx[best.l]<=0){
      # Set best.l to zero or if zero is not in xx, to the negative number that is closest to zero
      # so that O lies in the interval between best.l and the next point to the right
      best.l <- which.max(xx>0)-1
      a.l <- -Inf
      b.l <- 0
    }
  }


  if(all(evals.tilde<=0)){
    # Support is only over the negative real line
    limit.r <- 0
    if(xx[best.r]>=0){
      # Set best.r to zero or if zero is not in xx, to the positive number that is closest to zero
      # so that O lies in the interval between best.r and the next point to the left
      best.r <- which.max(xx>=0)
      a.r <- -Inf
      b.r <- 0
    }
  }


  # LA: change this to return an object which when called will evaluate the pdf
  # pointwise with trapezoidal interpolation.  Also override the base
  # integrate function so that when it sees this object type we intercept and
  # do the exact thing, but pass through all other standard function calls to
  # base R.

  # RC: might want to change best.l and best.r to be the outer rather than inner points of the extrapolation interval
  # in order to guarantee that the extrapolated line is in fact below the -log density even over the span of the extrapolation interval
  #


  list("x" = xx,
       "y" = if(neg.log==T) { ans2+log(SD.max) } else { ans/SD.max },
       "interval.width"=(b-a)/n,
       "limit.l" = limit.l,
       "index.l" = best.l,
       "index.r" = best.r,
       "limit.r" = limit.r,
       "a.l" = a.l,
       "b.l" = b.l,
       "a.r" = a.r,
       "b.r" = b.r)
}
