##########################################
# Essential functions for performing FFT #
# and extrapolation                      #
##########################################

#' @importFrom RcppRoll roll_mean roll_sd
#' @importFrom stats integrate fft dnorm pgamma pnorm sd splinefun uniroot

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
  # This function

  # this function returns the index of the r2 that is the most reliable extrapolation point minus 1
  # so that the result can be directly added to the index that began r2

  # win.select is the starting index of each window.  It records the start of the most linear segments along the region
  # We initialize win.select so that the initial windows are half overlapping
  win.select<-seq(1,num.baseline.intervals*floor(interval.length/2),by=floor(interval.length/2))
  win.mean<-mean(r2[win.select])
  win.se<-sd(r2[win.select])/sqrt(num.baseline.intervals)

  i<-max(win.select)+1

  # whenever slide the window out into the tail.  If we encounter a segment a segment that is two standard deviations less linear than the average of the segments we've
  # seen so far, we stop.  If we encounter a segment that is more linear than the current mean, we get rid of the least linear segment indexed by win.select
  # and add the new segment to the win.select set.
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
  #browser()


  sort(win.select,decreasing = TRUE)[2]-1
}

extrapolate_tail<-function(log.cdf,xx.standardized,start,best,num.windows,right.side=TRUE){

  side<-2*as.integer(right.side)-1

  # Choose the interval length so that the windows will be half overlapping
  # And such that the intervals half over lap and extend roughly from 1e-2 to 1e-5
  # but if the current more reliable extrapolation locus, best, is not as far our as 1e-5, choose the window
  # size based on that point.
  interval.length<-ifelse(right.side,floor(2*(min(best,which.max(log.cdf>-log(1e-5)))-start)/(num.windows+3)),
                          -floor(2*(max(best,which.max(-log.cdf>log(1e-5)))-1-start)/(num.windows+3)))

  if(interval.length<10){warning(paste("Too few evaluation points in the",
                                       ifelse(right.side,"right","left"),"tail to make reliable extrapolation: result returned.
                                       Try increasing the number of CDF domain grid points n (optional argument)",sep=" "))}

  if(interval.length<4){
    warning(paste("Too few evaluation points in the", ifelse(right.side,"right","left"),"tail to extrapolate: NAs returned for extrapolation.
                  Try increasing the number of CDF domain grid points n (optional argument)",sep=" "))
    return(list("b" = NA,"best" = NA,"successful"=F))
  }

  # RC: IN FUTURE VERSION, CONSIDER USING A.STANDARDIZED AND B.STANDARDIZED AND THE APPROXIMATE MODE OF 0.2, 0.8 QUANTILES OF THE CDF
  # TO ZOOM IN ON A TAIL AND RUN ANOTHER FFT JUST ON THAT TAIL.  THIS COULD MEAN THAT FOR VERY EXTREME DISTRIBUTIONS WE DO A GENERAL FFT AND THEN ONE IN
  # EACH TAIL TO GET THE EXTRAPOLATION / SCALE RIGHT, AND THEN STITCH THOSE 3 FFTS TOGETHER...BUT THAT'S MORE WORK THAN NEEDED FOR NOW.

  #browser()
  r2 <- calc.R2(log.cdf[start:best], xx.standardized[start:best], interval.length)
  res <- fft.croppper3(side*r2$corr,interval.length,num.windows)
  #browser()

  best <- start+side*res
  b <- r2$corr[res+1]*r2$y.sd[res+1]/r2$x.sd[res+1]
  #a <- log.cdf[best]-xx.standardized[best]*b

  list("b" = b,"best" = best,"successful"=TRUE)
  }

log.rho.Q.easy.centered<-function(t,evals.s,ncps,a,mu.s,sigma.s){
  complex(imaginary = -t*(a+mu.s),real = -0.5*(sigma.s*t)^2)+sum( complex(imaginary=ncps*t*evals.s)/complex(real=1,imaginary=-2*t*evals.s)
                                        - 0.5*log(complex(real=1,imaginary=-2*t*evals.s))  )
}




calc.QFcdf <- function(evals, ncps=rep(0,length(evals)), sigma = 0, n = 2^16-1, parallel.sapply = sapply) {
  # This function estimates the CDF of the truncated distribution with the identity function
  # as the function of interest

  # Sanity Checks
  ####################

  if(!(is.integer(evals) | is.numeric(evals))){stop("evals must be integer or numeric")}
  evals <- as.numeric(evals)

  if(!(is.integer(ncps) | is.numeric(ncps))){stop("ncps must be integer or numeric")}
  ncps <- as.numeric(ncps)

  if(!(is.integer(sigma) | is.numeric(sigma))){stop("sigma must be integer or numeric")}
  sigma <- as.numeric(sigma)

  if(length(sigma)!=1){stop("sigma must have length 1")}

  if(length(evals)!=length(ncps)){stop("length of evals does not equal length of ncps")}

  if(any(ncps<0)){stop("ncps must be non-negative.")}

  if(is.infinite(sigma) | any(is.infinite(evals)) | any(is.infinite(ncps))){stop("All evals and ncps must be finite.")}



  if(all(evals==0)){stop("All evals are 0.")}
  if(any(evals==0)){
    evals <- evals[which(evals!=0)]
    ncps <- ncps[which(evals!=0)]
  }

  if(length(unique(evals))==1){warning("All of the eigenvalues are the same.
                                       Proceeding to result but consider using pchisq instead.")}


  all.pos <- all(evals>=0)

  all.neg <- all(evals<=0)

  pos.support <- all.pos & sigma==0
  neg.support <- all.neg & sigma==0

  if(pos.support){type <- "pos"}
  if(neg.support){type <- "neg"}

  # If all of the eigenvalues are negative, flip them around to make them all positive and flip them back over later.
  if(all.neg){evals <- - evals}


  # Determine Rescaling and Centering for Numerical Stability
  # We center the distribution but standardize by the component with the largest variance, not the total variance
  #########################

  mu <- sum((ncps+1)*evals)
  v <- 2*(1+2*ncps)*(evals^2)
  Q.sd <- sqrt(sigma^2 + sum(v))
  if(is.infinite(Q.sd)){stop("The variance of the target distribution is too large to be numerically represented on this machine.  Please consider rescaling the target and try again.")}
  s <- sqrt(max(c(v,sigma^2)))

  evals.s <- evals / s
  mu.s <- mu / s
  sigma.s <- sigma / s


  # Determine on how fine a grid we need to FFT (will determine the interval on which we estimate the density)
  ####################################

  fft.conc.ineq <- function(z,nu,resid.op.norm.bd,bound) {
    # This is our concentration inequality
    ifelse(z <= nu / (4*resid.op.norm.bd),
           0.5*(z^2) / nu,
           z / (4*resid.op.norm.bd) - 0.5*nu / ((4*resid.op.norm.bd)^2) ) / log(10) - bound
  }

  nu <- 8*sum(ncps*(evals.s^2))+4*sum(evals.s^2)+sigma.s^2


  b.standardized <- uniroot(fft.conc.ineq, lower = 0, upper = 1e3*(Q.sd/s),tol = .Machine$double.eps,
                            nu = nu, resid.op.norm.bd = max(evals.s), bound = 17,
                            extendInt = "upX")$root


  if(pos.support | neg.support){

    a.standardized <- -mu.s

  }else{

    if( (all.pos | all.neg) & sigma!=0){

      a.standardized  <- qnorm(1e-17,0,sigma.s) # Only the Gaussian component opposes the chi square terms

    }else{

      # These are cases where either there are simply evals of mixed signs (with or without sigma)
      a.standardized <- -uniroot(fft.conc.ineq, lower = 0, upper = 1e3*(Q.sd/s),tol = .Machine$double.eps,
                                 nu = nu, resid.op.norm.bd = abs(min(evals.s)), bound = 17,
                                 extendInt = "upX")$root
    }
  }


  # Evaluate CDF using FFT
  ######################

  rho <- parallel.sapply(1:n, function(k, a.standardized, b.standardized, n, evals.s, ncps,mu.s,sigma.s) {
    exp(log.rho.Q.easy.centered((pi*n/(b.standardized-a.standardized))*(2*((k-1)/n)-1),
                                evals.s,ncps, a.standardized,mu.s,sigma.s)) / (pi*(2*((k-1)/n)-1))
  }, a.standardized = a.standardized, b.standardized = b.standardized, n = n, evals.s = evals.s, ncps = ncps, mu.s = mu.s, sigma.s = sigma.s)

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


  # Begin search for extrapolation points
  ########################################

  # Focus in on the tail
  ctstart.r <- which.max(fft_cdf>0.99)
  ctstart.l <- which.max(fft_cdf>0.01)-1

  # First: Eliminate Parts of the Estimated CDF that deviate below 0 or above 1 (allowed to hit 0 or 1 at the
  # last point if all of the eigenvalues have the same sign and the domain of the CDF has a bound on that side)


  r0 <- which(fft_cdf[ctstart.r:n] >= 1)

  best.r <- ifelse(length(r0)==0, n, ctstart.r-2+min(r0)) # One point in from the trouble point

  if(pos.support | neg.support){ # All eigenvalues are zero or positive and sigma == 0
    fft_cdf[1] <- 0
    l0 <- which(fft_cdf[2:ctstart.l] <= 0)+1
  }else{
    l0 <- which(fft_cdf[1:ctstart.l] <= 0)
  }

  best.l <- ifelse(length(l0)==0, 1, max(l0)+1) # One point in from the trouble point


  # Second: Eliminate points as unstable if/where the CDF stops being monotonic

  # Take -log density
  log.cdf.l <- suppressWarnings(-log(fft_cdf))
  log.cdf.r <- suppressWarnings(-log1p(-fft_cdf))

  rblips <- which(diff(log.cdf.r[ctstart.r:best.r]) < -sqrt(.Machine$double.eps))
  best.r <- ifelse(length(rblips)==0,best.r,ctstart.r-2+min(rblips))
  # Again, we take the point one in from the trouble
  lblips <- which(diff(log.cdf.l[best.l:ctstart.l]) > sqrt(.Machine$double.eps))
  best.l <- ifelse(length(lblips)==0,best.l,best.l+max(lblips))
  # Again, we take the point one in from the trouble


  # Finally: refine extrapolation points by cropping off the density when it starts becoming curvy due to numerical instability

  if(pos.support | neg.support){
    limit.l <- 0
    limit.r <- Inf
    right.tail<-extrapolate_tail(log.cdf.r,xx.standardized,ctstart.r,best.r,num.windows=20,right.side=TRUE)

    if(right.tail$successful){
      best.r<-right.tail$best
      b.r<-right.tail$b/s
      a.r <- log.cdf.r[best.r]-xx[best.r]*b.r
    }else{
      b.r <- a.r <- NA
    }

    if( best.l != 1 ){
      left.tail<-extrapolate_tail(log.cdf.l,log(xx.standardized+mu.s),ctstart.l,best.l,num.windows=20,right.side=FALSE)

      if(left.tail$successful){
        best.l<-left.tail$best
        b.l<-left.tail$b
        a.l <- log.cdf.l[best.l]-log(xx[best.l])*b.l
      }else{
        b.l <- a.l <- NA
      }

    }else{
      b.l <- NULL
      a.l <- NULL
    }
  }

  # If we have both positive and negative eigenvalues or sigma !=0
  if(!(pos.support | neg.support)){

    type <-"mixed"
    limit.r <- Inf
    limit.l <- -Inf

    left.tail<-extrapolate_tail(log.cdf.l,xx.standardized,ctstart.l,best.l,num.windows=20,right.side=FALSE)
    right.tail<-extrapolate_tail(log.cdf.r,xx.standardized,ctstart.r,best.r,num.windows=20,right.side=TRUE)

    if(left.tail$successful){
      best.l<-left.tail$best
      b.l<-left.tail$b/s
      a.l <- log.cdf.l[best.l]-xx[best.l]*b.l
    }else{
      b.l <- a.l <- NA
    }

    if(right.tail$successful){
      best.r<-right.tail$best
      b.r<-right.tail$b/s
      a.r <- log.cdf.r[best.r]-xx[best.r]*b.r
    }else{
      b.r <- a.r <- NA
    }
  }


  # If a.l and b.l return NA, that means extrapolation failed but should have been done
  # If a.l and b.l return NULL, that means that no extrapolation for them was necessary because
  # x covers limit.l


  xx = xx[best.l:best.r]
  fft_cdf = fft_cdf[best.l:best.r]

  # Flip around if neg.support
  if(all.neg){
    mu <- -mu
    xx <- -rev(xx)
    fft_cdf <- 1-rev(fft_cdf)
    limit.r <- ifelse(sigma==0,0,Inf)
    limit.l <- -Inf
    a.l.2 <- a.r
    b.l.2 <- -b.r
    a.r <- a.l
    if(sigma==0){ b.r <- b.l }else{ b.r <- -b.l }# If sigma==0, then b.l was estimated in a pos.support framework.  If sigma!=0, b.l
    # was estimated in a mixed support framework, which means to flip sides, we must also flip the sign.
    a.l <- a.l.2
    b.l <- b.l.2
    # best.l <- n - best.l + 1
    # best.r <- n - best.r + 1
  }


  list("x" = xx,
       "y" = fft_cdf,
       "type"= type,
       "n" = abs(best.r-best.l)+1,
       "interval.width"=(b-a)/n,
       "limit.l" = limit.l,
       "a.l" = a.l,
       "b.l" = b.l,
       "limit.r" = limit.r,
       "a.r" = a.r,
       "b.r" = b.r,
       "mu" = mu,
       "Q.sd" = Q.sd)
  }


##############################################
# Functions for returning a CDF/PDF function #
##############################################

eval.cdf.pos <- function(q, cdf, density = FALSE, lower.tail = TRUE, log.p = FALSE){

  if(density){

    if(log.p){

      return(ifelse(q < 0, 0, ifelse(q<cdf$x[1] & (!is.null(cdf$b.l)),suppressWarnings( log(-cdf$b.l)-(cdf$a.l + (cdf$b.l+1)*log(q)) ),
                                     ifelse(q>cdf$x[cdf$n],log(cdf$b.r)-(cdf$a.r + cdf$b.r*q),
                                            log(splinefun(cdf$x,cdf$y,method="mono")(q,deriv=1))))))
    }else{

      return(ifelse(q < 0, 0, ifelse(q<cdf$x[1] & (!is.null(cdf$b.l)),suppressWarnings( -cdf$b.l*exp(-(cdf$a.l + (cdf$b.l+1)*log(q))) ),
                                     ifelse(q>cdf$x[cdf$n],cdf$b.r*exp(-(cdf$a.r + cdf$b.r*q)),
                                            splinefun(cdf$x,cdf$y,method="mono")(q,deriv=1)))))
    }

  }else{

    if(lower.tail & !log.p){
      return(ifelse(q < 0, 0, ifelse(q<cdf$x[1] & (!is.null(cdf$b.l)),suppressWarnings( exp(-(cdf$a.l + cdf$b.l*log(q))) ),
                                     ifelse(q>cdf$x[cdf$n],suppressWarnings( -expm1(-(cdf$a.r + cdf$b.r*q)) ),
                                            splinefun(cdf$x,cdf$y,method="mono")(q)))))
    }

    if(!lower.tail & !log.p){
      return(ifelse(q < 0, 1, ifelse(q<cdf$x[1] & (!is.null(cdf$b.l)),suppressWarnings( -expm1(-(cdf$a.l + cdf$b.l*log(q))) ),
                                     ifelse(q>cdf$x[cdf$n],suppressWarnings( exp(-(cdf$a.r + cdf$b.r*q)) ),
                                            1-splinefun(cdf$x,cdf$y,method="mono")(q)))))
    }

    if(lower.tail & log.p){
      return(ifelse(q < 0, -Inf, ifelse(q<cdf$x[1] & (!is.null(cdf$b.l)),suppressWarnings( -(cdf$a.l + cdf$b.l*log(q)) ),
                                        ifelse(q>cdf$x[cdf$n],suppressWarnings( log(-expm1(-(cdf$a.r + cdf$b.r*q))) ),
                                               suppressWarnings(log(splinefun(cdf$x,cdf$y,method="mono")(q)))))))
    }
    if(!lower.tail & log.p){
      return(ifelse(q < 0, 0, ifelse(q<cdf$x[1] & (!is.null(cdf$b.l)),suppressWarnings( log(-expm1(-(cdf$a.l + cdf$b.l*log(q))) )),
                                     ifelse(q>cdf$x[cdf$n],-(cdf$a.r + cdf$b.r*q),
                                            suppressWarnings(log1p(-splinefun(cdf$x,cdf$y,method="mono")(q)))))))
    }
  }
}

eval.cdf.neg <- function(q, cdf, density = FALSE, lower.tail = TRUE, log.p = FALSE){

  if(density){

    if(log.p){

      return(ifelse(q > 0, 0, ifelse(q < cdf$x[1],log(-cdf$b.l)-(cdf$a.l + cdf$b.l*q),
                                     ifelse(q > cdf$x[cdf$n] & (!is.null(cdf$b.r)),suppressWarnings( log(-cdf$b.r)-(cdf$a.r + (cdf$b.r+1)*log(-q))),
                                            log(splinefun(cdf$x,cdf$y,method="mono")(q, deriv=1))))))

    }else{

      return(ifelse(q > 0, 0, ifelse(q < cdf$x[1],-cdf$b.l*exp(-(cdf$a.l + cdf$b.l*q)),
                                     ifelse(q > cdf$x[cdf$n] & (!is.null(cdf$b.r)),suppressWarnings( -cdf$b.r*exp(-(cdf$a.r + (cdf$b.r+1)*log(-q))) ),
                                            splinefun(cdf$x,cdf$y,method="mono")(q, deriv=1)))))
    }

  }else{

    if(lower.tail & !log.p){
      return(ifelse(q > 0, 1, ifelse(q < cdf$x[1],exp(-(cdf$a.l + cdf$b.l*q)),
                                     ifelse(q > cdf$x[cdf$n] & (!is.null(cdf$b.r)),suppressWarnings( -expm1(-(cdf$a.r + cdf$b.r*log(-q))) ),
                                            splinefun(cdf$x,cdf$y,method="mono")(q)))))
    }

    if(!lower.tail & !log.p){
      return(ifelse(q > 0, 0, ifelse(q < cdf$x[1],suppressWarnings(-expm1(-(cdf$a.l + cdf$b.l*q))),
                                     ifelse(q > cdf$x[cdf$n] & (!is.null(cdf$b.r)),suppressWarnings( exp(-(cdf$a.r + cdf$b.r*log(-q))) ),
                                            1-splinefun(cdf$x,cdf$y,method="mono")(q)))))
    }

    if(lower.tail & log.p){
      return(ifelse(q > 0, 0, ifelse(q < cdf$x[1],-(cdf$a.l + cdf$b.l*q),
                                     ifelse(q > cdf$x[cdf$n] & (!is.null(cdf$b.r)),suppressWarnings( log(-expm1(-(cdf$a.r + cdf$b.r*log(-q))) )),
                                            suppressWarnings(log(splinefun(cdf$x,cdf$y,method="mono")(q)))))))
    }
    if(!lower.tail & log.p){
      return(ifelse(q > 0, -Inf, ifelse(q < cdf$x[1],suppressWarnings(log(-expm1(-(cdf$a.l + cdf$b.l*q)))),
                                        ifelse(q > cdf$x[cdf$n] & (!is.null(cdf$b.r)),suppressWarnings( -(cdf$a.r + cdf$b.r*log(-q))),
                                               suppressWarnings(log1p(-splinefun(cdf$x,cdf$y,method="mono")(q)))))))
    }
  }
}

eval.cdf.mixed <- function(q, cdf, density = FALSE, lower.tail = TRUE, log.p = FALSE){

  if(density){

    if(log.p){

      return(ifelse(q<cdf$x[1],log(-cdf$b.l)-(cdf$a.l + cdf$b.l*q),
                    ifelse(q>cdf$x[cdf$n],log(cdf$b.r)-(cdf$a.r + cdf$b.r*q),
                           log(splinefun(cdf$x,cdf$y,method="mono")(q,deriv=1)))))
    }else{

      return(ifelse(q<cdf$x[1],-cdf$b.l*exp(-(cdf$a.l + cdf$b.l*q)),
                    ifelse(q>cdf$x[cdf$n],cdf$b.r*exp(-(cdf$a.r + cdf$b.r*q)),
                           splinefun(cdf$x,cdf$y,method="mono")(q,deriv=1))))
    }

  }else{

    if(lower.tail & !log.p){
      return(ifelse(q<cdf$x[1],exp(-(cdf$a.l + cdf$b.l*q)),
                    ifelse(q>cdf$x[cdf$n],suppressWarnings(-expm1(-(cdf$a.r + cdf$b.r*q))),
                           splinefun(cdf$x,cdf$y,method="mono")(q))))
    }

    if(!lower.tail & !log.p){
      return(ifelse(q<cdf$x[1],suppressWarnings(-expm1(-(cdf$a.l + cdf$b.l*q))),
                    ifelse(q>cdf$x[cdf$n],exp(-(cdf$a.r + cdf$b.r*q)),
                           1-splinefun(cdf$x,cdf$y,method="mono")(q))))
    }

    if(lower.tail & log.p){
      return(ifelse(q<cdf$x[1],-(cdf$a.l + cdf$b.l*q),
                    ifelse(q>cdf$x[cdf$n],suppressWarnings(log(-expm1(-(cdf$a.r + cdf$b.r*q)))),
                           suppressWarnings(log(splinefun(cdf$x,cdf$y,method="mono")(q))))))
    }
    if(!lower.tail & log.p){
      return(ifelse(q<cdf$x[1],suppressWarnings(log(-expm1(-(cdf$a.l + cdf$b.l*q)))),
                    ifelse(q>cdf$x[cdf$n],-(cdf$a.r + cdf$b.r*q),
                           suppressWarnings(log1p(-splinefun(cdf$x,cdf$y,method="mono")(q))))))
    }
  }
}

wrap.QFcdf <- function(cdf){
  # Returns a function that will evaluate the CDF pointwise
  switch(cdf$type,
         mixed =   function(x, density = FALSE, lower.tail = TRUE, log.p = FALSE) eval.cdf.mixed(x, cdf, density = density, lower.tail = lower.tail, log.p = log.p),
         pos = function(x, density = FALSE, lower.tail = TRUE, log.p = FALSE) eval.cdf.pos(x, cdf, density = density, lower.tail = lower.tail, log.p = log.p),
         neg = function(x, density = FALSE, lower.tail = TRUE, log.p = FALSE) eval.cdf.neg(x, cdf, density = density, lower.tail = lower.tail, log.p = log.p)
  )
}
