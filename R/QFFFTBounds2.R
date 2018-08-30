#' Quadratic form bounds
#'
#' Compute bounds on the CDF of a quadratic form.
#'
#' The above was the one line summary which appears at the top of the documentation.  The gory detailed description of what is does goes here.
#'
#' Those gory details can span as many paragraphs as you want.
#'
#' @param obs scalar; the observed value of the quadratic form: that is, \eqn{y^T M y}.
#' @param evals vector of eigen-values of the matrix \eqn{M}.  These need not be all eigen values.
#' @param ncps description here
#' @param E_R description here
#' @param nu description here
#' @param N description here
#' @param resid.op.norm.bd description here
#' @param if.insuff.eigs string indicating what action to take if there are insufficient eigen-values to produce an accurate bound.  If \code{"trivial"} then the bounds \eqn{[0,1]} are returned; if \code{"missing"} then \code{NA} is returned for both bounds.
#' @param lower.tail logical; if \code{TRUE} (default), probability is \eqn{P(y^T M y \le obs)}, otherwise \eqn{P(y^T M y > obs)}
#' @param log logical; if \code{TRUE}, probability \eqn{p} is given as \eqn{log(p)}
#'
#' @return A data frame containing the variables lower and upper which provide the bounds on the CDF of the quadratic form.
#'
#' @examples
#' # Some code here which runs a self-contained example
#'
QFFFTBounds2 <- function(obs, evals, ncps, E_R, nu, N, resid.op.norm.bd, if.insuff.eigs = "trivial", lower.tail = TRUE, log = FALSE, parallel = FALSE, extended.acc = FALSE) {

  if(class(parallel)[1] != "logical") {
    if("cluster" %in% class(parallel)) {
      cl <- parallel
      stop.cl <- FALSE
    } else {
      stop.cl <- cl <- makePSOCKcluster(parallel)
    }
  } else {
    if(parallel) {
      stop.cl <- cl <- makePSOCKcluster(detectCores(TRUE))
    } else {
      stop.cl <- cl <- FALSE
    }
  }
  if("cluster" %in% class(cl)) {
    qfft.apply <- function(...) {
      parSapply(cl, ...)
    }
  } else {
    qfft.apply <- sapply
  }

  # PRELIMINARY Crm QFFFALCULATIONS

  # Rescale by N to keep calculations within a range of reasonable precision
  obs <- obs/N
  E_R <- E_R/N
  nu <- nu/(N^2)
  evals <- evals/N
  resid.op.norm.bd <- resid.op.norm.bd/N

  # Calculate the mean and variance of the truncated part of the quadratic form
  # E_T <- sum(evals*(1+ncps))
  Var_T <- 2*sum((evals^2)*(1+2*ncps))
  Var_R <- nu/2

  # COMPUTE FFT
  if(extended.acc) {
    stop("Extended accuracy not implemented yet")
  } else {
    n <- 2^20  # Number of desired intervals (computational budget)
    desired.Mod.rho.bound <- 4  #how thin we want to tail of rho to be at the end of our integration interval in units of -log10.

    # Note that the cgf is evaluated from -pi*(n/(b-a)) to pi*(n/(b-a)) in steps of 2*pi/(b-a).  We assume that a = -b.  So
    # We are evaluating from -(pi/2)*(n/b) to (pi/2)*(n/b) in steps of pi/b .
    b <- uniroot(function(b, bound) {
      Re(log.phi(-pi*n/(2*b), lambda, ncps))/log(10)+bound
    }, lower = sqrt(.Machine$double.xmin), upper = 1e6, extendInt = "upX", tol = .Machine$double.eps,
    bound = desired.Mod.rho.bound)$root
    a <- -b


    rho <- qfft.apply(1:n, function(k, a, b, n, evals, ncps) {
      exp(log.phi( (pi*n/(b-a))*(2*((k-1)/n)-1), evals, ncps) - complex(imaginary = pi*a*(n/(b-a))*(2*((k-1)/n)-1)) )
    }, a = a, b = b, n = n, evals = evals, ncps = ncps)

    dFy <- Re(fft(rho)*(-1)^(0:(n-1))/(b-a))
  }

  # FIND LINEAR EXTENSIONS
  x <- seq(a, b-(b-a)/n, length.out = n)
  ldFy <- suppressWarnings(-log10(dFy))
  lmsz <- 100 # size of window to fit linear model

  # Left tail
  r <- which.min(ldFy)
  if(any(!is.finite(ldFy[1:r]))) {
    l <- max(which(!is.finite(ldFy[1:r])))+1
  } else {
    l <- 1
  }
  x.bar <- roll_mean(x[l:r], lmsz)
  y.bar <- roll_mean(ldFy[l:r], lmsz)
  x.sd <- roll_sd(x[l:r], lmsz)
  y.sd <- roll_sd(ldFy[l:r], lmsz)
  xy.sum <- lmsz*roll_mean(x[l:r]*ldFy[l:r], lmsz)
  corr <- (xy.sum - lmsz*x.bar*y.bar)/(lmsz*x.sd*y.sd)
  best.l <- which.min(corr)
  b.l <- corr[best.l]*y.sd[best.l]/x.sd[best.l]
  a.l <- y.bar[best.l]-b.l*x.bar[best.l]
  # Choose left end of left interval ... if there is slight upward curvature,
  # the left most point will be above the line for sure, whereas choosing right
  # end would result in some of the observations in the middle of the interval
  # falling below
  cut.l <- x[best.l+l-1]

  # Right tail
  l <- which.min(ldFy)
  if(any(!is.finite(ldFy[l:length(ldFy)]))) {
    r <- min(which(!is.finite(ldFy[l:length(ldFy)])))+(l-1)-1
  } else {
    r <- length(ldFy)
  }
  x.bar <- roll_mean(x[l:r], lmsz)
  y.bar <- roll_mean(ldFy[l:r], lmsz)
  x.sd <- roll_sd(x[l:r], lmsz)
  y.sd <- roll_sd(ldFy[l:r], lmsz)
  xy.sum <- lmsz*roll_mean(x[l:r]*ldFy[l:r], lmsz)
  corr <- (xy.sum - lmsz*x.bar*y.bar)/(lmsz*x.sd*y.sd)
  best.r <- which.max(corr)
  b.r <- corr[best.r]*y.sd[best.r]/x.sd[best.r]
  a.r <- y.bar[best.r]-b.r*x.bar[best.r]
  # Choose right end of right interval ... if there is slight upward curvature,
  # the right most point will be above the line for sure, whereas choosing left
  # end would result in some of the observations in the middle of the interval
  # falling below
  cut.r <- x[best.r+l+lmsz-2]

  if(TRUE) { # Diagnostic
    tmp <- (which.min(abs(x-cut.l))-20):(which.min(abs(x-cut.l))+20)
    plot(x[tmp], ldFy[tmp], type = "l")
    lines(x[tmp], Vectorize(function(z) {
      if(z <= cut.l) {
        return(a.l + b.l*z)
      }
      if(z >= cut.r) {
        return(a.r + b.r*z)
      }
      return(ldFy[which.min(abs(x-z))])
    })(x[tmp]), col = "green")
    abline(v=cut.l, lty = 2)
    tmp <- (which.min(abs(x-cut.r))-20):(which.min(abs(x-cut.r))+20)
    plot(x[tmp], ldFy[tmp], type = "l")
    lines(x[tmp], Vectorize(function(z) {
      if(z <= cut.l) {
        return(a.l + b.l*z)
      }
      if(z >= cut.r) {
        return(a.r + b.r*z)
      }
      return(ldFy[which.min(abs(x-z))])
    })(x[tmp]), col = "green")
    abline(v=cut.r, lty = 2)
    plot(x, ldFy, type = "l")
    lines(x, Vectorize(function(z) {
      if(z <= cut.l) {
        return(a.l + b.l*z)
      }
      if(z >= cut.r) {
        return(a.r + b.r*z)
      }
      return(ldFy[which.min(abs(x-z))])
    })(x), col = "green")
    abline(v=cut.l, lty = 2)
    abline(v=cut.r, lty = 2)
  }


  # COMPUTE INTEGRAL


  # WRAP UP AND TRANSFORM OUTPUT AS REQUIRED
  if("cluster" %in% class(stop.cl)) {
    stopCluster(stop.cl)
  }
  return(rho)

  if(lower.tail) {
    return(data.frame(lower = ifelse(log, lower, exp(lower)),
                      upper = ifelse(log, upper, exp(upper))))
  } else {
    return(data.frame(lower = ifelse(log, log(-expm1(upper)), -expm1(upper)),
                      upper = ifelse(log, log(-expm1(lower)), -expm1(lower))))
  }
}

log.phi <- function(t, evals, ncps) {
  sum(
    complex(imaginary = ncps*t*evals)/complex(real = 1, imaginary = -2*t*evals) -
      0.5*log(complex(real = 1,imaginary = -2*t*evals))
  )
}




H <- function(eps, nu, resid.op.norm.bd) {
  exp(ifelse(eps <= nu / (4*resid.op.norm.bd),
             -0.5*(eps^2) / nu,
             0.5*nu / ((4*resid.op.norm.bd)^2) - eps / (4*resid.op.norm.bd)))
}
