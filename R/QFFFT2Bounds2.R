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
QFFFT2Bounds2 <- function(obs, evals, ncps, E_R, nu, N, resid.op.norm.bd, if.insuff.eigs = "trivial", lower.tail = TRUE, log = FALSE, parallel = FALSE, extended.acc = FALSE) {

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
    t.pdf <- Tpdf(evals, ncps, qfft.apply = qfft.apply)
  }

  if(FALSE) { # Diagnostic
    # Full density plot
    plot(t.pdf$x, -log(t.pdf$d), type = "l")
    abline(v = t.pdf$x[t.pdf$index.l])
    abline(v = t.pdf$x[t.pdf$index.r])
    abline(t.pdf$a.l, t.pdf$b.l, col = "blue")
    abline(t.pdf$a.r, t.pdf$b.r, col = "orange")
    lines(t.pdf$x, Vectorize(function(z) {
      if(z <= t.pdf$x[t.pdf$index.l]) {
        return(t.pdf$a.l + t.pdf$b.l*z)
      }
      if(z >= t.pdf$x[t.pdf$index.r]) {
        return(t.pdf$a.r + t.pdf$b.r*z)
      }
      return(-log(t.pdf$d[which.min(abs(t.pdf$x-z))]))
    })(t.pdf$x), col = "green")

    # Zoom left extrapolation
    tmp <- (t.pdf$index.l-20):(t.pdf$index.l+20)
    plot(t.pdf$x[tmp], -log(t.pdf$d[tmp]), type = "l")
    lines(t.pdf$x[tmp], Vectorize(function(z) {
      if(z <= t.pdf$x[t.pdf$index.l]) {
        return(t.pdf$a.l + t.pdf$b.l*z)
      }
      if(z >= t.pdf$x[t.pdf$index.r]) {
        return(t.pdf$a.r + t.pdf$b.r*z)
      }
      return(-log(t.pdf$d[which.min(abs(t.pdf$x-z))]))
    })(t.pdf$x[tmp]), col = "green")
    abline(v = t.pdf$x[t.pdf$index.l])

    # Zoom right extrapolation
    tmp <- (t.pdf$index.r-20):(t.pdf$index.r+20)
    plot(t.pdf$x[tmp], -log(t.pdf$d[tmp]), type = "l")
    lines(t.pdf$x[tmp], Vectorize(function(z) {
      if(z <= t.pdf$x[t.pdf$index.l]) {
        return(t.pdf$a.l + t.pdf$b.l*z)
      }
      if(z >= t.pdf$x[t.pdf$index.r]) {
        return(t.pdf$a.r + t.pdf$b.r*z)
      }
      return(-log(t.pdf$d[which.min(abs(t.pdf$x-z))]))
    })(t.pdf$x[tmp]), col = "green")
    abline(v = t.pdf$x[t.pdf$index.r])
  }


  # COMPUTE INTEGRAL


  # WRAP UP AND TRANSFORM OUTPUT AS REQUIRED
  if("cluster" %in% class(stop.cl)) {
    stopCluster(stop.cl)
  }
  return(t.pdf)

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
