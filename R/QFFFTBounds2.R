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
QFFFTBounds2 <- function(obs, evals, ncps, E_R, nu, resid.op.norm.bd, lower.tail = TRUE, log = FALSE) {
  lower <- 0
  upper <- 1

  upper.res <-

    integrate(QFIntBounds.ineq, lower = 0, upper = Inf,
                         obs = obs, evals = evals, ncps = ncps, E_R = E_R, nu = nu, resid.op.norm.bd = resid.op.norm.bd, upper.bound = TRUE,
                         stop.on.error = FALSE)
  if(upper.res$message != "OK") {
    warning("Numerical integration failed to achieve required accuracy for upper bound, returning trivial bound.")
    if(log)
      upper <- 0
    else
      upper <- 1
  } else {
    if(log) {
      if(lower.tail)
        upper <- log(upper.res$value)
      else
        lower <- log(1-upper.res$value)
    } else {
      if(lower.tail)
        upper <- upper.res$value
      else
        lower <- 1-upper.res$value
    }
  }

  lower.res <- integrate(QFIntBounds.ineq, lower = 0, upper = Inf,
                         obs = obs, evals = evals, ncps = ncps, E_R = E_R, nu = nu, resid.op.norm.bd = resid.op.norm.bd, upper.bound = FALSE,
                         stop.on.error = FALSE)
  if(lower.res$message != "OK") {
    warning("Numerical integration failed to achieve required accuracy for lower bound, returning trivial bound.")
    if(log)
      lower <- -.Machine$double.xmax
    else
      lower <- 0
  } else {
    if(log) {
      if(lower.tail)
        lower <- log(lower.res$value)
      else
        upper <- log(1-lower.res$value)
    } else {
      if(lower.tail)
        lower <- lower.res$value
      else
        upper <- 1-lower.res$value
    }
  }

  data.frame(lower = lower, upper = upper)
}


HFFT <- function(z, nu, resid.op.norm.bd) {
  z<-ifelse(z>0,0,z)
  exp(ifelse(z <= nu / (4*resid.op.norm.bd),
             -0.5*(z^2) / nu,
             0.5*nu / ((4*resid.op.norm.bd)^2) - z / (4*resid.op.norm.bd)))
}



QFIntBounds.ineq <- function(z, obs, evals, ncps, E_R, nu, resid.op.norm.bd, upper.bound = TRUE) {
  side_indicator <- 2*as.integer(upper.bound)-1

  HInt(z,nu,resid.op.norm.bd) *
    SurvivalFuncInt(q = -obs + E_R - side_indicator*z,
                    lambda = -evals,
                    delta = ncps,
                    lim = 20000,
                    acc = 1e-12)
}



