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
QFIntBounds2 <- function(obs, evals, ncps, E_R, nu, resid.op.norm.bd,log) {

  require(CompQuadForm)

  #browser()

  upper<-integrate(QFIntBounds.ineq, lower=0, upper=Inf,obs=obs, evals=evals, ncps=ncps, E_R=E_R, nu=nu, resid.op.norm.bd,upper.bound=T)$value

  lower<-integrate(QFIntBounds.ineq, lower=0, upper=Inf,obs=obs, evals=evals, ncps=ncps, E_R=E_R, nu=nu, resid.op.norm.bd,upper.bound=F)$value
  if(log==T){return(c(log(lower),log(upper)))}else{return(c(lower,upper))}

  # TRANSFORM OUTPUT AS REQUIRED

  # if(lower.tail) {
  #   return(data.frame(lower = ifelse(log, lower, exp(lower)),
  #                     upper = ifelse(log, upper, exp(upper))))
  # } else {
  #   return(data.frame(lower = ifelse(log, log(-expm1(upper)), -expm1(upper)),
  #                     upper = ifelse(log, log(-expm1(lower)), -expm1(lower))))
  # }
}



HInt<- function(z, nu, resid.op.norm.bd) {
  ifelse(z <= nu / (4*resid.op.norm.bd),
         (z/nu)*exp(-0.5*(z^2)/ nu),exp(nu/(32*resid.op.norm.bd^2)-z/(4*resid.op.norm.bd))/(4*resid.op.norm.bd))
}

QFIntBounds.ineq<- function(z,obs, evals, ncps, E_R, nu, resid.op.norm.bd,upper.bound=T) {
  side_indicator<-2*as.integer(upper.bound)-1
  return(HInt(z,nu,resid.op.norm.bd)-HInt(z,nu,resid.op.norm.bd)*davies(q=obs-E_R+ side_indicator*z,lambda=evals,delta=ncps,lim = 20000, acc = 1e-12)$Qq)
}


