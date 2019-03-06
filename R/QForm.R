#'  Quadratic Form Bounds
#'
#' Compute upper and lower bounds on the CDF of a quadratic form in normal random variables.
#'
#' Gory details....
#'
#' @param obs vector; observed values of the quadratic form for which upper and lower bounds on the CDF
#' @param qf.cdf QForm cdf object; output of QFcdf for some set of coefficients and ncps
#' @param fun character string; function to be applied to coefficients.  Current options: "identity","power","exponential"
#' @param fun.args list; list of function arguments required for the corresponding function type provided in fun
#' @param lower.tail logical; if \code{TRUE} (default), probability is \eqn{P(y^T M y \le obs)}, otherwise \eqn{P(y^T M y > obs)}
#' @param log logical; if \code{TRUE}, probability \eqn{p} is given as \eqn{log(p)}
#'
#' @return A data frame containing the variables lower and upper which provide the bounds on the CDF of the quadratic form.
#'
#' @examples
#' # Some code here which runs a self-contained example
#'
#' @export

QForm <- function(lambda, delta = rep(0,length(lambda)), n = 2^16-1){

  if(n%%2==0){stop("n must be odd")}
  cdf <- calc.QFcdf(evals = lambda, ncps = delta^2, n = n)
  cdf.func <- wrap.QFcdf(cdf)
  attr(cdf.func,"tail.features") <- list("lambda.signs" = cdf$type,
                                         "extrapolation.point.l" = cdf$x[1],
                                         "extrapolation.point.r" = cdf$x[cdf$n],
                                         "a.l" = cdf$a.l,
                                         "b.l" = cdf$b.l,
                                         "a.r" = cdf$a.r,
                                         "b.r" = cdf$b.r)
  cdf.func
}





