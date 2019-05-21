#' CDF/PDF of a Quadratic Form in Gaussians
#'
#' Returns the CDF and PDF for random variables \eqn{T_f}{T_f} of the form \deqn{T_f = \sum\limits_i f\left(\eta_i \right) \left(Z_i + \delta_i)^2}{T_f = \Sigma_i f (\eta_i) (Z_i + \delta_i)^2} where \eqn{Z_i \sim N(0,1).}{Z_i ~ N(0,1).}
#'
#' The returned function has three optional, logical arguments.  The first is a \code{density}, which when \code{TRUE}, prompts the function to evaluate the PDF rather than the CDF.  \code{density} defaults to FALSE.  \code{lower.tail} returns 1 minus the CDF when \code{TRUE} (not used if \code{density}==TRUE) and is highly recommended for those interested in the upper tail of \eqn{T_f}.  \code{log.p} returns the desired probabilities in log space.
#'
#' \code{n} is the number of sub-intervals used in the left-sided Reimann integral approximation of the Fourier transform carried out by \code{stats::fft}.  The default 2^16-1 should work for the vast majority of cases, but n may need to be increased to achieve accurate CDF estimation when \eqn{T_f} has many terms (when \code{f.eta} is long).
#'
#' Since \code{stats::fft} can only evaluate the CDF up to double precision, we extrapolate the tails of \eqn{T_f}.  QForm automatically detects the region where the estimated CDF begins to lose precision.  A log-linear function is used for tails that go out to infinity and a log-monomial functions is used for tails truncated at 0 (when all of the \code{f.eta} have the same sign).  These extrapolated tails, motivated by the form of the characteristic function, provide accurate approximations in most cases when compared against a quad-precision implementation (not yet included in \code{QForm}).
#'
#' Our current tail extrapolation scheme can struggle in cases where the target distribution is extremely skewed, since that leaves fewer points from the FFT with which to perform extrapolation on one of the tails.  We plan to address this in future versions by deploying a second FFT when needed.
#'
#' A note on unbounded densities: The density of \eqn{T_f} if guaranteed to be bounded if length(f.eta) > 2 and there is no trouble in density estimation posed by asymptotes.  In the length(f.eta)==2 case, if the two components of f.eta are of opposite signs, then the density of \eqn{T_f}{T_f} may have an asymptote at some value \eqn{t}{t}.  While the density in the neighborhood around that \eqn{t} should be accurately calculated, due to the FFT and spline interpolation approach used, the density reported at \eqn{t}{t} may be reported as some finite rather than as Inf.  In the length(f.eta)==1 case, QFGauss resorts to dchisq and the density at 0 is accurately reported as Inf.
#'
#'
#' @param f.eta vector; real-valued coefficients, \eqn{f(\eta_i)}, (may be positive or negative)
#' @param delta vector; mean shifts for each \eqn{Z_i}
#' @param n numeric; number of points at which to evaluate the characteristic function of \eqn{T_f}, must be odd (see Details).
#'
#' @return A function that evaluates the CDF or PDF of \eqn{T_f}.
#'
#' @examples
#' f.eta <- c(-12, -7, 5, 7, -9, 10, 8)
#' delta <- c(2, 10, -4, 3, 8, -5, -12)
#'
#' cdf <- QFGauss(f.eta, delta)
#'
#' # Inspect computed CDF
#' plot(cdf)
#'
#' # Plot computed CDF at desired points
#' x <- seq(-1500, 2000, len = 1e3)
#' plot(x,cdf(x),type="l") # ECDF
#'
#' # Compare computed CDF to empirical CDF of target distribution based on 10,000 samples
#' test(cdf)
#'
#' @export

QFGauss <- function(f.eta, delta = rep(0,length(f.eta)), n = 2^16-1){

  if(n%%2==0){stop("n must be odd")}

  if(length(f.eta)!=length(delta)){stop("delta does not have the same length as f.eta.")}

  if(all(f.eta==0)){stop("All f.eta are zero.")}

  f.eta.order <- order(abs(f.eta),decreasing=T)
  f.eta <- f.eta[f.eta.order]
  delta <- delta[f.eta.order]


  if(all(f.eta==f.eta[1])){

    df <- length(f.eta)
    C <- abs(f.eta[1])

      cdf.func <- function(q, density = FALSE, lower.tail = TRUE, log.p = FALSE){
        if(density){
          if(log.p){
            return( dchisq((2*(f.eta[1]>0)-1)*q/C,df,sum(delta^2),TRUE)-log(C) )
          }else{
            return( dchisq((2*(f.eta[1]>0)-1)*q/C,df,sum(delta^2))/C )
          }
        }else{
          if(f.eta[1] > 0){
            return( pchisq(q/C,df,sum(delta^2),lower.tail,log.p) )
          }else{
            return( pchisq(-q/C,df,sum(delta^2),!lower.tail,log.p) )

          }
        }
      }
      attr(cdf.func,"fft_used") <- FALSE
      attr(cdf.func,"f.eta") <- f.eta
      attr(cdf.func,"delta") <- delta
      attr(cdf.func,"mu") <- sum(f.eta*(1+delta^2))
      attr(cdf.func,"sigma") <- sum(2*(1+2*delta^2)*f.eta^2)
      attr(cdf.func,"tail.features") <- list("lambda.signs" = ifelse(f.eta[1] > 0,"pos","neg"),
                                             "extrapolation.point.l" = ifelse(f.eta[1] > 0,0,-Inf), #C*qchisq(1e-16,df,sum(delta^2))
                                             "extrapolation.point.r" = ifelse(f.eta[1] > 0,Inf,0),
                                             "a.l" = NULL,
                                             "b.l" = NULL,
                                             "a.r" = NULL,
                                             "b.r" = NULL)
      class(cdf.func) <- c("QFGaussCDF",class(cdf.func))
    return(cdf.func)
  }

  cdf <- calc.QFcdf(evals = f.eta, ncps = delta^2, n = n)
  cdf.func <- wrap.QFcdf(cdf)

  attr(cdf.func,"fft_used") <- TRUE
  attr(cdf.func,"f.eta") <- f.eta
  attr(cdf.func,"delta") <- delta
  attr(cdf.func,"mu") <- cdf$mu
  attr(cdf.func,"sigma") <- cdf$sigma
  attr(cdf.func,"tail.features") <- list("lambda.signs" = cdf$type,
                                         "extrapolation.point.l" = cdf$x[1],
                                         "extrapolation.point.r" = cdf$x[cdf$n],
                                         "a.l" = cdf$a.l,
                                         "b.l" = cdf$b.l,
                                         "a.r" = cdf$a.r,
                                         "b.r" = cdf$b.r)

  class(cdf.func) <- c("QFGaussCDF",class(cdf.func))
  cdf.func
}



#' Plot method for a QFGaussCDF object
#'
#' Plots the CDF computed by QFGauss.
#'
#' @param cdf a QFGaussCDF
#' @return There is nothing returned.
#'
#' @export
plot.QFGaussCDF <- function(cdf,...){
  if(class(cdf)[1]!="QFGaussCDF"){stop("cdf must be of class QFGaussCDF")}

  fft_used <- attr(cdf,"fft_used")
  f.eta <- attr(cdf,"f.eta")
  delta <- attr(cdf,"delta")
  mu <- attr(cdf,"mu")
  sigma <- attr(cdf,"sigma")

  tf <- attr(cdf,"tail.features")
  lambda.signs <- tf$lambda.signs
  ep.l <- tf$extrapolation.point.l
  ep.r <- tf$extrapolation.point.r
  a.l <- tf$a.l
  b.l <- tf$b.l
  a.r <- tf$a.r
  b.r <- tf$b.r


  if(any(is.na(c(a.l,b.l,a.r,b.r)))){stop("test cannot be run because tail extrapolation in QFGauss failed for at least one tail, resulting in NA value for a.l, b.l, a.r, or a.r")}

  if(!fft_used){
    df <- length(f.eta)
    C <- abs(f.eta[1])
    ep.l <- C*qchisq(1e-16,length(f.eta),sum(delta^2))
    ep.r <- C*qchisq(1e-16,length(f.eta),sum(delta^2),lower.tail = FALSE)
  }

  if(lambda.signs == "mixed"){
    x.max <-uniroot(function(z) {- cdf(z,lower.tail = F,log.p = T) / log(10) - 20},
                    lower = ep.r, upper = ep.r + 0.1*sigma,tol = .Machine$double.eps, extendInt = "upX")$root
    x.min <-uniroot(function(z) {- cdf(z,lower.tail = T,log.p = T) / log(10) - 20},
                    lower = ep.l - 0.1*sigma, upper = ep.l,tol = .Machine$double.eps, extendInt = "downX")$root
    x <- seq(x.min,x.max,len=1e5)
  }

  if(lambda.signs == "pos"){
    x.max <-uniroot(function(z) {- cdf(z,lower.tail = F,log.p = T) / log(10) - 20},
                    lower = ep.r, upper = ep.r + 0.1*sigma,tol = .Machine$double.eps, extendInt = "upX")$root
    x <- seq(0,x.max,len=1e5)
  }

  if(lambda.signs == "neg"){
    x.min <-uniroot(function(z) {- cdf(z,lower.tail = T,log.p = T) / log(10) - 20},
                    lower = ep.l - 0.1*sigma, upper = ep.l,tol = .Machine$double.eps, extendInt = "downX")$root
    x <- seq(x.min,0,len=1e3)
  }

  old.par <- par(no.readonly = T)
  plot(x, cdf(x), type = "l", lwd=1.5, ylab = expression(CDF)) # plot lower tail of CDF
  par(old.par)
}


#' Test function for a QFGaussCDF object
#'
#' Compares the CDF inferred by QFGauss to an empicical CDF.
#'
#' Four plots are produced.  The top-left plot overlays the CDF computed by QFGauss (in black) and the empirical CDF (in red) based on 10,000 samples.
#' The top-right plot shows the distance between the empirical and QFGauss-computed CDF and the corresponding ks.test p-value (two-sided alternative).
#' The two bottom plots compare the empirical CDF (in red) with the computed CDF (in black) in each tail.
#'
#' @param cdf a QFGaussCDF
#' @param n.samps number of draws from the target distribution with which to construct the empirical CDF
#' @return There is nothing returned.
#'
#' @export
TestQFGauss <- function(cdf, n.samps = 1e4){

  if(class(cdf)[1]!="QFGaussCDF"){stop("cdf must be of class QFGaussCDF")}

  fft_used <- attr(cdf,"fft_used")
  f.eta <- attr(cdf,"f.eta")
  delta <- attr(cdf,"delta")
  mu <- attr(cdf,"mu")
  sigma <- attr(cdf,"sigma")

  tf <- attr(cdf,"tail.features")
  lambda.signs <- tf$lambda.signs
  ep.l <- tf$extrapolation.point.l
  ep.r <- tf$extrapolation.point.r
  a.l <- tf$a.l
  b.l <- tf$b.l
  a.r <- tf$a.r
  b.r <- tf$b.r


  if(any(is.na(c(a.l,b.l,a.r,b.r)))){stop("test cannot be run because tail extrapolation in QFGauss failed for at least one tail, resulting in NA value for a.l, b.l, a.r, or a.r")}


  if(!fft_used){
    df <- length(f.eta)
    C <- abs(f.eta[1])
    ep.l <- C*qchisq(1e-16,length(f.eta),sum(delta^2))
    ep.r <- C*qchisq(1e-16,length(f.eta),sum(delta^2),lower.tail = FALSE)
  }

  if(lambda.signs == "mixed"){
    x.max <-uniroot(function(z) {- cdf(z,lower.tail = F,log.p = T) / log(10) - 20},
                    lower = ep.r, upper = ep.r + 0.1*sigma,tol = .Machine$double.eps, extendInt = "upX")$root
    x.min <-uniroot(function(z) {- cdf(z,lower.tail = T,log.p = T) / log(10) - 20},
                    lower = ep.l - 0.1*sigma, upper = ep.l,tol = .Machine$double.eps, extendInt = "downX")$root
    x <- seq(x.min,x.max,len=1e5)
  }

  if(lambda.signs == "pos"){
    x.max <-uniroot(function(z) {- cdf(z,lower.tail = F,log.p = T) / log(10) - 20},
                    lower = ep.r, upper = ep.r + 0.1*sigma,tol = .Machine$double.eps, extendInt = "upX")$root
    x <- seq(0,x.max,len=1e5)
  }

  if(lambda.signs == "neg"){
    x.min <-uniroot(function(z) {- cdf(z,lower.tail = T,log.p = T) / log(10) - 20},
                    lower = ep.l - 0.1*sigma, upper = ep.l,tol = .Machine$double.eps, extendInt = "downX")$root
    x <- seq(x.min,0,len=1e5)
  }

  old.par <- par(no.readonly = T)
  par(mfrow=c(2,2))
  samps <- c(colSums(f.eta * (matrix(rnorm(n.samps * length(f.eta)), nrow = length(f.eta)) + delta)^2))
  qecdf <- ecdf(samps)

  plot(x, cdf(x), type = "l", lwd=1.5, ylab = expression(CDF)) # plot lower tail of CDF
  lines(x, qecdf(x), lwd=1,col="red") # plot lower tail of CDF
  #legend("topleft",legend=c("QForm CDF","ECDF"),col=c("black","blue"),lty=c(1,3))

  plot(x,cdf(x)-qecdf(x), type="l",main=paste("ks.test p-value:",ks.test(samps,cdf)$p.value),lwd=1.5)
  abline(h=0)

  plot(x, -cdf(x, log.p = TRUE)/log(10), type = "l",
       ylab = expression(-log[10](CDF)) )
  lines(x,-log10(qecdf(x)),col="red")
  #legend("topright",legend=c("QForm CDF","ECDF"),col=c("black","blue"),lty=c(1,3))


  plot(x, -cdf(x, lower.tail = FALSE, log.p = TRUE)/log(10), type = "l",
       ylab = expression(-log[10](1 - CDF))) # plot upper tail of CDF
  lines(x,-log1p(-qecdf(x))/log(10),col="red")
  #legend("topleft",legend=c("QForm CDF","ECDF"),col=c("black","blue"),lty=c(1,3))

  par(old.par)
}







