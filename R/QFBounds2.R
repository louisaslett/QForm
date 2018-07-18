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
QFBounds2 <- function(obs, evals, ncps, E_R, nu, N, resid.op.norm.bd, if.insuff.eigs = "trivial", lower.tail = TRUE, log = FALSE) {

  # PRELIMINARY CALCULATIONS

  # Input Sanity Check
  if(!(if.insuff.eigs %in% c("missing","trivial"))) {
    stop("if.insuff.eigs must be set to return \"trivial\" or return \"missing\" bounds.")
  }

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

  # Find the machine precision we care about (the precision at which we could detect
  # changes in our objective funtion duing optimization)
  opt_tol <- min(c(machine_eps(abs(E_R)), machine_eps(abs(obs))))


  # CHECK FOR SUFFICIENT EIGENVALUES

  # We need to make sure that the number of eigen values we have allows us to get a sensible
  # result out of davies and squeeze the right interval in if not to within the
  # tolerance that optimize can work with

  OK_to_optimize <- TRUE

  # Ensure that the truncated distribution puts mass at obseravtion when epsilon is zero
  if(any(c(0,1) == SurvivalFunc(q = obs-E_R,
                                lambda = evals,
                                delta = ncps))) {

    OK_to_optimize <- FALSE

    if(if.insuff.eigs == "trivial") {
      warning("Insufficient eigenvalues: truncated distribution does not have tails that extend out to the observation given: returning trivial 0,1 bounds.")
      upper <- 0
      lower <- -.Machine$double.xmax
    } else {
      warning("Insufficient eigenvalues: truncated distribution does not have tails that extend out to the observation given: returning missing bounds.")
      upper <- NA
      lower <- NA
    }
  }

  # This check not only ensures that our bounding approach has a shot of yielding a tight bound, but also that
  # we will not drastically overshoot the support of the truncated distribution while optimizing our bounds.
  if(Var_T/Var_R < 100) {

    OK_to_optimize <- FALSE

    if(if.insuff.eigs == "trivial") {
      warning("Insufficient eigenvalues: truncated distribution does not have at least 10x the scale of the residual distribution: returning trivial 0,1 bounds.")
      upper <- 0
      lower <- -.Machine$double.xmax
    } else {
      warning("Insufficient eigenvalues: truncated distribution does not have at least 10x the scale of the residual distribution: returning missing bounds.")
      upper <- NA
      lower <- NA
    }
  }


  # BEGIN OPTIMIZATION

  if(OK_to_optimize) {
    # Find search interval:
    int_right <- uniroot(f = RemainderBoundSupportFinder,
                         interval = c(nu/(4*resid.op.norm.bd), sqrt(Var_R)*100),
                         resid.op.norm.bd = resid.op.norm.bd,
                         nu = nu,
                         extendInt = "downX",
                         tol = 1e-20)$root

    # This browser is for checking that the objective/optimization is working
    #browser()

    lower <- as.numeric(optimize(QFBounds.ineq.lower,
                                 interval = c(0, int_right),
                                 obs = obs,
                                 evals = evals,
                                 ncps = ncps,
                                 E_R = E_R,
                                 nu = nu,
                                 resid.op.norm.bd = resid.op.norm.bd,
                                 maximum = TRUE,
                                 tol = opt_tol)$objective)

    upper <- as.numeric(optimize(QFBounds.ineq.upper,
                                 interval = c(0, int_right),
                                 obs = obs,
                                 evals = evals,
                                 ncps = ncps,
                                 E_R = E_R,
                                 nu = nu,
                                 resid.op.norm.bd = resid.op.norm.bd,
                                 maximum = FALSE,
                                 tol = opt_tol)$objective)

    # If upper in prob space is all geq 1, then we do not have enough eigenvalues to find a non-trivial bound

    if(upper >= 0 && if.insuff.eigs == "missing") {
      upper <- NA
      warning("Insufficient eigenvalues: returning NA for upper bound on CDF.")
    }
    if(upper >= 0 && if.insuff.eigs == "trivial") {
      upper <- 0
      warning("Insufficient eigenvalues: returning 0 as the trivial upper bound on the log of the CDF.")
    }

    # If lower in log space is -.Machine$double.xmax, then we do not have enough eigenvalues to find a non-trivial bound
    if(lower == -.Machine$double.xmax && if.insuff.eigs == "missing") {
      lower <- NA
      warning("Insufficient eigenvalues: returning NA for lower bound on CDF.")
    }
    if(lower == -.Machine$double.xmax && if.insuff.eigs == "trivial") {
      warning("Insufficient eigenvalues: returning -.Machine$double.xmax as the trivial lower bound on the log of the CDF.")
    }
  }


  # TRANSFORM OUTPUT AS REQUIRED

  if(lower.tail) {
    return(data.frame(lower = ifelse(log, lower, exp(lower)),
                      upper = ifelse(log, upper, exp(upper))))
  } else {
    return(data.frame(lower = ifelse(log, log(-expm1(upper)), -expm1(upper)),
                      upper = ifelse(log, log(-expm1(lower)), -expm1(lower))))
  }
}


QFBounds.ineq.lower <- function(eps, obs, evals, ncps, E_R, nu, resid.op.norm.bd) {

  survival_func_est <- SurvivalFunc(q = c(obs-E_R-eps),
                                    lambda = evals,
                                    delta = ncps)


  F_H <- -survival_func_est - H(eps, nu, resid.op.norm.bd)

  # F_H <=-1 happens whenever CDF(a-eps)==0 and in other cases.
  # The if else statement below constrains optimization to the support of the truncated distribution

  if(F_H <= -1) {
    return(-.Machine$double.xmax)
  } else {
    return(log1p(F_H))
  }
}


QFBounds.ineq.upper <- function(eps, obs, evals, ncps, E_R, nu, resid.op.norm.bd) {

  survival_func_est <- SurvivalFunc(q = c(obs-E_R+eps),
                                    lambda = evals,
                                    delta = ncps)

  # If CDF(a+eps)==1 then set to CDF to 2 in order to constraint optimization to the support of the truncated distribution
  if(survival_func_est <= 0) {
    survival_func_est <- -1
  }

  F_H <- -survival_func_est + H(eps, nu, resid.op.norm.bd)

  # Upper can return upper bounds that correspond to values larger than one.  If these are found, then they are
  # caught later in QFBounds2
  return(log1p(F_H))
}


RemainderBoundSupportFinder <- function(eps, resid.op.norm.bd, nu) {
  exp(0.5*nu / (abs(4*resid.op.norm.bd)^2) - eps / abs(4*resid.op.norm.bd)) - 1e-19
}


SurvivalFunc <- function(q, lambda, delta, lim = 20000, acc = 1e-12) {

  result <- davies(q = q,
                   lambda = lambda,
                   delta = delta,
                   lim = lim,
                   acc = acc)

  if(result$ifault != 0) {
    if(result$ifault == 1) {
      out_warn <- "requested accuracy could not be obtained."
    } else if(result$ifault == 2) {
      out_warn <- "round-off error possibly significant."
    } else if(result$ifault == 3) {
      out_warn <- "invalid parameters."
    } else if(result$ifault == 4){
      out_warn <- "unable to locate integration parameters."
    } else {
      out_warn <- "davies produced unrcognized fault (not 1, 2, 3, or 4)."
    }
    warning(paste("fault in evaluating davies:", out_warn))
  }

  if(!is.numeric(result$Qq)) {
    stop("davies returned non-numeric survival function estimate for Qq.")
  }

  if(result$Qq < 0 || result$Qq > 1) {
    stop("davies returned survival function estimate for Qq outside of [0,1].")
  }

  return(result$Qq)
}


H <- function(eps, nu, resid.op.norm.bd) {
  exp(ifelse(eps <= nu / (4*resid.op.norm.bd),
             -0.5*(eps^2) / nu,
             0.5*nu / ((4*resid.op.norm.bd)^2) - eps / (4*resid.op.norm.bd)))
}
