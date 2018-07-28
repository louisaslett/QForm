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
QFIntBounds2 <- function(obs, evals, ncps, E_R, nu, resid.op.norm.bd, lower.tail = TRUE, log = FALSE, llin.approx = FALSE) {
  lower <- 0
  upper <- 1
  did.approx <- FALSE

  upper.res <- integrate(QFIntBounds.ineq, lower = 0, upper = Inf,
                         obs = obs, evals = evals, ncps = ncps, E_R = E_R, nu = nu, resid.op.norm.bd = resid.op.norm.bd, upper.bound = TRUE,
                         stop.on.error = FALSE, abs.tol = 0)
  if(is.environment(llin.approx) && (!(upper.res$value>0) || upper.res$abs.error/upper.res$value>.Machine$double.eps^0.25)) {
    # Failed to achieve accuracy, so do log-linear extrapolation
    upper.res <- list(message = "OK",
                      value = QFLLinApprox(llin.approx, obs, evals, ncps, E_R, nu, resid.op.norm.bd, upper.bound = TRUE))
    did.approx <- TRUE
  }
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

  # Easy to check if no chance of lower bound working ...
  if(is.environment(llin.approx) && !(SurvivalFuncInt(q = -obs + E_R,
                                                      lambda = -evals,
                                                      delta = ncps,
                                                      lim = 20000,
                                                      acc = 1e-12) > 0)) {
    lower.res <- list(message = "OK",
                      value = QFLLinApprox(llin.approx, obs, evals, ncps, E_R, nu, resid.op.norm.bd, upper.bound = FALSE))
    did.approx <- TRUE
  } else {
    lower.res <- integrate(QFIntBounds.ineq, lower = 0, upper = Inf,
                           obs = obs, evals = evals, ncps = ncps, E_R = E_R, nu = nu, resid.op.norm.bd = resid.op.norm.bd, upper.bound = FALSE,
                           stop.on.error = FALSE, abs.tol = 0)
    if(is.environment(llin.approx) && (!(lower.res$value>0) || lower.res$abs.error/lower.res$value>.Machine$double.eps^0.25)) {
      # Failed to achieve accuracy, so do log-linear extrapolation
      lower.res <- list(message = "OK",
                        value = QFLLinApprox(llin.approx, obs, evals, ncps, E_R, nu, resid.op.norm.bd, upper.bound = FALSE))
      did.approx <- TRUE
    }
  }
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

  data.frame(lower = lower, upper = upper, llin.approx = did.approx)
}

QFLLinApprox <- function(llin.approx, obs, evals, ncps, E_R, nu, resid.op.norm.bd, upper.bound) {
  if(upper.bound && exists("model.u", llin.approx)) {
    return(exp(predict(llin.approx$model.u, data.frame(x = obs))[[1]]))
  }
  if(!upper.bound && exists("model.l", llin.approx)) {
    return(exp(predict(llin.approx$model.l, data.frame(x = obs))[[1]]))
  }
  x <- c()
  y <- c()

  # Find a reasonable starting point for linear extrapolation
  k <- -14
  new.tail <- 0
  while(!(new.tail > 1e-10)) {
    end <- uniroot(function(obs, lambda, delta, E_R, k) {
      res <- SurvivalFuncInt(q = -obs + E_R,
                             lambda = lambda,
                             delta = delta,
                             lim = 20000,
                             acc = 1e-12)
      if(10^k < res && res < 10^(k+1)) {
        return(0)
      } else {
        return(res - 5.5*10^k)
      }
    }, c(0, obs), lambda = -evals, delta = ncps, E_R = E_R, k = k)$root
    k <- k+1
    new.tail <- integrate(QFIntBounds.ineq, lower = 0, upper = Inf,
                          obs = end, evals = evals, ncps = ncps, E_R = E_R, nu = nu, resid.op.norm.bd = resid.op.norm.bd, upper.bound = upper.bound,
                          stop.on.error = FALSE, abs.tol = 0)$value
  }
  for(q in seq(1e-10, 1e-12, length.out = 8)) {
    match.q <- uniroot(function(new.obs, q, fn, evals, ncps, E_R, nu, resid.op.norm.bd, upper.bound) {
      integrate(fn, lower = 0, upper = Inf,
                obs = new.obs, evals = evals, ncps = ncps, E_R = E_R, nu = nu, resid.op.norm.bd = resid.op.norm.bd, upper.bound = upper.bound,
                stop.on.error = FALSE, abs.tol = 0)$value -
        q
    }, c(obs, end), q = q, fn = QFIntBounds.ineq, evals = evals, ncps = ncps, E_R = E_R, nu = nu, resid.op.norm.bd = resid.op.norm.bd, upper.bound = upper.bound)
    x <- c(x, match.q$root)
    y <- c(y, match.q$f.root+q)
    end <- match.q$root
  }

  # # Push points out into tail until we lose the integral
  # new.obs <- new.obs + sign(new.obs)*sqrt(nu)
  # new.tail <- integrate(QFIntBounds.ineq, lower = 0, upper = Inf,
  #                       obs = new.obs, evals = evals, ncps = ncps, E_R = E_R, nu = nu, resid.op.norm.bd = resid.op.norm.bd, upper.bound = upper.bound,
  #                       stop.on.error = FALSE, abs.tol = 0)$value
  # while(new.tail > 0) {
  #   x <- c(x, new.obs)
  #   y <- c(y, new.tail)
  #   new.obs <- new.obs + sign(new.obs)*sqrt(nu)
  #   new.tail <- integrate(QFIntBounds.ineq, lower = 0, upper = Inf,
  #                         obs = new.obs, evals = evals, ncps = ncps, E_R = E_R, nu = nu, resid.op.norm.bd = resid.op.norm.bd, upper.bound = upper.bound,
  #                         stop.on.error = FALSE, abs.tol = 0)$value
  # }
  #
  # # Fill the rest with points the other side of starting point
  # while(length(x) < 8) {
  #   new.obs <- new.obs + sign(new.obs)*sqrt(nu)
  #   new.tail <- integrate(QFIntBounds.ineq, lower = 0, upper = Inf,
  #                         obs = new.obs, evals = evals, ncps = ncps, E_R = E_R, nu = nu, resid.op.norm.bd = resid.op.norm.bd, upper.bound = upper.bound,
  #                         stop.on.error = FALSE, abs.tol = 0)$value
  #   x <- c(x, new.obs)
  #   y <- c(y, new.tail)
  # }

  # Store this model for future log-linear approximations
  if(upper.bound) {
    llin.approx$model.u <- lm(log(y) ~ x, data.frame(x = x, y = y))
    if(summary(llin.approx$model.u)$r.squared < 0.995) {
      warning(glue("Log-linear approximation for upper bound may be dubious, R^2 = {summary(llin.approx)$r.squared}"))
    }
  } else {
    llin.approx$model.l <- lm(log(y) ~ x, data.frame(x = x, y = y))
    if(summary(llin.approx$model.l)$r.squared < 0.995) {
      warning(glue("Log-linear approximation for lower may be dubious, R^2 = {summary(llin.approx)$r.squared}"))
    }
  }

  # Recursively return prediction
  return(QFLLinApprox(llin.approx, obs, evals, ncps, E_R, nu, resid.op.norm.bd, upper.bound))
}

HInt <- function(z, nu, resid.op.norm.bd) {
  ifelse(z <= nu / (4*resid.op.norm.bd),
         (z/nu)*exp(-0.5*(z^2)/nu),
         exp(nu/(32*resid.op.norm.bd^2)-z/(4*resid.op.norm.bd))/(4*resid.op.norm.bd))
}

SurvivalFuncInt <- function(q, lambda, delta, lim = 20000, acc = 1e-12) {

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

QFIntBounds.ineq <- function(z, obs, evals, ncps, E_R, nu, resid.op.norm.bd, upper.bound = TRUE) {
  side_indicator <- 2*as.integer(upper.bound)-1

  HInt(z,nu,resid.op.norm.bd) *
    SurvivalFuncInt(q = -obs + E_R - side_indicator*z,
                    lambda = -evals,
                    delta = ncps,
                    lim = 20000,
                    acc = 1e-12)
}





################################################################################

QFIntBounds2.ryan <- function(obs, evals, ncps, E_R, nu, resid.op.norm.bd,log) {

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



HInt.ryan<- function(z, nu, resid.op.norm.bd) {
  ifelse(z <= nu / (4*resid.op.norm.bd),
         (z/nu)*exp(-0.5*(z^2)/ nu),exp(nu/(32*resid.op.norm.bd^2)-z/(4*resid.op.norm.bd))/(4*resid.op.norm.bd))
}

QFIntBounds.ineq.ryan<- function(z,obs, evals, ncps, E_R, nu, resid.op.norm.bd,upper.bound=T) {
  side_indicator<-2*as.integer(upper.bound)-1
  return(HInt(z,nu,resid.op.norm.bd)-HInt(z,nu,resid.op.norm.bd)*davies(q=obs-E_R+ side_indicator*z,lambda=evals,delta=ncps,lim = 20000, acc = 1e-12)$Qq)
}
