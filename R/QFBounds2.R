QFBounds2 <- function(obs, evals, ncps, E_R, nu2, N, lower.tail = TRUE, log = FALSE) {
  obs <- obs/N

  # What is the machine precision we care about
  opt_tol <- min(c(machine_eps(abs(E_R)), machine_eps(abs(obs))))

  # Make sure that the number of eigen values we have allows us to get a sensible
  # result out of davies and squeeze the right interval in if not to within the
  # tolerance that optimize will work with
  if(any(c(0,1) == davies(q = obs,
                          lambda = evals/N,
                          delta = ncps,
                          acc = 1e-12)$Qq)) {
    stop("Too few eigenvalues for truncated distribution to place probability mass at the required quantile.")
  }

  int_right <- max(abs(obs), 10)
  while(int_right >= opt_tol && any(c(0,1) == davies(q = obs,
                                                     lambda = evals/N,
                                                     delta = ncps,
                                                     acc = 1e-12)$Qq)) {
    int_right <- int_right * 0.5
  }

  lower <- as.numeric(optimize(QFBounds.ineq.lower,
                               interval = c(0, int_right),
                               obs = obs,
                               evals = evals/N,
                               ncps = ncps,
                               E_R = E_R,
                               nu2 = nu2,
                               maximum = TRUE,
                               tol = opt_tol)$objective)

  upper <- as.numeric(optimize(QFBounds.ineq.upper,
                               interval = c(0, int_right),
                               obs = obs,
                               evals = evals/N,
                               ncps = ncps,
                               E_R = E_R,
                               nu2 = nu2,
                               maximum = FALSE,
                               tol = opt_tol)$objective)
  upper <- min(upper, 0)

  if(isTRUE(all.equal(lower, upper))) {
    lower <- -Inf
    upper <- 0
    warning("Double precision floating point lack sufficient precision to accurately optimize bounds.  Returning trivial bounds.")
  }

  if(lower.tail) {
    return(data.frame(lower = ifelse(log, lower, exp(lower)),
                      upper = ifelse(log, upper, exp(upper))))
  } else {
    return(data.frame(lower = ifelse(log, log(-expm1(upper)), -expm1(upper)),
                      upper = ifelse(log, log(-expm1(lower)), -expm1(lower))))
  }
}

QFBounds.ineq.lower <- function(eps, obs, evals, ncps, E_R, nu2) {
  a <- eps - E_R
  resid_operator_norm_bound <- abs(evals[length(evals)])
  F_H <- -davies(q = c(obs-eps),
                 lambda = evals,
                 delta = ncps,
                 acc = 1e-12)$Qq -
    exp(ifelse(a <= nu2 / abs(4*resid_operator_norm_bound),
               -0.5*(a^2) / nu2,
               0.5*nu2 / (abs(4*resid_operator_norm_bound)^2) - a / abs(4*resid_operator_norm_bound)))
  if(F_H <= -1 || F_H == 0) {
    return(-.Machine$double.xmax)
  } else {
    return(log1p(F_H))
  }
}

QFBounds.ineq.upper <- function(eps, obs, evals, ncps, E_R, nu2) {
  a <- eps - E_R
  resid_operator_norm_bound <- abs(evals[length(evals)])
  F_H <- -davies(q = c(obs+eps),
                 lambda = evals,
                 delta = ncps,
                 acc = 1e-12)$Qq +
    exp(ifelse(a <= nu2 / abs(4*resid_operator_norm_bound),
               -0.5*(a^2) / nu2,
               0.5*nu2 / (abs(4*resid_operator_norm_bound)^2) - a / abs(4*resid_operator_norm_bound)))
  if(F_H <= -1) {
    return(.Machine$double.xmax)
  } else {
    return(log1p(F_H))
  }
}
