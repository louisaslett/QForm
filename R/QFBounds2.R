QFBounds2 <- function(obs, evals, ncps, E_R, nu2, N, resid_operator_norm_bound, lower.tail = TRUE, log = FALSE,if.insufficient.eigs="trivial") {

  # Input Sanity Check
  if(!(if.insufficient.eigs%in%c("missing","trivial"))){stop("if.insufficient.eigs must be set to return ``trivial`` or return ``missing`` bounds ")}


  # Rescale by N to keep calculations within a range of reasonable precision
  obs <- obs/N
  E_R<-E_R/N
  nu2<-nu2/(N^2)
  evals<-evals/N
  resid_operator_norm_bound<-resid_operator_norm_bound/N

  # Calculate the mean and variance of the truncated part of the quadratic form
  # E_T <- sum(evals*(1+ncps))
  Var_T<-2*sum((evals^2)*(1+2*ncps))
  Var_R<-nu2/2


  # Find the machine precision we care about (the precision at which we could detect
  # changes in our objective funtion duing optimization)
  opt_tol <- min(c(machine_eps(abs(E_R)), machine_eps(abs(obs))))


  # INPUT CHECKS

  # Make sure that the number of eigen values we have allows us to get a sensible
  # result out of davies and squeeze the right interval in if not to within the
  # tolerance that optimize will work with

  # Ensure that the truncated distribution puts mass at obseravtion when epsilon is zero
  if(any(c(0,1) == SurvivalFunc(q = obs-E_R,
                                lambda = evals,
                                delta = ncps))) {
    stop("Insufficient eigenvalues: truncated distribution does not have tails that extend out to the observation given")
  }



  # This check not only insures that our bounding approach has a shot of yielding a tight bound, but also that
  # we will not drastically overshoot the support of the truncated distribution while optimizing our bounds.
  if(Var_T/Var_R < 100){
    stop("Insufficient eigenvalues: truncated distribution does not have at least 10x the scale of the residual distribution.")
  }


  # Find search interval:
  int_right<-uniroot(f=RemainderBoundSupportFinder,interval=c(nu2/(4*resid_operator_norm_bound),sqrt(Var_R)*100),resid_operator_norm_bound=resid_operator_norm_bound,nu2=nu2,extendInt = "downX",tol=1e-20)$root

  # This browser is for checking that the objective/optimization is working
  browser()

  lower <- as.numeric(optimize(QFBounds.ineq.lower,
                               interval = c(0, int_right),
                               obs = obs,
                               evals = evals,
                               ncps = ncps,
                               E_R = E_R,
                               nu2 = nu2,
                               resid_operator_norm_bound=resid_operator_norm_bound,
                               maximum = TRUE,
                               tol = opt_tol)$objective)

  upper <- as.numeric(optimize(QFBounds.ineq.upper,
                               interval = c(0, int_right),
                               obs = obs,
                               evals = evals,
                               ncps = ncps,
                               E_R = E_R,
                               nu2 = nu2,
                               resid_operator_norm_bound=resid_operator_norm_bound,
                               maximum = FALSE,
                               tol = opt_tol)$objective)

  # If upper in prob space is all geq 1, then we do not have enough eigenvalues to find a non-trivial bound

  if(upper>=0 && if.insufficient.eigs=="missing"){
    upper<-NA
    warning("Insufficient eigenvalues: returning NA for upper bound on CDF")
  }
  if(upper>=0 && if.insufficient.eigs=="trivial"){
    upper<-0
    warning("Insufficient eigenvalues: returning 0 as the trivial upper bound on the log of the CDF")
  }

  # If lower in log space is -.Machine$double.xmax, then we do not have enough eigenvalues to find a non-trivial bound
  if(lower==-.Machine$double.xmax && if.insufficient.eigs=="missing"){
    lower<-NA
    warning("Insufficient eigenvalues: returning NA for lower bound on CDF")
  }
  if(lower==-.Machine$double.xmax && if.insufficient.eigs=="trivial"){
    warning("Insufficient eigenvalues: returning -.Machine$double.xmax as the trivial lower bound on the log of the CDF")
  }


  if(lower.tail) {
    return(data.frame(lower = ifelse(log, lower, exp(lower)),
                      upper = ifelse(log, upper, exp(upper))))
  } else {
    return(data.frame(lower = ifelse(log, log(-expm1(upper)), -expm1(upper)),
                      upper = ifelse(log, log(-expm1(lower)), -expm1(lower))))
  }
}


QFBounds.ineq.lower <- function(eps, obs, evals, ncps, E_R, nu2,resid_operator_norm_bound) {

  survival_func_est<-SurvivalFunc(q = c(obs-E_R-eps),
                                  lambda = evals,
                                  delta = ncps)


  F_H <- - survival_func_est - H(eps,nu2,resid_operator_norm_bound)


  # F_H <=-1 happens whenever CDF(a-eps)==0 and in other cases.
  # The if else statement below constrains optimization to the support of the truncated distribution

  if(F_H <=-1){
    return(-.Machine$double.xmax)
  }else{
    return(log1p(F_H))
  }
}

QFBounds.ineq.upper <- function(eps, obs, evals, ncps, E_R, nu2,resid_operator_norm_bound) {

  survival_func_est<-SurvivalFunc(q = c(obs-E_R+eps),
                                  lambda = evals,
                                  delta = ncps)

  # If CDF(a+eps)==1 then set to CDF to 2 in order to constraint optimization to the support of the truncated distribution
  if(survival_func_est<=0){survival_func_est<--1}

  F_H <- -survival_func_est + H(eps,nu2,resid_operator_norm_bound)

  # Upper can return upper bounds that correspond to values larger than one.  If these are found, then they are
  # caught later in QFBounds2
  return(log1p(F_H))
}



RemainderBoundSupportFinder <- function(eps, resid_operator_norm_bound, nu2) {
  exp(0.5*nu2 / (abs(4*resid_operator_norm_bound)^2) - eps / abs(4*resid_operator_norm_bound)) - 1e-19
}



SurvivalFunc<-function(q,lambda,delta,lim=20000,acc=1e-12){

  result<-davies(q = q,
                 lambda = lambda,
                 delta = delta,
                 lim =lim,
                 acc = acc)

  if (result$ifault!=0) {
    if (result$ifault==1) {
      out_warn<-"requested accuracy could not be obtained"
    } else if (result$ifault==2) {
      out_warn<-"round-off error possibly significant"
    } else if (result$ifault==3) {
      out_warn<-"invalid parameters"
    } else if (result$ifault==4){
      out_warn<-"unable to locate integration parameters"
    } else {
      out_warn<-"davies produced unrcognized fault (not 1, 2, 3, or 4)"
    }
    warning(paste("fault in evaluating davies:",out_warn))
  }

  if(!is.numeric(result$Qq)){
    stop("davies returned non-numeric survival function estimate for Qq")
  }

  if(result$Qq<0 || result$Qq >1){
    stop("davies returned survival function estimate for Qq outside of [0,1]")
  }

  return(result$Qq)

}


H<-function(eps,nu2,resid_operator_norm_bound){
  exp(ifelse(eps <= nu2 / (4*resid_operator_norm_bound),
             -0.5*(eps^2) / nu2,
             0.5*nu2 / ((4*resid_operator_norm_bound)^2) - eps / (4*resid_operator_norm_bound)))
}

# ErrorHandleDavies<-function(q,lambda,delta,lim=20000,acc=1e-12)
#
#   result<-tryCatch({
#
#   }, warning = function(war){
#
#   }, error = function(err){
#
#   }, finally = {}
#   )
#
#   return(result)
# }

#tryCatch({survival_at_obs<-})

