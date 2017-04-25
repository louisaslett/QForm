QFBounds3 <- function(obs, evals, ncps, E_R, nu2, N, lower.tail = TRUE, log = FALSE, if.insufficient.eigs="one-sided") {

  if(!(if.insufficient.eigs%in%c("one-sided","abort"))){stop("if.insufficient.eigs must be set to ``one-sided`` or ``abort``")}

  # Make sure that the number of eigen values we have allows us to get a sensible
  # result out of davies and squeeze the right interval in if not to within the
  # tolerance that optimize will work with
  if(any(c(0,1) == davies(q = obs,
                          lambda = evals/N,
                          delta = ncps,
                          acc = 1e-12)$Qq)) {
    stop("Too few eigenvalues for truncated distribution to place probability mass at the required quantile.")
  }

  # Calculate the mean and variance of the truncated part of the quadratic form
  E_T <- sum(evals*(1+ncps))
  Var_T<-2*sum((evals^2)*(1+2*ncps))
  Var_R<-nu2/2

  # This check not only insures that our bounding approach has a shot of yielding a tight bound, but also that
  # we will not drastically overshoot the support of the truncated distribution while optimizing our bounds.
  if(Var_T/Var_R < 100){
    stop("Too few eigenvalues for truncated distribution to have at least 10x the scale of the residual distribution.")
  }

  resid_operator_norm_bound<-min(abs(evals))

  # Find search interval:
  int_right<-uniroot(f=RemainderBoundSupportFinder,interval=c(nu2/(4*resid_operator_norm_bound),sqrt(Var_R)*100)+E_R,f.upper=1,f.lower=- 1,resid_operator_norm_bound=resid_operator_norm_bound,ncps=ncps,E_R=E_R,nu2=nu2,extendInt = "downX",tol=1e-20)$root


  # Rescale observation to keep computation in a more numerically precise range
  obs <- obs/N

  # What is the machine precision we care about
  opt_tol <- min(c(machine_eps(abs(E_R)), machine_eps(abs(obs))))


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

  # if(isTRUE(all.equal(lower, upper))) {
  #   lower <- -Inf
  #   upper <- 0
  #   warning("Double precision floating point lack sufficient precision to accurately optimize bounds.  Returning trivial bounds.")
  # }

  if(lower.tail) {
    return(data.frame(lower = ifelse(log, lower, exp(lower)),
                      upper = ifelse(log, upper, exp(upper))))
  } else {
    return(data.frame(lower = ifelse(log, log(-expm1(upper)), -expm1(upper)),
                      upper = ifelse(log, log(-expm1(lower)), -expm1(lower))))
  }
}







test_func<-function(x ,evals ,ncps){
y<-davies(q = x,
                  lambda = evals,
                  delta = ncps,
                  acc = 1e-20,lim = 40000)$Qq
if(y<=1e-19){y<--1}
return(y)
}

load("/Users/ryan/Documents/lab/thesis/Chapter3/chapter3_R_code/abo_min_mat_and_p_for_concentration_bounds.RData")
require(CompQuadForm)
require(RSpectra)
evals<-eigs(abo_min_mat,k=25)$values
evals<-evals[order(abs(evals),decreasing = T)]
test<-uniroot(f=test_func,interval=c(0,35e6),f.upper=1,f.lower=-1,evals=evals,ncps=rep(0,length(evals)),extendInt = "downX",trace = ,tol=1e-18)
test
test_func(test$root,evals,rep(0,length(evals)))

test_func(test$root-1,evals,rep(0,length(evals)))

# need to build in mean as lower bound and potentially 10 sigma as upper bound?
# then do the same for the other tail

# then optimize both sides.


concentration_bounder <- function(x, evals, ncps, E_R, nu2) {
  a <- x - E_R
  resid_operator_norm_bound <- abs(evals[length(evals)])
  0.5*nu2 / (abs(4*resid_operator_norm_bound)^2) - a / abs(4*resid_operator_norm_bound) - 1e-19
}


nu2<-4 * sum(abo_min_mat^2)

# setting the upper limit to nu2*40
test<-uniroot(f=concentration_bounder,interval=c(nu2/(4*abs(evals[length(evals)])),sqrt((nu2/2))*100),f.upper=1,f.lower=- 1,evals=evals,ncps=rep(0,length(evals)),E_R=0,nu2=nu2,extendInt = "downX",tol=1e-20)
concentration_bounder(test$root+10,evals=evals,ncps=rep(0,length(evals)),E_R=0,nu2=nu2)
test$root









