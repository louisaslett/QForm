

QFIntegrate2<-function(q, t.pdf, log.conc.ineq.left,log.conc.ineq.right, L, d, psi, nu){
  
  upper.components <- lower.components <- rep(-Inf,4)
  
  if(q<t.pdf$x[t.pdf$index.l]){
    
    upper.components[1] <- ExpInt(t.pdf$a.l,t.pdf$b.l,-Inf,q) # upper bound on CDF
    upper.components[2] <- AnalyticBoundInt(t.pdf$a.l,t.pdf$b.l,q,t.pdf$x[t.pdf$index.l]) # gap from q to x[index.l] * bound
    upper.components[3] <- AnalyticBoundInt(t.pdf$a.r,t.pdf$b.r,t.pdf$x[t.pdf$index.r],Inf) # x[index.r] to Inf * bound
    upper.components[4] <- boole(t.pdf$x[t.pdf$index.l], t.pdf$x[t.pdf$index.l:t.pdf$index.r], exp(log(t.pdf$y[t.pdf$index.l:t.pdf$index.r])+log.conc.ineq.right(t.pdf$x[t.pdf$index.l:t.pdf$index.r],q,L,d,psi,nu)), int.to.right = T) # Middle * bound
    
    lower.components[1] <- boole(t.pdf$x[t.pdf$index.l],t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r],int.to.right = T,one.minus = T) # Middle
    lower.components[2] <- ExpInt(t.pdf$a.l,t.pdf$b.l,q,t.pdf$x[t.pdf$index.l]) # Left #
    lower.components[3] <- ExpInt(t.pdf$a.r,t.pdf$b.r,t.pdf$x[t.pdf$index.r],Inf) # Right  # WE COULD REPLEACE THESE THREE TERMS WITH A LOWER BOUND ON CDF ON LEFT TAIL
    lower.components[4] <- AnalyticBoundInt(t.pdf$a.l,t.pdf$b.l,-Inf,q,right=F) # -Inf to q * bound
    
    }
  
  if(q>t.pdf$x[t.pdf$index.r]){
    
    upper.components[1] <- ExpInt(t.pdf$a.l,t.pdf$b.l,-Inf,t.pdf$x[t.pdf$index.l]) # Left
    upper.components[2] <- boole(t.pdf$x[t.pdf$index.l],t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r]) # Middle
    upper.components[3] <- ExpInt(t.pdf$a.r,t.pdf$b.r,t.pdf$x[t.pdf$index.r],q) # Right
    upper.components[4] <- AnalyticBoundInt(t.pdf$a.r,t.pdf$b.r,q,Inf) # q to Inf * bound
    
    lower.components[1] <- boole(t.pdf$x[t.pdf$index.l],t.pdf$x[t.pdf$index.l:t.pdf$index.r],exp(log(t.pdf$y[t.pdf$index.l:t.pdf$index.r])+log.conc.ineq.left(t.pdf$x[t.pdf$index.l:t.pdf$index.r],q,L,d,psi,nu)),int.to.right = T,one.minus = T) # Middle* bound
    lower.components[2] <- ExpInt(t.pdf$a.r,t.pdf$b.r,q,Inf) # Right
    lower.components[3] <- AnalyticBoundInt(t.pdf$a.r,t.pdf$b.r,t.pdf$x[t.pdf$index.r],q,right=F) # Right * bound
    lower.components[4] <- AnalyticBoundInt(t.pdf$a.l,t.pdf$b.l,-Inf,t.pdf$x[t.pdf$index.l],right=F) # Left * bound
    
  }
  
  if( t.pdf$x[t.pdf$index.l]<=q & q<=t.pdf$x[t.pdf$index.r] ){
    
    upper.components[1] <- ExpInt(t.pdf$a.l,t.pdf$b.l,-Inf,t.pdf$x[t.pdf$index.l]) # Left
    upper.components[2] <- boole(q,t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r],int.to.right = F) # Middle Left
    upper.components[3] <- boole(q,t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r]*exp(log.conc.ineq.right(t.pdf$x[t.pdf$index.l:t.pdf$index.r],q,L,d,psi,nu)),int.to.right = T) # Middle Right * Bound
    upper.components[4] <- AnalyticBoundInt(t.pdf$a.r,t.pdf$b.r,t.pdf$x[t.pdf$index.r],Inf) # Right * Bound
    
    lower.components[1] <- boole(q,t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r],int.to.right = T,one.minus = T) # Middle Right
    lower.components[2] <- ExpInt(t.pdf$a.r,t.pdf$b.r,t.pdf$x[t.pdf$index.r],Inf) # Right
    lower.components[3] <- boole(q,t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r]*exp(log.conc.ineq.left(t.pdf$x[t.pdf$index.l:t.pdf$index.r],q,L,d,psi,nu)),int.to.right = F) # Middle Left
    lower.components[4] <- AnalyticBoundInt(t.pdf$a.l,t.pdf$b.l,-Inf,t.pdf$x[t.pdf$index.l],right=F) # Left * bound
    
  }
  
  
  u.sign.vec <- rep(1,4)
  u.norm.const <- max(upper.components)
  u.res <- neumaierSum(u.sign.vec*exp(upper.components-u.norm.const))*exp(u.norm.const)
  if(u.res < 0){u.res <- 0}
  if(u.res > 1){u.res <- 1}
  
  l.sign.vec <- c(1,-1,-1,-1)
  l.norm.const <- max(lower.components)
  l.res <- neumaierSum(l.sign.vec*exp(lower.components-l.norm.const))*exp(l.norm.const)
  if(l.res < 0){l.res <- 0}
  if(l.res > 1){l.res <- 1}
  
  c(l.res,u.res)
  
}




# 
# 
# QFIntegrate<-function(q, t.pdf, log.conc.ineq.left,log.conc.ineq.right, L, d, psi, nu){
# 
#   upper.components<-lower.components<-rep(0,4)
#   
#   if(q<t.pdf$x[t.pdf$index.l]){
# 
#     upper.components[1] <- ExpInt(t.pdf$a.l,t.pdf$b.l,-Inf,q) # upper bound on CDF
#     upper.components[2] <- ExpBoundInt(t.pdf$a.l,t.pdf$b.l,q,t.pdf$x[t.pdf$index.l],log.conc.ineq.right) # gap from q to x[index.l] * bound
#     upper.components[3] <- ExpBoundInt(t.pdf$a.r,t.pdf$b.r,t.pdf$x[t.pdf$index.r],Inf,log.conc.ineq.right) # x[index.r] to Inf * bound
#     upper.components[4] <- boole(t.pdf$x[t.pdf$index.l], t.pdf$x[t.pdf$index.l:t.pdf$index.r], exp(log(t.pdf$y[t.pdf$index.l:t.pdf$index.r])+log.conc.ineq.right(t.pdf$x[t.pdf$index.l:t.pdf$index.r],q,L,d,psi,nu)), int.to.right = T) # Middle * bound
# 
#     lower.components[1] <- boole(t.pdf$x[t.pdf$index.l],t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r],int.to.right = T,one.minus = T) # Middle
#     lower.components[2] <- -ExpInt(t.pdf$a.l,t.pdf$b.l,q,t.pdf$x[t.pdf$index.l]) # Left #
#     lower.components[3] <- -ExpInt(t.pdf$a.r,t.pdf$b.r,t.pdf$x[t.pdf$index.r],Inf) # Right  # WE COULD REPLEACE THESE THREE TERMS WITH A LOWER BOUND ON CDF ON LEFT TAIL
#     lower.components[4] <- -ExpBoundInt(t.pdf$a.l,t.pdf$b.l,-Inf,q,log.conc.ineq.left) # -Inf to q * bound
#     
#   }
# 
#   if(q>t.pdf$x[t.pdf$index.r]){
#     
#     upper.components[1] <- ExpInt(t.pdf$a.l,t.pdf$b.l,-Inf,t.pdf$x[t.pdf$index.l]) # Left
#     upper.components[2] <- boole(t.pdf$x[t.pdf$index.l],t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r],int.to.right = T) # Middle
#     upper.components[3] <- ExpInt(t.pdf$a.r,t.pdf$b.r,t.pdf$x[t.pdf$index.r],q) # Right
#     upper.components[4] <- ExpBoundInt(t.pdf$a.r,t.pdf$b.r,q,Inf,log.conc.ineq.right) # q to Inf * bound
# 
# 
#     lower.components[1] <- boole(t.pdf$x[t.pdf$index.l],t.pdf$x[t.pdf$index.l:t.pdf$index.r],exp(log(t.pdf$y[t.pdf$index.l:t.pdf$index.r])+log.conc.ineq.left(t.pdf$x[t.pdf$index.l:t.pdf$index.r],q,L,d,psi,nu)),int.to.right = T,one.minus = T) # Middle* bound
#     lower.components[2] <- -ExpInt(t.pdf$a.r,t.pdf$b.r,q,Inf) # Right
#     lower.components[3] <- -ExpBoundInt(t.pdf$a.r,t.pdf$b.r,t.pdf$x[t.pdf$index.r],q,log.conc.ineq.left) # Right * bound
#     lower.components[4] <- -ExpBoundInt(t.pdf$a.l,t.pdf$b.l,-Inf,t.pdf$x[t.pdf$index.l],log.conc.ineq.left) # Left * bound
# 
#   }
# 
#   if( t.pdf$x[t.pdf$index.l]<=q & q<=t.pdf$x[t.pdf$index.r] ){
#     
# 
#     upper.components[1] <- ExpInt(t.pdf$a.l,t.pdf$b.l,-Inf,t.pdf$x[t.pdf$index.l]) # Left
#     upper.components[2] <- boole(q,t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r],int.to.right = F) # Middle Left
#     upper.components[3] <- boole(q,t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r]*exp(log.conc.ineq.right(t.pdf$x[t.pdf$index.l:t.pdf$index.r],q,L,d,psi,nu)),int.to.right = T) # Middle Right * Bound
#     upper.components[4] <- ExpBoundInt(t.pdf$a.r,t.pdf$b.r,t.pdf$x[t.pdf$index.r],Inf,log.conc.ineq.right) # Right * Bound
# 
#     lower.components[1] <- boole(q,t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r],int.to.right = T,one.minus = T) # Middle Right
#     lower.components[2] <- -ExpInt(t.pdf$a.r,t.pdf$b.r,t.pdf$x[t.pdf$index.r],Inf) # Right
#     lower.components[3] <- -boole(q,t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r]*exp(log.conc.ineq.left(t.pdf$x[t.pdf$index.l:t.pdf$index.r],q,L,d,psi,nu)),int.to.right = F) # Middle Left
#     lower.components[4] <- -ExpBoundInt(t.pdf$a.l,t.pdf$b.l,-Inf,t.pdf$x[t.pdf$index.l],log.conc.ineq.left) # Left * bound
# 
#   }
# 
#   c(max(neumaierSum(lower.components), 0),min(neumaierSum(upper.components), 1))
# 
# }


# Boole Quadrature Method
boole<-function(a,x,y,int.to.right=T,one.minus=F){

  # for int.to.right = T , this function integrates the function from a to max(x)
  # for int.to.right = F , this function integrates from min(x) to a
  # to obtain integral over the entire domain, just set a to min(x) and use int.to.right=T

  n <- length(x)

  # Swivel problem around the y axis if we need to go in the other direction
  if(int.to.right==F){a <- - a ; x <- -rev(x) ; y <- rev(y)}

  if(a>x[n]){stop("There are no evaluations of the function on the choosen interval of integration")}

  left.remainder <- right.remainder <- inner.integral <- 0

  # a.index is the index of the smallest element in x that is greater than or equal to a.
  a.index <- findInterval(a,x,left.open = T)+1

  if(x[a.index]!=a){
    # Make trapezoid rule approximation to integrate from observed value to nearest larger grid point
    left.remainder <- (y[a.index-1]+y[a.index])*(x[a.index]-a)/2
  }

  if(a.index!=n){

    num.right.ragged.terms <- (n-a.index)%%4

    if(num.right.ragged.terms==0){
      right.remainder <- 0
      last.index <- n
    }else if(num.right.ragged.terms==1){
      right.remainder <- (y[n]+y[n-1])*(x[n]-x[n-1])/2
      last.index <- n-1
    }else if(num.right.ragged.terms==2){
      right.remainder <- (y[n]+4*y[n-1]+y[n-2])*(x[n]-x[n-2])/6
      last.index <- n-2
    }else if(num.right.ragged.terms==3){
      right.remainder <- 3*(y[n]+3*y[n-1]+3*y[n-2]+y[n-3])*(x[n]-x[n-3])/8
      last.index <- n-3
    }
  }

  left.remainder.scaled <- left.remainder/((x[a.index+1]-x[a.index])*(2/45))
  right.remainder.scaled <- right.remainder/((x[a.index+1]-x[a.index])*(2/45))

  
  w<-c(14,32,12,32)
  
  if(one.minus==F){

    if(n-a.index >= 4){
      vec<-c(w*y[a.index:(last.index-1)],7*y[last.index],-7*y[a.index],left.remainder.scaled,right.remainder.scaled)
      norm.const <- 1/max(c(.Machine$double.eps, vec))
      res <- (x[a.index+1]-x[a.index])*(2/45)*neumaierSum(c(vec*norm.const))/norm.const
    }else{
      vec<-c(left.remainder,right.remainder)
      res <- neumaierSum(vec)
    }

  }else{

    if(n-a.index >= 4){
      vec<-c(45/(2*(x[a.index+1]-x[a.index])),-w*y[a.index:(last.index-1)],-7*y[last.index],7*y[a.index],
             -left.remainder.scaled,-right.remainder.scaled)
      norm.const <- 1/max(c(.Machine$double.eps, vec))
      res <- (x[a.index+1]-x[a.index])*(2/45)*neumaierSum(c(vec/norm.const))*norm.const
    }else{
      vec<-c(1,-left.remainder,-right.remainder)
      res <- neumaierSum(vec)
    }


  }
  
  log(max(0,res))

  # RC: We could consider trying to do the log sum exp trick with boole, restandardizing, with the hope of getting better precision
  # if that doesn't work, maybe need to consider quad precision as an input to boole?
  
  
}


# # Testing Boole
# print(boole(t.pdf$x[t.pdf$index.l],t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r],int.to.right = T),digits=22)
# print(boole(t.pdf$x[60022],t.pdf$x,t.pdf$y,int.to.right = T)+boole(t.pdf$x[60022],t.pdf$x,t.pdf$y,int.to.right = F),digits=22)
#
# print(boole(t.pdf$x[t.pdf$index.l],t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r],int.to.right = T,one.minus = T),digits=22)
# print(boole(t.pdf$x[60022],t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r],int.to.right = T,one.minus = T)+
#         boole(t.pdf$x[60022],t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r],int.to.right = F,one.minus = T),digits=22)

# Seems that sometimes getting sums on either side of 1....not sure if error is coming from
# a) FFT
# b) Quadrature rule
# c) sum / lack of precision

# Maybe we could do integrals from the closer side and subtract 1 instead....be more robust to errors in FFT / simpsons....
# but then just need a way to lower bound the tails---get David's result working?

# Wait, I should probably go back and break up some of these terms and check if they're greater than or equal to 1 at some point
# because we know they shouldn't be and could have some internal regularization checks?



ExpInt<-function(a,b,from,to){

  # This performs the integral $ \int_\from^\to exp(a - b*z) dz $ for any sign of a,b.

  if(from==to){return(-Inf)}
  if(to<from){stop("from must be less than or equal to to")}
  if(b==0){stop("b cannot be zero or else this is just an integral of a constant function")}
  if(all(is.infinite(c(from,to)))){stop("both from and psi cannot be infinite")}
  if(is.infinite(c(from)) & from>0){stop("from cannot be positive infinity")}
  if(is.infinite(c(to)) & to<0){stop("from cannot be negative infinity")}

  if(is.infinite(from)){
    if(b>0){ stop("b cannot be positive when from is negative infinity --> integral does not converge")}
    # b is negative
    return(a - log(abs(b)) - b * to)
  }

  if(is.infinite(to)){
    if(b<0){ stop("b cannot be negative when to is infinity --> integral does not converge")}
    # b is positive
    return(a - log(abs(b)) - b * from)
  }
  
  if(b>0){temp <- -from}else{temp <- to}

  a-log(abs(b))+log(-expm1(-abs(b)*(to-from)))+abs(b)*temp
}
# 
# 
# # Now ExpInt is fixed
-ExpInt(1,-1.2,from=-Inf,to=-2000)/log(10)
-log(integrate(function(z,a,b,k,nu){exp(a-b*z)},-Inf,-2000,a=1,b=-1.2)$value)/log(10)



GaussInt<-function(a,b,from,to,k,nu){
  
  if(from==to){return(-Inf)}
  if(to<from){stop("from must be less than or equal to to")}
  if(b==0){stop("b cannot be zero or else this is just an integral of a constant function")}
  if(all(is.infinite(c(from,to)))){stop("both from and psi cannot be infinite")}
  if(is.infinite(c(from)) & from>0){stop("from cannot be positive infinity")}
  if(is.infinite(c(to)) & to<0){stop("from cannot be negative infinity")}
  
  prefix <- a + b^2*nu/2 - k*b + 0.5*log(2*pi*nu)
  
  if(is.infinite(from)){
    if(b>0){ stop("b cannot be positive when from is negative infinity --> integral does not converge")}
    # b is negative
    return( prefix + pnorm((to-(k-b*nu))/sqrt(nu),log.p = T) )
  }
  
  if(is.infinite(to)){
    if(b<0){ stop("b cannot be negative when to is infinity --> integral does not converge")}
    # b is positive
    return( prefix + pnorm((from-(k-b*nu))/sqrt(nu),log.p = T,lower.tail = F) )
  }
  
  temp<-pnorm((to-(k-b*nu))/sqrt(nu),log.p = T)
  
  prefix + temp + log(-expm1( pnorm((from-(k-b*nu))/sqrt(nu),log.p = T) - pnorm((to-(k-b*nu))/sqrt(nu),log.p = T) )) 
}

# # Looks like GaussInt is correct
# xxx<-seq(-1000,1000,len=100)
# yyy<-xxx
# for(i in 1:length(xxx)) yyy[i] <- -GaussInt(1,1.2,from=xxx[i],to=Inf,0,20)/log(10)
# plot(xxx,yyy)
# 
# -GaussInt(1,-1.2,from=-Inf,to=-1800,0,20)/log(10)
# -log(integrate(function(z,a,b,k,nu){exp(-0.5*(k-z)^2/nu+a-b*z)},-Inf,-800,a=1,b=-1.2,k=0,nu=20)$value)/log(10)



AnalyticBoundInt<-function(a,b,from,to,right=T,q.=q,L.=L,d.=d,psi.=psi,nu.=nu){
  
  if(from==to){return(-Inf)}
  if(to<from){stop("from must be less than or equal to to")}
  if(b==0){stop("b cannot be zero or else this is just an integral of a constant function")}
  if(all(is.infinite(c(from,to)))){stop("both from and psi cannot be infinite")}
  if(is.infinite(c(from)) & from>0){stop("from cannot be positive infinity")}
  if(is.infinite(c(to)) & to<0){stop("psi cannot be negative infinity")}
  
  if(is.infinite(from) & b>0){ stop("b cannot be positive when from is negative infinity --> integral does not converge")}
  
  if(is.infinite(to) & b<0){stop("b cannot be negative when to is infinity --> integral does not converge")}
  
  # This must integrate exp(a-bz) * bound where bound may be in one of three stages (flat, gauss, or exp)
  
  k.<-q.-psi.
  
  rho <- as.integer(right)
  
  if(right==T){
    const.interval <- c(-Inf, k.)
    gauss.interval <- c(k., nu./(4*d.)+k.)
    exp.interval <- c(nu./(4*d.)+k., Inf)
  }else{
    const.interval <- c(k., Inf)
    gauss.interval <- c(k.-nu./(4*d.), k.)
    exp.interval <- c(-Inf, k.-nu./(4*d.))
  }
  
  const.component <- gauss.component <- exp.component <- -Inf
  
  ##### Constant Component #####
  
  # Interval inside region of integration
  if(from <= const.interval[1] & const.interval[2]<=to){
    const.component <- ExpInt(a,b,const.interval[1],const.interval[2])
  }
  
  # Interval overlaps left part of region of integration
  if(const.interval[1] < from & from < const.interval[2] & const.interval[2] <= to){
    const.component <- ExpInt(a,b,from,const.interval[2])
  }
  
  # Interval covers entire region of integration
  if(const.interval[1] < from & to < const.interval[2] ){
    const.component <- ExpInt(a,b,from,to)
  }
  
  # Interval overlaps right part of region of integration
  if(from <= const.interval[1] & const.interval[1] < to & to < const.interval[2]){
    const.component <- ExpInt(a,b,const.interval[1],to)
  }
  
  ##### Gauss Component #####
  
  # Interval inside region of integration
  if(from <= gauss.interval[1] & gauss.interval[2]<=to){
    gauss.component <- GaussInt(a,b,gauss.interval[1],gauss.interval[2],k.,nu.)
  }
  
  # Interval overlaps left part of region of integration
  if(gauss.interval[1] < from & from < gauss.interval[2] & gauss.interval[2] <= to){
    gauss.component <- GaussInt(a,b,from,gauss.interval[2],k.,nu.)
  }
  
  # Interval covers entire region of integration
  if(gauss.interval[1] < from & to < gauss.interval[2] ){
    gauss.component <- GaussInt(a,b,from,to,k.,nu.)
  }
  
  # Interval overlaps right part of region of integration
  if(from <= gauss.interval[1] & gauss.interval[1] < to & to < gauss.interval[2]){
    gauss.component <- GaussInt(a,b,gauss.interval[1],to,k.,nu.)
  }
  
  
  ##### Exp Component #####
  
  # Interval inside region of integration
  if(from <= exp.interval[1] & exp.interval[2]<=to){
    exp.component <- ExpInt(a+0.5*nu./((4*d.)^2)+rho*k./(4*d.),b+rho/(4*d.),exp.interval[1],exp.interval[2])
  }
  
  # Interval overlaps left part of region of integration
  if(exp.interval[1] < from & from < exp.interval[2] & exp.interval[2] <= to){
    exp.component <- ExpInt(a+0.5*nu./((4*d.)^2)+rho*k./(4*d.),b+rho/(4*d.),from,exp.interval[2])
  }
  
  # Interval covers entire region of integration
  if(exp.interval[1] < from & to < exp.interval[2] ){
    exp.component <- ExpInt(a+0.5*nu./((4*d.)^2)+rho*k./(4*d.),b+rho/(4*d.),from,to)
  }
  
  # Interval overlaps right part of region of integration
  if(from <= exp.interval[1] & exp.interval[1] < to & to < exp.interval[2]){
    exp.component <- ExpInt(a+0.5*nu./((4*d.)^2)+rho*k./(4*d.),b+rho/(4*d.),exp.interval[1],to)
  }
  
  # Components are here in log space
  temp <- c(const.component,gauss.component,exp.component)
  norm.const <- max(temp)
  # Return result in the probability space
  neumaierSum(c(log(neumaierSum(exp(temp-norm.const))),norm.const))
}


# xx<-seq(0,10000,len=1000)
# yy<-rep(0,length(xx))
# for(i in 1:length(xx)) yy[i]<--log10(ExpBoundInt(t.pdf$a.l,t.pdf$b.l,-100,xx[i],log.conc.ineq = log.conc.ineq.l))
# plot(xx,yy,col="red",type="l")
# 
# xx<-seq(-1540,-1480,len=1000)
# yy<-rep(0,length(xx))
# for(i in 1:length(xx)) yy[i]<--AnalyticBoundInt(t.pdf$a.l,t.pdf$b.l,-1550,xx[i],right = F)/log(10)
# plot(xx,yy,type="l",col="blue")
# abline(v=-1501)
# abline(v=-1526)
# 
# 
# print(-AnalyticBoundInt(t.pdf$a.r,t.pdf$b.r,2410.521,Inf,right=T,q.=q,L.=L,d.=d,psi.=psi,nu.=nu)/log(10),digits=22)
# 
# xx<-seq(2410.521,4000,len=1000)
# yy<-rep(0,length(xx))
# for(i in 1:length(xx)) yy[i]<--AnalyticBoundInt(t.pdf$a.r,t.pdf$b.r,xx[i],Inf,right = T)/log(10)
# plot(xx,yy,type="l")
# 
# xx<-seq(-10000,10000,len=1000)
# yy<-rep(0,length(xx))
# for(i in 1:length(xx)) yy[i]<--AnalyticBoundInt(t.pdf$a.l,t.pdf$b.l,-Inf,xx[i],right = F)/log(10)
# plot(xx,yy,type="l")
# 
# 



# 
# ExpBoundInt<-function(a,b,from,to,log.conc.ineq,q.=q,L.=L,d.=d,psi.=psi,nu.=nu){
# 
#   # This performs the integral $ \int_\from^\to exp(a - b*z + log.conc.ineq(z,q,L,d,psi,nu)) dz $ for any sign of a,b.
#   
#   if(from==to){return(0)}
#   if(to<from){stop("from must be less than or equal to to")}
#   if(b==0){stop("b cannot be zero or else this is just an integral of a constant function")}
#   if(all(is.infinite(c(from,to)))){stop("both from and psi cannot be infinite")}
#   if(is.infinite(c(from)) & from>0){stop("from cannot be positive infinity")}
#   if(is.infinite(c(to)) & to<0){stop("psi cannot be negative infinity")}
# 
#   if(is.infinite(from) & b>0){ stop("b cannot be positive when from is negative infinity --> integral does not converge")}
# 
#   if(is.infinite(to) & b<0){stop("b cannot be negative when to is infinity --> integral does not converge")}
# 
#   integrate(f=function(z){
#     exp(a - b*z + log.conc.ineq(z,q.,L.,d.,psi.,nu.))
#   }, lower = from, upper = to, rel.tol = .Machine$double.eps)$value
# 
#   
# }
# 

# This is h_k(t) meant to be used in obtaining upper bounds
log.conc.ineq.r <-function(z,q,L,d,psi,nu){
  k<-q-psi
  ifelse(z <= k ,0,ifelse(z <= nu/(4*d)+k,-0.5*((k-z)^2)/nu, 0.5*nu/((4*d)^2) + (k-z)/(4*d)))
}

# This is h_-k(-t) meant to be used in obtaining lower bounds
log.conc.ineq.l<-function(z,q,L,d,psi,nu){
  log.conc.ineq.r(-z,-q,L,d,-psi,nu)
}

# Concentration Inequality bounds look right
# xx<-seq(q-psi-30,q-psi+nu/(4*d)+30,len=1000)
# plot(xx,-log.conc.ineq.r(xx,q,L,d,psi,nu)/log(10),type="l")
# abline(v=q-psi)
# abline(v=q-psi+nu/(4*d))

# 
# plot(xx,-log.conc.ineq.l(xx,q,L,d,psi,nu)/log(10),type="l")
# abline(v=(q-psi))
# abline(v= c(q-psi)-nu/(4*d))


# # Analytic Exponential Integrals Look Good
# xx<-seq(-2000,2000,len=1000)
# yy<-xx
# for(i in 1:length(xx)) yy[i]<-log10(ExpInt(t.pdf$a.r,t.pdf$b.r,xx[i],Inf))
# plot(xx,yy)




######
# Example Code Here
#######


# Gauss Quadrature Exponential * Bound Integral looks....bad
# 
# require(QForm)
# require(PreciseSums)
# 
# evals<- c(-10:-1,1:10)
# ncps<-1:20
# 
# q<--1500
# L<-1
# d<-1
# psi<-1
# nu<-1000
# 
# t.cdf<-Tcdf_new(evals,ncps)
# 
# plot(t.pdf$x,t.pdf$y,type="l")
# 
# 
# t.pdf<-Tpdf(evals,ncps)
# 
# QFIntegrate2(-4000,t.pdf,log.conc.ineq.l,log.conc.ineq.r,L = L,d=d,psi=psi,nu=nu)

#
# 
# xx<-seq(-2000,2000,len=1000)
# yy<-matrix(nrow=length(xx),ncol=2)
# system.time(for(i in 1:length(xx)){
#   yy[i,]<-QFIntegrate2(xx[i],t.pdf,log.conc.ineq.l,log.conc.ineq.r,L = L,d=d,psi=psi,nu=nu)
# })
# 
# 
# plot(xx,yy[,1],type="l",col="red")
# points(xx,yy[,2],type="l",col="blue")
# 
# 
# plot(xx,-log1p(-yy[,1])/log(10),type="l",col="red",ylim=c(0,5))
# lines(xx,-log1p(-yy[,2])/log(10),col="blue")
# 
# plot(xx,-log(yy[,1])/log(10),type="l",col="red",ylim=c(0,5))
# lines(xx,-log(yy[,2])/log(10),col="blue")
# 
# 
# 
# xxx<-seq(-4000,4000,len=100)
# yyy<-matrix(nrow=length(xxx),ncol=3)
# system.time(for(i in 1:length(xxx)){
#   yyy[i,]<-as.numeric(as.vector(QFIntBounds2(obs = xxx[i],evals = evals,ncps=ncps,E_R=psi,nu = nu,resid.op.norm.bd = d)))
#     })
# 
# plot(xxx,yyy[,1],type="l",col="red")
# points(xxx,yyy[,2],type="l",col="blue")
# 
# 
# plot(xxx,-log1p(-yyy[,1])/log(10),type="l",col="red",ylim=c(0,18))
# lines(xxx,-log1p(-yyy[,2])/log(10),col="blue")
# 
# plot(xxx,-log(yyy[,1])/log(10),type="l",col="red",ylim=c(0,18))
# lines(xxx,-log(yyy[,2])/log(10),col="blue")
# 
# # Confirm that it's not a bug in how we combine results or a numerical stability problem combining results
# # Confirm it's not a prob in Analytic integral logic
# # See if there's a problem in a particular term
# # If not any of these then we have a problem I think in boole:
# # check for a problem in boole...where is it hitting 1 or 0 ?
# # Could rewrite boole function to do integral of exp(g) where we standardize g to keep exp in a better interval
# # Could go with quad precision
# 
# 
# 
# 
# 
# 
# 
# 
# xx<-seq(-500,4000,len=1000)
# yy<-rep(0,length(yy))
# for(i in 1:length(xx)) yy[i]<--log10(ExpBoundInt(t.pdf$a.r,t.pdf$b.r,xx[i],Inf,log.conc.ineq = log.conc.ineq.l))
# plot(xx,yy)
# 
# 
# xx<-seq(-500,4000,len=1000)
# yy<-rep(0,length(yy))
# for(i in 1:length(xx)) yy[i]<--log10(AnalyticBoundInt(t.pdf$a.r,t.pdf$b.r,xx[i],Inf,right = F))
# plot(xx,yy)
# 
# xx<-seq(-1000,4000,len=1000)
# yy<-rep(0,length(yy))
# for(i in 1:length(xx)) yy[i]<--log10(ExpBoundInt(t.pdf$a.l,t.pdf$b.l,-Inf,xx[i],log.conc.ineq = log.conc.ineq.r))
# plot(xx,yy)
# 
# ######
# ### Looks like Analytic Bound does the wrong thing in this case...need to check this out
# #######
# xx<-seq(-1000,4000,len=1000)
# yy<-rep(0,length(yy))
# for(i in 1:length(xx)) yy[i]<- AnalyticBoundInt(t.pdf$a.l,t.pdf$b.l,-Inf,xx[i])
# plot(xx,yy)
# 
# xx<-seq(-2000,2000,len=1000)
# yy<-rep(0,length(yy))
# for(i in 1:length(xx)) yy[i]<--log10(ExpBoundInt(t.pdf$a.l,t.pdf$b.l,-Inf,xx[i],log.conc.ineq = log.conc.ineq.l))
# plot(xx,yy,type="l")
# 
# ############
# # This one is also different!!!
# ##########
# 
# xx<-seq(-2000,2000,len=1000)
# yy<-rep(0,length(yy))
# for(i in 1:length(xx)) yy[i]<--log10(AnalyticBoundInt(t.pdf$a.l,t.pdf$b.l,-Inf,xx[i],right = F))
# plot(xx,yy,type="l")
# 
# 
# 
# #######
# 
# 
# 
# 
# QFIntegrate2(q,t.pdf,log.conc.ineq.l,log.conc.ineq.r,L = L,d=d,psi=psi,nu=nu)
# 






