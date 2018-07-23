#####
# Evaluating QForm Bounds using FFT
#####

# Open Questions:
# 1) Work to get better resolution in the tail....maybe use safer methods?
# 2) Make sure numerically safe
# 3) Implement automatic determination of a, b and n based on quantiles desired and bounds on the modulus

# KALIS = kit for acelerated li and stephens
# Alis = acelerated li and Stephens

#lambda<-c(1e10,-1000:1000,-1e12)
lambda<-c(200,-300,100:110)
ncps<-1:length(lambda)

# If you don't make a and b large enough, you get nonsense


a<- -1000000
b<- 1000000
n<-2^16
# There appears to be a limit...the relative error with respect to davies doesn't improve with increasing n.
log.phi<-function(t,lambda,ncps){
  sum(complex(imaginary=ncps*t*lambda)/complex(real=1, imaginary=-2*t*lambda) - 0.5*log(complex(real=1,imaginary=-2*t*lambda)))
}


rho<-1:n
for(k in 1:n)  rho[k]<-exp(log.phi( (pi*n/(b-a))*(2*((k-1)/n)-1), lambda, ncps) - complex(imaginary = pi*a*(n/(b-a))*(2*((k-1)/n)-1)) )

xx<-seq(a,b-(b-a)/n,length.out = n)

plot(xx,-log10(rho),type="l")

ans<-Re(fft(rho)*(-1)^(0:(n-1))/(b-a))

plot(xx,-log10(ans),type="l")
lines(xx,-log10(ans2),col="red")


dF<-ifelse(ans<10^(-16),10^(-16),ans)

plot(xx, dF,type="l",xlab="Q",ylab="density")
plot(xx, -log10(dF),type="l",xlab="Q",ylab="density")

# Here we verify that we're getting the same answer as CompQuadForm

require(CompQuadForm)
true<-xx
for( i in 1:length(xx)){
true[i]<-1-davies(q =xx[i],lambda = lambda,delta = ncps)$Qq
}
plot(true,type="l",col="blue")
lines(cumsum(ans)*((b-a)/n),type="l",col="red")
plot((cumsum(ans)*((b-a)/n)-true)/true,type="l")
plot(true-(cumsum(ans)*((b-a)/n))/(cumsum(ans)*((b-a)/n)),type="l")


# Now lets integrate this with respect to our bound

# Boole Quadrature Method
boole<-function(x,a,b){
  n<-length(x)
  if((n-1)%%4!=0) stop("x must be of length 4N+1 for some natural number N")
  w<-c(14,32,12,32)
  ((b-a)/(n-1))*(2/45)*(sum(w*x[-n])+7*x[n]-7*x[1])
}


# User provides
N<-8
a<-0
b<-pi

# Calculate end point
b.prime<-b+3*(b-a)/((2^N)-4)

# Evaluate function on grid with 2^N points and throw last 3 away
xx<-seq(a,b.prime,length.out = 2^N)
yy<-dcauchy(xx,scale = 10)[-c(2^N-0:2)]

# Check the error
print(boole(yy,0,pi)-(pcauchy(b,scale=10)-pcauchy(a,scale=10)),digits=20)
