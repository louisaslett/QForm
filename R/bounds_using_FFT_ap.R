
#lambda<-c(1e10,-1000:1000,-1e12)
lambda<-mpfr(c(200,-300,100:110),acc)
ncps<-mpfr(1:length(lambda),acc)


a<- mpfr(-1000000,acc)
b<- mpfr(1000000,acc)
n.forR<-2^16
n<-mpfr(n.forR,acc)


c.div <- function(x, y) {
  d<-rowSums(y^2)
  matrix(c(
    (x[,1]*y[,1] + x[,2]*y[,2])/d,
    (x[,2]*y[,1] - x[,1]*y[,2])/d
  ), nrow = nrow(x), ncol = 2)
}

log.phi.ap<-function(t,lambda,ncps){
  x<-cbind(mpfr(0,acc),ncps*t*lambda)
  y<-cbind(mpfr(1,acc),-2*t*lambda)
  matrix(colSums(c.div(x,y) - 0.5*c.log(y)),ncol=2)
}

rho.ap<-mpfr(matrix(0,nrow=n.forR,ncol=2),acc)
require(foreach)
require(doMC)
registerDoMC(cores = 8)

rho.ap<-foreach(k = 1:n.forR,.combine=rbind,.multicombine = TRUE, .maxcombine = 10000) %dopar% c.exp(log.phi.ap((Const("pi",acc)*n/(b-a))*(2*((k-1)/n)-1), lambda, ncps) - matrix(c(mpfr(0,acc),Const("pi",acc)*a*(n/(b-a))*(2*((k-1)/n)-1)),ncol=2))

rho<-complex(real=as.double(rho.ap[,1]),imaginary=as.double(rho.ap[,2]))

complex(imaginary = pi*a*(n/(b-a))*(2*((k-1)/n)-1))

xx<-seq(a,b-(b-a)/n.forR,length.out = n.forR)

plot(xx,-log10(rho.ap[,1]),type="l")
points(xx,-log10(Re(rho)),type="l", col="red")

ans<-Re(fft(rho)*(-1)^(0:(n.forR-1))/as.numeric(b-a))

ans2<-fft.ap2(rho.ap)[,1]*(-1)^(0:(n-1))/(b-a)

plot(xx,-log10(ans),type="l")
lines(xx,-log10(ans2),col="red")



ans2<-fft.ap(rho)[,1]*(-1)^(0:(n-1))/(b-a)

