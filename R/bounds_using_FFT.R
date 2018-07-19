#####
# Evaluating QForm Bounds using FFT
#####

#lambda<-c(1e10,-1000:1000,-1e12)
lambda<-c(-5:5,-10,20,30)
ncps<-1:length(lambda)



a<--10000
b<-10000
n<-200000

log.phi<-function(t,lambda,ncps){
  sum(complex(imaginary=ncps*t*lambda)/complex(real=1, imaginary=-2*t*lambda) - 0.5*log(complex(real=1,imaginary=-2*t*lambda)))
}

rho<-1:n
for(k in 1:n){
  rho[k]<-exp(log.phi( (pi*n/(b-a))*(2*((k-1)/n)-1), lambda, ncps) - complex(imaginary = pi*a*(n/(b-a))*(2*((k-1)/n)-1)) )
  print(k)
}

ans<-Re(fft(rho)*(-1)^(0:(n-1))/(b-a))
xx<-seq(a,b-(b-a)/n,length.out = n)
plot(xx, ans,type="l",xlab="Q",ylab="density")

true<-xx
for( i in 1:length(xx)){
true[i]<-1-davies(q =xx[i],lambda = lambda,delta = ncps)$Qq
}

plot(true,type="l",col="blue")
lines(cumsum(ans)*((b-a)/n),type="l",col="red")
plot((cumsum(ans)*((b-a)/n)-true),type="l")



require(CompQuadForm)

davies(q =  ,lambda = lambda,delta = ncps)$Qq

require(parallel)

sum.phi<-function(k,n,l,lambda,ncps){
  #cl<-makeCluster(detectCores() - 1, type="FORK")
  #clusterExport(cl,c("t","lambda","ncps","log.phi"))
  #ans<-parSapply(cl, X = -l:l, function(x){exp(log.phi((n/(b-a))*(2*pi*((k-1)/n)+x),lambda,ncps)-complex(imaginary = a*(n/(b-a))*(2*pi*((k-1)/n)+x)))})
  ans<-sapply(X = -l:l, function(x){exp(log.phi((n/(b-a))*(2*pi*((k-1)/n)+x),lambda,ncps)-complex(imaginary = a*(n/(b-a))*(2*pi*((k-1)/n)+x)))})
  #stopCluster(cl)
  sum(ans)
}

y<-1:n
for(k in 1:n){
y[k]<-sum.phi(k,n,1000,lambda,ncps)
print(k)
}

x<-fft(y,inverse=FALSE)/(b-a)
plot(Re(x))
#


#require(PreciseSums)

# logphi<-function(t){
#   v<-complex(imaginary=ncps*t*lambda)/complex(real=1, imaginary=-2*t*lambda) - 0.5*log(complex(real=1,imaginary=-2*t*lambda))
#   complex(real=fsum(Re(v)),imaginary=fsum(Im(v)))
# }
#
#
# sum.phi<-function(t,l){
#   v<--l:l
#   for(i in -l:l){
#     v[i] <- logphi(2*pi*t+i)
#   }
#   mrv<-max(Re(v))
#   miv<-max(Im(v))
#   w<-exp(v-complex(real=mrv,imaginary=miv))
#   exp(complex(real=mrv,imaginary=miv)+log(complex(real=fsum(Re(w)),fsum(Im(w)))))
# }

sum.phi(2,10000)




xx<-seq(0,10,by=0.001)
yy<-xx
for(i in 1:length(xx)){
  yy[i]<-phi(xx[i])
}

zz<-fft(yy,inverse = TRUE)/length(yy)
plot(xx,zz,type="l")

zz2<-fft(zz,inverse=FALSE)
plot(xx,zz2,type="l")


plot(xx,Re(yy),type="l")
#


complex(imaginary=ncps*t*lambda)
complex(real=c(1,2,3),imaginary=1)
