# 
# 
# ## Trial vignette
# 
# lambda <- rep(c(100,50,20,10,50,-90,-80,70,-60),2)
# cdf <- QForm(lambda=lambda)
# 
# xx <- seq(-1e4,1e4,len=250)
# plot(xx,cdf(xx),type="l")
# plot(xx,cdf(xx,density = T),type="l")
# 
# cdf.bounds<-QBounds(cdf = cdf,f="identity", bound.args = list("Er"=8,"nu"=1e4,"L"=10))
# # lambda.remain
# # lambda.remain.sq
# # bounds args is an option for advanced users..
# 
# #QForm 
# #QBounds
# 
# 
# 
# yy<-cdf.bounds(xx)
# 
# plot(xx,cdf(xx),type="l")
# lines(xx,yy[1,],col="blue",lty=2)
# lines(xx,yy[2,],col="red",lty=2)
# 
# 
# plot(xx,-cdf(xx,log.p = T)/log(10),type="l")
# lines(xx,-log10(yy[1,]),col="blue",lty=2)
# lines(xx,-log10(yy[2,]),col="red",lty=2)
# 
# 
# # Boole Quadrature Method
# require(PreciseSums)
# boole<-function(a,x,y,int.to.right=T,one.minus=F){
#   
#   # for int.to.right = T , this function integrates the function from a to max(x)
#   # for int.to.right = F , this function integrates from min(x) to a
#   # to obtain integral over the entire domain, just set a to min(x) and use int.to.right=T
#   
#   n <- length(x)
#   
#   # Swivel problem around the y axis if we need to go in the other direction
#   if(int.to.right==F){a <- - a ; x <- -rev(x) ; y <- rev(y)}
#   
#   if(a>x[n]){stop("There are no evaluations of the function on the choosen interval of integration")}
#   
#   left.remainder <- right.remainder <- inner.integral <- 0
#   
#   # a.index is the index of the smallest element in x that is greater than or equal to a.
#   a.index <- findInterval(a,x,left.open = T)+1
#   
#   if(x[a.index]!=a){
#     # Make trapezoid rule approximation to integrate from observed value to nearest larger grid point
#     left.remainder <- (y[a.index-1]+y[a.index])*(x[a.index]-a)/2
#   }
#   
#   if(a.index!=n){
#     
#     num.right.ragged.terms <- (n-a.index)%%4
#     
#     if(num.right.ragged.terms==0){
#       right.remainder <- 0
#       last.index <- n
#     }else if(num.right.ragged.terms==1){
#       right.remainder <- (y[n]+y[n-1])*(x[n]-x[n-1])/2
#       last.index <- n-1
#     }else if(num.right.ragged.terms==2){
#       right.remainder <- (y[n]+4*y[n-1]+y[n-2])*(x[n]-x[n-2])/6
#       last.index <- n-2
#     }else if(num.right.ragged.terms==3){
#       right.remainder <- 3*(y[n]+3*y[n-1]+3*y[n-2]+y[n-3])*(x[n]-x[n-3])/8
#       last.index <- n-3
#     }
#   }
#   
#   left.remainder.scaled <- left.remainder/((x[a.index+1]-x[a.index])*(2/45))
#   right.remainder.scaled <- right.remainder/((x[a.index+1]-x[a.index])*(2/45))
#   
#   
#   w<-c(14,32,12,32)
#   
#   if(one.minus==F){
#     
#     if(n-a.index >= 4){
#       vec<-c(w*y[a.index:(last.index-1)],7*y[last.index],-7*y[a.index],left.remainder.scaled,right.remainder.scaled)
#       norm.const <- 1/max(c(.Machine$double.eps, vec))
#       res <- (x[a.index+1]-x[a.index])*(2/45)*neumaierSum(c(vec*norm.const))/norm.const
#     }else{
#       vec<-c(left.remainder,right.remainder)
#       res <- neumaierSum(vec)
#     }
#     
#   }else{
#     
#     if(n-a.index >= 4){
#       vec<-c(45/(2*(x[a.index+1]-x[a.index])),-w*y[a.index:(last.index-1)],-7*y[last.index],7*y[a.index],
#              -left.remainder.scaled,-right.remainder.scaled)
#       norm.const <- 1/max(c(.Machine$double.eps, vec))
#       res <- (x[a.index+1]-x[a.index])*(2/45)*neumaierSum(c(vec/norm.const))*norm.const
#     }else{
#       vec<-c(1,-left.remainder,-right.remainder)
#       res <- neumaierSum(vec)
#       #res <- sum(sort(vec))
#     }
#     
#     
#   }
#   
#   res
# }
# # # # Testing Boole
# # testfunc<-function(x){sin(x*10)/(10*x)}
# # xx<-seq(1,3,len=10000)
# # boole(2,xx,testfunc(xx))
# # integrate(f = testfunc,lower = 2,upper = 3,rel.tol = 1e-14)
# #
# #
# 
# QFIntegrate3<-function(ql,qu,t.cdf,conc.ineqs,lower.tail = T, log.p = F){ # L= 1/(4*max(abs(lambda)))
#   
#   # Initialize
#   upper.components <- lower.components <- rep(0,5)
#   F.cdf <- wrap.QFcdf(t.cdf)
#   
#   ### Compute Upper Bound
#   
#   upper.components[1] <- (1-conc.ineqs$l(-conc.ineqs$c2))*F.cdf(qu)
#   
#   if(qu < t.cdf$x[1]){
#     
#     if(t.cdf$type=="ind"){
#       upper.components[2] <- conc.ineqs$int.h2.expx(qu,t.cdf$x[1],qu,t.cdf$a.l,t.cdf$b.l)
#       upper.components[4] <- conc.ineqs$int.h2.const(t.cdf$x[t.cdf$n],Inf,qu)
#       upper.components[5] <- -conc.ineqs$int.h2.expx(t.cdf$x[t.cdf$n],Inf,qu,t.cdf$a.r,t.cdf$b.r)
#     }
#     if(t.cdf$type=="psd"){
#       upper.components[2] <- conc.ineqs$int.h2.explogx(max(qu,0),t.cdf$x[1],qu,t.cdf$a.l,t.cdf$b.l)
#       upper.components[4] <- conc.ineqs$int.h2.const(t.cdf$x[t.cdf$n],Inf,qu)
#       upper.components[5] <- -conc.ineqs$int.h2.expx(t.cdf$x[t.cdf$n],Inf,qu,t.cdf$a.r,t.cdf$b.r) 
#     }
#     if(t.cdf$type=="nsd"){
#       upper.components[2] <- conc.ineqs$int.h2.expx(qu,t.cdf$x[1],qu,t.cdf$a.l,t.cdf$b.l)
#       upper.components[4] <- conc.ineqs$int.h2.const(t.cdf$x[t.cdf$n],Inf,qu) # This integrates out to infinity because we need to account for the constant plateau of the NSD density above 0
#       upper.components[5] <- -conc.ineqs$int.h2.explognegx(t.cdf$x[t.cdf$n],0,qu,t.cdf$a.r,t.cdf$b.r)
#     }
#     
#     #upper.components[3] <- boole(t.cdf$x[1],t.cdf$x, conc.ineqs$h2(qu,t.cdf$x)*t.cdf$y) #Boole integral from t.cdf$x[1] to t.cdf$x[n]
#     upper.components[3] <- GaussQuadCDF(t.cdf$x[1], t.cdf$x[t.cdf$n], F, F.cdf, ql, qu, conc.ineqs) #Gauss Quad integral from t.cdf$x[1] to t.cdf$x[n]
#     
#   }
#   
#   if((qu >= t.cdf$x[1]) & (qu < t.cdf$x[t.cdf$n])){
#     
#     if(t.cdf$type=="ind"){
#       upper.components[4] <- conc.ineqs$int.h2.const(t.cdf$x[t.cdf$n],Inf,qu)
#       upper.components[5] <- -conc.ineqs$int.h2.expx(t.cdf$x[t.cdf$n],Inf,qu,t.cdf$a.r,t.cdf$b.r)
#     }
#     if(t.cdf$type=="psd"){
#       upper.components[4] <- conc.ineqs$int.h2.const(t.cdf$x[t.cdf$n],Inf,qu)
#       upper.components[5] <- -conc.ineqs$int.h2.expx(t.cdf$x[t.cdf$n],Inf,qu,t.cdf$a.r,t.cdf$b.r) 
#     }
#     if(t.cdf$type=="nsd"){
#       upper.components[4] <- conc.ineqs$int.h2.const(t.cdf$x[t.cdf$n],Inf,qu) # This integrates out to infinity because we need to account for the constant plateau of the NSD density above 0
#       upper.components[5] <- -conc.ineqs$int.h2.explognegx(t.cdf$x[t.cdf$n],0,qu,t.cdf$a.r,t.cdf$b.r)
#     }
#     
#     #upper.components[3] <- boole(qu,t.cdf$x, conc.ineqs$h2(qu,t.cdf$x)*t.cdf$y) #Boole integral from qu to t.cdf$x[n]
#     upper.components[3] <- GaussQuadCDF(qu, t.cdf$x[t.cdf$n], F, F.cdf, ql, qu, conc.ineqs) #Gauss Quad integral from qu to t.cdf$x[n]
#     
#   }
#   
#   if(qu >= t.cdf$x[t.cdf$n]){
#     
#     if(t.cdf$type=="ind"){
#       upper.components[4] <- conc.ineqs$int.h2.const(qu,Inf,qu)
#       upper.components[5] <- -conc.ineqs$int.h2.expx(qu,Inf,qu,t.cdf$a.r,t.cdf$b.r)
#     }
#     if(t.cdf$type=="psd"){
#       upper.components[4] <- conc.ineqs$int.h2.const(qu,Inf,qu)
#       upper.components[5] <- -conc.ineqs$int.h2.expx(qu,Inf,qu,t.cdf$a.r,t.cdf$b.r) 
#     }
#     if(t.cdf$type=="nsd"){
#       upper.components[4] <- conc.ineqs$int.h2.const(qu,Inf,qu) # This integrates out to infinity because we need to account for the constant plateau of the NSD density above 0
#       upper.components[5] <- -conc.ineqs$int.h2.explognegx(qu,0,qu,t.cdf$a.r,t.cdf$b.r)
#     }
#     
#   }
#   
#   if(lower.tail){
#     upper.bound <- sum(upper.components)
#   }else{
#     upper.bound <- sum(1,-upper.components[1:4])-upper.components[5]
#   }
#   
#   
#   ### Compute Lower Bound
#   
#   lower.components[1] <- (1-conc.ineqs$u(-conc.ineqs$c1))*F.cdf(ql)
#   
#   if(ql > t.cdf$x[t.cdf$n]){
#     
#     if(t.cdf$type=="ind"){
#       lower.components[2] <- conc.ineqs$int.h1.expx(-Inf,t.cdf$x[1],ql,t.cdf$a.l,t.cdf$b.l)
#       lower.components[4] <- conc.ineqs$int.h1.const(t.cdf$x[t.cdf$n],ql,ql)
#       lower.components[5] <- -conc.ineqs$int.h1.expx(t.cdf$x[t.cdf$n],ql,ql,t.cdf$a.r,t.cdf$b.r)
#     }
#     if(t.cdf$type=="psd"){
#       lower.components[2] <- conc.ineqs$int.h1.explogx(0,t.cdf$x[1],ql,t.cdf$a.l,t.cdf$b.l)
#       lower.components[4] <- conc.ineqs$int.h1.const(t.cdf$x[t.cdf$n],ql,ql)
#       lower.components[5] <- -conc.ineqs$int.h1.expx(t.cdf$x[t.cdf$n],ql,ql,t.cdf$a.r,t.cdf$b.r) 
#     }
#     if(t.cdf$type=="nsd"){
#       lower.components[2] <- conc.ineqs$int.h1.expx(-Inf,t.cdf$x[1],ql,t.cdf$a.l,t.cdf$b.l)
#       lower.components[4] <- conc.ineqs$int.h1.const(t.cdf$x[t.cdf$n],ql,ql) # This integrates out to ql instead of min(ql,0) because we need to account for the constant plateau of the NSD density above 0
#       lower.components[5] <- -conc.ineqs$int.h1.explognegx(t.cdf$x[t.cdf$n],min(0,ql),ql,t.cdf$a.r,t.cdf$b.r)
#     }
#     
#     #lower.components[3] <- boole(t.cdf$x[1],t.cdf$x, conc.ineqs$h1(ql,t.cdf$x)*t.cdf$y) #Boole integral from t.cdf$x[1] to t.cdf$x[n]
#     lower.components[3] <- GaussQuadCDF(t.cdf$x[1], t.cdf$x[t.cdf$n], T, F.cdf, ql, qu, conc.ineqs) #Gauss Quad integral from t.cdf$x[1] to t.cdf$x[n]
#     
#   }
#   
#   if((ql > t.cdf$x[1]) & (ql <= t.cdf$x[t.cdf$n])){
#     
#     if(t.cdf$type=="ind"){
#       lower.components[2] <- conc.ineqs$int.h1.expx(-Inf,t.cdf$x[1],ql,t.cdf$a.l,t.cdf$b.l)
#     }
#     if(t.cdf$type=="psd"){
#       lower.components[2] <- conc.ineqs$int.h1.explogx(0,t.cdf$x[1],ql,t.cdf$a.l,t.cdf$b.l)
#     }
#     if(t.cdf$type=="nsd"){
#       lower.components[2] <- conc.ineqs$int.h1.expx(-Inf,t.cdf$x[1],ql,t.cdf$a.l,t.cdf$b.l)
#     }
#     
#     #lower.components[3] <- boole(ql,t.cdf$x, conc.ineqs$h1(ql,t.cdf$x)*t.cdf$y,int.to.right = F) #Boole integral from t.cdf$x[1] to t.cdf$x[n]
#     lower.components[3] <- GaussQuadCDF(t.cdf$x[1], ql, T, F.cdf, ql, qu, conc.ineqs) #Gauss Quad integral from t.cdf$x[1] to ql
#   }
#   
#   if(ql <= t.cdf$x[1]){
#     
#     if(t.cdf$type=="ind"){
#       lower.components[2] <- conc.ineqs$int.h1.expx(-Inf,ql,ql,t.cdf$a.l,t.cdf$b.l)
#     }
#     if(t.cdf$type=="psd"){
#       lower.components[2] <- conc.ineqs$int.h1.explogx(0,ql,ql,t.cdf$a.l,t.cdf$b.l)
#     }
#     if(t.cdf$type=="nsd"){
#       lower.components[2] <- conc.ineqs$int.h1.expx(-Inf,ql,ql,t.cdf$a.l,t.cdf$b.l)
#     }
#     
#   }
#   
#   if(lower.tail){
#     lower.bound <- sum(lower.components)
#   }else{
#     lower.bound <- sum(1,-lower.components[1:4])-lower.components[5]
#   }
#   
#   print(c(qu,upper.components))
#   print(c(ql,lower.components))
#   
#   # Make sure the results fall within [0,1]
#   
#   if(lower.bound < 0){lower.bound <- 0}
#   if(lower.bound > 1){lower.bound <- 1}
#   
#   if(upper.bound < 0){upper.bound <- 0}
#   if(upper.bound > 1){upper.bound <- 1}
#   
#   
#   if(log.p){
#     ans <- log(c(lower.bound,upper.bound))
#   }else{
#     ans <- c(lower.bound,upper.bound)
#   }
#   ans
# }
# 
# 
# 
# 
# 
# nu <- 1000
# L <- 10
# Er <- 20
# #conc.ineqs<-WrapConcIneq.identity(Er,Er,nu,L)
# 
# 
# 
# 
# 
# 
# evals <- c(1000,50,20,10,50,-90,-80,70,-60)
# 
# t.cdf <- calc.QFcdf(evals)
# conc.ineqs<-WrapConcIneq.identity(Er,Er,nu,L)
# 
# 
# xx <- seq(-1e5,1e5,len=250)
# 
# cdf <- QFcdf(evals=evals)
# test<-QFBounds(cdf = cdf,bound.args = list("Er"=Er,"nu"=nu,"L"=L))
# yy<-test(xx-10,xx+10,log.p = F,lower.tail = T) # if provided, ql must be the same length as qu.  by default, qu = rep(NULL,length(ql))
# 
# #yy<-sapply(xx,function(x){QFIntegrate3(x,x,t.cdf,conc.ineqs,lower.tail = T,log.p = F)})
# 
# plot(xx,cdf(xx),type="l")
# lines(xx,yy[1,],col="blue",lty=2)
# lines(xx,yy[2,],col="red",lty=2)
# plot(xx,cdf(xx,density = T,log.p = F),type="l")
# 
# 
# plot(xx,-cdf(xx,log.p = T,lower.tail = T)/log(10),type="l")
# points(xx,-log(yy[2,])/log(10),type="l",col="red")
# points(xx,-log(yy[1,])/log(10),type="l",col="blue")
# 
# plot(xx,-cdf(xx,log.p = T,lower.tail = F)/log(10),type="l",ylim=c(0,20))
# points(xx,-log1p(-yy[2,])/log(10),type="l",col="red")
# points(xx,-log1p(-yy[1,])/log(10),type="l",col="blue")
# 
# 
# 
# conc.ineqs$int.h2.const(qu,Inf,qu)-conc.ineqs$int.h2.expx(qu,Inf,qu,t.cdf$a.r,t.cdf$b.r) 
# 
# print(-log1p(-QFIntegrate3(xx[which.min(diff(yy[1,]))+1],xx[which.min(diff(yy[1,]))+1],t.cdf,conc.ineqs))/log(10),digits=22)
# 
# xxx<-seq(xx[which.min(diff(yy[1,]))-4],xx[which.min(diff(yy[1,]))+4],len=500)
# zzz<-sapply(xxx,function(x){-log1p(-(GaussQuadCDF(t.cdf$x[1], t.cdf$x[t.cdf$n], T, F.cdf,x, x, conc.ineqs)+conc.ineqs$int.h1.const(t.cdf$x[t.cdf$n],x,x)-conc.ineqs$int.h1.expx(t.cdf$x[t.cdf$n],x,x,t.cdf$a.r,t.cdf$b.r)))/log(10) })
# 
# zzz1<-sapply(xxx,function(x){-log1p(-(GaussQuadCDF(t.cdf$x[1], t.cdf$x[t.cdf$n], T, F.cdf,x, x, conc.ineqs)))/log(10) })
# 
# zzz2<-sapply(xxx,function(x){-log1p(-(conc.ineqs$int.h1.const(t.cdf$x[t.cdf$n],x,x)-conc.ineqs$int.h1.expx(t.cdf$x[t.cdf$n],x,x,t.cdf$a.r,t.cdf$b.r)))/log(10) })
# 
# zzz3<-sapply(xxx,function(x){GaussQuadCDF(t.cdf$x[1], t.cdf$x[t.cdf$n], T, F.cdf,x, x, conc.ineqs)+(conc.ineqs$int.h1.const(t.cdf$x[t.cdf$n],x,x)-conc.ineqs$int.h1.expx(t.cdf$x[t.cdf$n],x,x,t.cdf$a.r,t.cdf$b.r)) })
# 
# #zzz<-sapply(xxx,function(x){-log10(GaussQuadCDF(t.cdf$x[1], t.cdf$x[t.cdf$n], T, F.cdf,x, x, conc.ineqs))})
# plot(xxx,zzz,type="l",ylim=c(0,12))
# lines(xxx,zzz1,col="green")
# lines(xxx,zzz2,col="orange")
# lines(xxx,zzz3,col="blue")
# plot(xxx,zzz3)
# 
# zzz4<-sapply(xx,function(x){conc.ineqs$int.h2.const(x,Inf,x)})
# zzz5<-sapply(xx,function(x){conc.ineqs$int.h2.expx(x,Inf,x,t.cdf$a.r,t.cdf$b.r)})
# 
# zzz6<-sapply(seq(0,10000,len=200),function(x){ -conc.ineqs$int.h2.expx(x,Inf,x,t.cdf$a.r,t.cdf$b.r)})# + conc.ineqs$int.h2.const(x,Inf,x) + (1-conc.ineqs$l(-conc.ineqs$c2))*F.cdf(x) })
# 
# plot(seq(0,10000,len=200),-log10(-zzz6))
# 
# 
# 
# abline(h=-log10(GaussQuadCDF(t.cdf$x[1], t.cdf$x[t.cdf$n], T, F.cdf,xx[which.min(diff(yy[1,]))+1], xx[which.min(diff(yy[1,]))+1], conc.ineqs)))
# abline(v=xx[which.min(diff(yy[1,]))+1])
# yyy<-sapply(xxx,function(x){GaussQuadCDF(t.cdf$x[1], t.cdf$x[t.cdf$n], T, F.cdf, x, x, conc.ineqs)})
# plot(-log10(yyy))
# plot(xxx,yyy,type="l")
# 
# yyy<-sapply(xxx,function(x){conc.ineqs$int.h1.explogx(0,t.cdf$x[1],x,t.cdf$a.l,t.cdf$b.l)})
# plot(xxx,yyy)
# 
# plot(xxx,yyy+zzz)
# abline()
# 
# 
# lines(xx,yy[2,],col="red")
# which(yy[2,]<F.cdf(xx))
# which(yy[1,]>F.cdf(xx))
# 
# 
# 
# 
# 
# # 
# u <- function(x,nu,L) ifelse(x<0,1,ifelse(x < nu*L, exp(-0.5*(x^2)/nu), exp(0.5*nu*L^2-L*x)))
# u.prime <- function(x,nu,L) ifelse(x<0,0,ifelse(x < nu*L, -x*exp(-0.5*(x^2)/nu)/nu,-L*exp(0.5*nu*L^2-L*x)))
# # h1 <- function(t,a,nu,L) -u.prime(a-t,nu,L)
# # h2 <- function(t,a,nu,L) -u.prime(t-a,nu,L)
# l <- function(x,nu,L) ifelse(x>0,1,ifelse(x > -nu*L, exp(-0.5*(x^2)/nu), exp(0.5*nu*L^2+L*x)))
# # l.prime <- function(x,nu,L) ifelse(x>0,0,ifelse(x > -nu*L, -x*exp(-0.5*(x^2)/nu)/nu, L*exp(0.5*nu*L^2+L*x)))
# 
# lower.bound <- function(q,nu,L,Er,F.cdf) (1-u(-Er,nu,L))*F.cdf(q) - integrate(function(x,nu,L,Er){F.cdf(q-Er-x)*u.prime(x,nu,L)},lower=max(-Er,0),upper=Inf,nu=nu,L=L,Er=Er,rel.tol = 1e-12)$value
# upper.bound <- function(q,nu,L,Er,F.cdf) (1-u(Er,nu,L))*F.cdf(q) - integrate(function(x,nu,L,Er){F.cdf(q-Er+x)*u.prime(x,nu,L)},lower=max(Er,0),upper=Inf,nu=nu,L=L,Er=Er,rel.tol=1e-12)$value
# 
# lower.bound2 <- function(ql,nu,L,c1,F.cdf){
#   components <- rep(0,3)
#   components[1] <- (1-conc.ineqs$u(-c1))*F.cdf(ql)
#   components[2] <- integrate( function(t,ql){conc.ineqs$h1(ql,t)*F.cdf(t)}, lower=-Inf,upper=min(ql,ql-c1-sqrt(nu)),ql=ql)$value
#   if(ql > (ql-c1-sqrt(nu)) ){
#     components[3] <- integrate( function(t,ql){conc.ineqs$h1(ql,t)*F.cdf(t)}, lower=ql-c1-sqrt(nu),upper=min(ql,ql-c1),ql=ql)$value
#   }
#   sum(components)
# }
# 
# upper.bound2 <- function(qu,nu,L,c2,F.cdf){
#   components <- rep(0,3)
#   components[1] <- (1-conc.ineqs$l(-c2))*F.cdf(qu)
#   if(qu < (qu-c2+sqrt(nu)) ){
#     components[2] <- integrate( function(t,qu){conc.ineqs$h2(qu,t)*F.cdf(t)}, lower=max(qu,qu-c2),upper=(qu-c2+sqrt(nu)),qu=qu)$value
#   }
#     components[3] <- integrate( function(t,qu){conc.ineqs$h2(qu,t)*F.cdf(t)}, lower=(qu-c2+sqrt(nu)),upper=Inf,qu=qu)$value
#   sum(components)
# }
# 
# # system.time(lower.bound(-30,nu,L))
# # system.time(upper.bound(-30,nu,L))
# #
# # test(10)
# # (1-u(-Er))*test(-20)
# # unlist(integrate(function(x){test(x)*u.prime((-20-Er)-x)},lower=-Inf,upper=-20,rel.tol = 1e-10))
# 
# wrapped.lower.bound<-function(x,nu,L,Er,F.cdf){tryCatch(
#   {lower.bound(x,nu,L,Er,F.cdf)},
#   error = function(e) NA,
#   warning = function(w) NA)}
# 
# wrapped.upper.bound<-function(x,nu,L,Er,F.cdf){tryCatch(
#   {upper.bound(x,nu,L,Er,F.cdf)},
#   error = function(e) NA,
#   warning = function(w) NA)}
# 
# yyy.upper <- sapply(xx,wrapped.upper.bound,nu=nu,L=L,Er=Er, F.cdf = F.cdf)
# yyy.lower <- sapply(xx,wrapped.lower.bound,nu=nu,L=L,Er=Er, F.cdf = F.cdf)
# 
# 
# 
# wrapped.lower.bound2<-function(ql,nu,L,c1,F.cdf){tryCatch(
#   {lower.bound2(ql,nu,L,c1,F.cdf)},
#   error = function(e) NA,
#   warning = function(w) NA)}
# 
# wrapped.upper.bound2<-function(qu,nu,L,c2,F.cdf){tryCatch(
#   {upper.bound2(qu,nu,L,c2,F.cdf)},
#   error = function(e) NA,
#   warning = function(w) NA)}
# 
# yyy.upper2 <- sapply(xx,wrapped.upper.bound2,nu=nu,L=L,c2=Er, F.cdf = F.cdf)
# yyy.lower2 <- sapply(xx,wrapped.lower.bound2,nu=nu,L=L,c1=Er, F.cdf = F.cdf)
# 
# lower.bound2(20,nu,L,Er,F.cdf)
# GaussQuadCDF(-Inf,20,lower.bound = T,F.cdf,20,20,conc.ineqs)
# 
# 
# plot(xx,F.cdf(xx),type="l")
# points(xx,yyy.upper2,type="l",col="red",lty=2)
# points(xx,yyy.lower2,type="l",col="blue",lty=2)
# 
# plot(xx,-F.cdf(xx,log.p = T)/log(10),type="l")
# points(xx,-log10(yyy.upper2),type="l",col="red",lty=2)
# points(xx,-log10(yyy.lower2),type="l",col="blue",lty=2)
# 
# plot(xx,-F.cdf(xx,log.p = T,lower.tail = F)/log(10),type="l")
# points(xx,-log10(1-yyy.upper2),type="l",col="red",lty=2)
# points(xx,-log10(1-yyy.lower2),type="l",col="blue",lty=2)
# 
# any(F.cdf(xx)<yyy.lower2)
# any(F.cdf(xx)>yyy.upper2)
# 
# plot(xxx,log10(test(xxx)-yyy.lower),col="blue",ylab="log10 Bound Distance from T CDF")
# points(xxx,log10(yyy.upper-test(xxx)),col="red")
# 
# #Testing analytic ints
# 
# evals <- abs(c(50,20,10,-9,-8,7,-6,5,4,3,2,1))
# 
# F.cdf <- QFcdf(evals=evals,n=2^16-1)
# 
# nu <- 10
# L <- 1
# Er <- 0.5
# conc.ineqs<-WrapConcIneq.identity(Er,Er,nu,L)
# 
# xx <- seq(-1000,2000,len=1000)
# plot(xx,F.cdf(xx),type="l")
# 
# xx<-seq(0,7,len=100)
# 
# yy<-conc.ineqs$int.h1.expx(lower = 0, upper = 7,q = 3,a=t.cdf$a.l,b = t.cdf$b.l)
# 
# plot(xx,yy,type="l")
# 
# plot(xx,conc.ineqs$d.const(xx,12),type="l")
# plot(xx,conc.ineqs$d.expx(xx,upper=12,a=0.1,b=-0.2),type="l")
# plot(xx,conc.ineqs$d.explogx(xx,upper=12,a=0.1,b=-0.2),type="l")
# 
# 
# # Maybe no longer need this: If t.cdf is NSD, convert it to PSD (just need to flip results at the end)
# was.nsd <- F
# if(t.cdf$type=="nsd"){
#   was.nsd <- T
#   t.cdf$x <- -rev(t.cdf$x)
#   t.cdf$y <- 1-rev(t.cdf$y)
#   t.cdf$limit.r <- Inf
#   t.cdf$limit.l <- 0
#   a.r.2 <- t.cdf$a.l
#   b.r.2 <- -t.cdf$b.l
#   t.cdf$a.l <- t.cdf$a.r
#   t.cdf$b.l <- t.cdf$b.r
#   t.cdf$a.r <- a.r.2
#   t.cdf$b.r <- b.r.2
# }
# 
# 
# 
# 
# 
# 
# #
# #
# # QFIntegrate<-function(q, t.pdf, log.conc.ineq.left,log.conc.ineq.right, L, d, psi, nu){
# #
# #   upper.components<-lower.components<-rep(0,4)
# #
# #   if(q<t.pdf$x[t.pdf$index.l]){
# #
# #     upper.components[1] <- ExpInt(t.pdf$a.l,t.pdf$b.l,-Inf,q) # upper bound on CDF
# #     upper.components[2] <- ExpBoundInt(t.pdf$a.l,t.pdf$b.l,q,t.pdf$x[t.pdf$index.l],log.conc.ineq.right) # gap from q to x[index.l] * bound
# #     upper.components[3] <- ExpBoundInt(t.pdf$a.r,t.pdf$b.r,t.pdf$x[t.pdf$index.r],Inf,log.conc.ineq.right) # x[index.r] to Inf * bound
# #     upper.components[4] <- boole(t.pdf$x[t.pdf$index.l], t.pdf$x[t.pdf$index.l:t.pdf$index.r], exp(log(t.pdf$y[t.pdf$index.l:t.pdf$index.r])+log.conc.ineq.right(t.pdf$x[t.pdf$index.l:t.pdf$index.r],q,L,d,psi,nu)), int.to.right = T) # Middle * bound
# #
# #     lower.components[1] <- boole(t.pdf$x[t.pdf$index.l],t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r],int.to.right = T,one.minus = T) # Middle
# #     lower.components[2] <- -ExpInt(t.pdf$a.l,t.pdf$b.l,q,t.pdf$x[t.pdf$index.l]) # Left #
# #     lower.components[3] <- -ExpInt(t.pdf$a.r,t.pdf$b.r,t.pdf$x[t.pdf$index.r],Inf) # Right  # WE COULD REPLEACE THESE THREE TERMS WITH A LOWER BOUND ON CDF ON LEFT TAIL
# #     lower.components[4] <- -ExpBoundInt(t.pdf$a.l,t.pdf$b.l,-Inf,q,log.conc.ineq.left) # -Inf to q * bound
# #
# #   }
# #
# #   if(q>t.pdf$x[t.pdf$index.r]){
# #
# #     upper.components[1] <- ExpInt(t.pdf$a.l,t.pdf$b.l,-Inf,t.pdf$x[t.pdf$index.l]) # Left
# #     upper.components[2] <- boole(t.pdf$x[t.pdf$index.l],t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r],int.to.right = T) # Middle
# #     upper.components[3] <- ExpInt(t.pdf$a.r,t.pdf$b.r,t.pdf$x[t.pdf$index.r],q) # Right
# #     upper.components[4] <- ExpBoundInt(t.pdf$a.r,t.pdf$b.r,q,Inf,log.conc.ineq.right) # q to Inf * bound
# #
# #
# #     lower.components[1] <- boole(t.pdf$x[t.pdf$index.l],t.pdf$x[t.pdf$index.l:t.pdf$index.r],exp(log(t.pdf$y[t.pdf$index.l:t.pdf$index.r])+log.conc.ineq.left(t.pdf$x[t.pdf$index.l:t.pdf$index.r],q,L,d,psi,nu)),int.to.right = T,one.minus = T) # Middle* bound
# #     lower.components[2] <- -ExpInt(t.pdf$a.r,t.pdf$b.r,q,Inf) # Right
# #     lower.components[3] <- -ExpBoundInt(t.pdf$a.r,t.pdf$b.r,t.pdf$x[t.pdf$index.r],q,log.conc.ineq.left) # Right * bound
# #     lower.components[4] <- -ExpBoundInt(t.pdf$a.l,t.pdf$b.l,-Inf,t.pdf$x[t.pdf$index.l],log.conc.ineq.left) # Left * bound
# #
# #   }
# #
# #   if( t.pdf$x[t.pdf$index.l]<=q & q<=t.pdf$x[t.pdf$index.r] ){
# #
# #
# #     upper.components[1] <- ExpInt(t.pdf$a.l,t.pdf$b.l,-Inf,t.pdf$x[t.pdf$index.l]) # Left
# #     upper.components[2] <- boole(q,t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r],int.to.right = F) # Middle Left
# #     upper.components[3] <- boole(q,t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r]*exp(log.conc.ineq.right(t.pdf$x[t.pdf$index.l:t.pdf$index.r],q,L,d,psi,nu)),int.to.right = T) # Middle Right * Bound
# #     upper.components[4] <- ExpBoundInt(t.pdf$a.r,t.pdf$b.r,t.pdf$x[t.pdf$index.r],Inf,log.conc.ineq.right) # Right * Bound
# #
# #     lower.components[1] <- boole(q,t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r],int.to.right = T,one.minus = T) # Middle Right
# #     lower.components[2] <- -ExpInt(t.pdf$a.r,t.pdf$b.r,t.pdf$x[t.pdf$index.r],Inf) # Right
# #     lower.components[3] <- -boole(q,t.pdf$x[t.pdf$index.l:t.pdf$index.r],t.pdf$y[t.pdf$index.l:t.pdf$index.r]*exp(log.conc.ineq.left(t.pdf$x[t.pdf$index.l:t.pdf$index.r],q,L,d,psi,nu)),int.to.right = F) # Middle Left
# #     lower.components[4] <- -ExpBoundInt(t.pdf$a.l,t.pdf$b.l,-Inf,t.pdf$x[t.pdf$index.l],log.conc.ineq.left) # Left * bound
# #
# #   }
# #
# #   c(max(neumaierSum(lower.components), 0),min(neumaierSum(upper.components), 1))
# #
# # }
# 
# 
# 
# 
# GaussInt<-function(a,b,from,to,k,nu){
#   
#   if(from==to){return(-Inf)}
#   if(to<from){stop("from must be less than or equal to to")}
#   if(b==0){stop("b cannot be zero or else this is just an integral of a constant function")}
#   if(all(is.infinite(c(from,to)))){stop("both from and psi cannot be infinite")}
#   if(is.infinite(c(from)) & from>0){stop("from cannot be positive infinity")}
#   if(is.infinite(c(to)) & to<0){stop("from cannot be negative infinity")}
#   
#   prefix <- a + b^2*nu/2 - k*b + 0.5*log(2*pi*nu)
#   
#   if(is.infinite(from)){
#     if(b>0){ stop("b cannot be positive when from is negative infinity --> integral does not converge")}
#     # b is negative
#     return( prefix + pnorm((to-(k-b*nu))/sqrt(nu),log.p = T) )
#   }
#   
#   if(is.infinite(to)){
#     if(b<0){ stop("b cannot be negative when to is infinity --> integral does not converge")}
#     # b is positive
#     return( prefix + pnorm((from-(k-b*nu))/sqrt(nu),log.p = T,lower.tail = F) )
#   }
#   
#   
#   prefix + pnorm((to-(k-b*nu))/sqrt(nu),log.p = T)
#   + log(-expm1( pnorm((from-(k-b*nu))/sqrt(nu),log.p = T) - pnorm((to-(k-b*nu))/sqrt(nu),log.p = T) ))
# }
# 
# # # Looks like GaussInt is not correct
# xxx<-seq(-100,100,len=1000)
# yyy<-sapply(xxx,function(z,a,b,k,nu){exp(-0.5*((k-z)^2)/nu+a-b*z)},a=1,b=4,k=0,nu=20)
# 
# plot(xxx,yyy)
# 
# exp(GaussInt(1,4,from=-100,to=100,0,20))
# integrate(function(z,a,b,k,nu){exp(-0.5*((k-z)^2)/nu+a-b*z)},-100,100,a=1,b=4,k=0,nu=20)$value
# 
# 
# 
# AnalyticBoundInt<-function(a,b,from,to,right=T,q.=q,L.=L,d.=d,psi.=psi,nu.=nu){
#   
#   if(from==to){return(-Inf)}
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
#   # This must integrate exp(a-bz) * bound where bound may be in one of three stages (flat, gauss, or exp)
#   
#   k.<-q.-psi.
#   
#   rho <- as.integer(right)
#   
#   if(right==T){
#     const.interval <- c(-Inf, k.)
#     gauss.interval <- c(k., nu./(4*d.)+k.)
#     exp.interval <- c(nu./(4*d.)+k., Inf)
#   }else{
#     const.interval <- c(k., Inf)
#     gauss.interval <- c(k.-nu./(4*d.), k.)
#     exp.interval <- c(-Inf, k.-nu./(4*d.))
#   }
#   
#   const.component <- gauss.component <- exp.component <- -Inf
#   
#   ##### Constant Component #####
#   
#   # Interval inside region of integration
#   if(from <= const.interval[1] & const.interval[2]<=to){
#     const.component <- ExpInt(a,b,const.interval[1],const.interval[2])
#   }
#   
#   # Interval overlaps left part of region of integration
#   if(const.interval[1] < from & from < const.interval[2] & const.interval[2] <= to){
#     const.component <- ExpInt(a,b,from,const.interval[2])
#   }
#   
#   # Interval covers entire region of integration
#   if(const.interval[1] < from & to < const.interval[2] ){
#     const.component <- ExpInt(a,b,from,to)
#   }
#   
#   # Interval overlaps right part of region of integration
#   if(from <= const.interval[1] & const.interval[1] < to & to < const.interval[2]){
#     const.component <- ExpInt(a,b,const.interval[1],to)
#   }
#   
#   ##### Gauss Component #####
#   
#   # Interval inside region of integration
#   if(from <= gauss.interval[1] & gauss.interval[2]<=to){
#     gauss.component <- GaussInt(a,b,gauss.interval[1],gauss.interval[2],k.,nu.)
#   }
#   
#   # Interval overlaps left part of region of integration
#   if(gauss.interval[1] < from & from < gauss.interval[2] & gauss.interval[2] <= to){
#     gauss.component <- GaussInt(a,b,from,gauss.interval[2],k.,nu.)
#   }
#   
#   # Interval covers entire region of integration
#   if(gauss.interval[1] < from & to < gauss.interval[2] ){
#     gauss.component <- GaussInt(a,b,from,to,k.,nu.)
#   }
#   
#   # Interval overlaps right part of region of integration
#   if(from <= gauss.interval[1] & gauss.interval[1] < to & to < gauss.interval[2]){
#     gauss.component <- GaussInt(a,b,gauss.interval[1],to,k.,nu.)
#   }
#   
#   
#   ##### Exp Component #####
#   
#   # Interval inside region of integration
#   if(from <= exp.interval[1] & exp.interval[2]<=to){
#     exp.component <- ExpInt(a+0.5*nu./((4*d.)^2)+rho*k./(4*d.),b+rho/(4*d.),exp.interval[1],exp.interval[2])
#   }
#   
#   # Interval overlaps left part of region of integration
#   if(exp.interval[1] < from & from < exp.interval[2] & exp.interval[2] <= to){
#     exp.component <- ExpInt(a+0.5*nu./((4*d.)^2)+rho*k./(4*d.),b+rho/(4*d.),from,exp.interval[2])
#   }
#   
#   # Interval covers entire region of integration
#   if(exp.interval[1] < from & to < exp.interval[2] ){
#     exp.component <- ExpInt(a+0.5*nu./((4*d.)^2)+rho*k./(4*d.),b+rho/(4*d.),from,to)
#   }
#   
#   # Interval overlaps right part of region of integration
#   if(from <= exp.interval[1] & exp.interval[1] < to & to < exp.interval[2]){
#     exp.component <- ExpInt(a+0.5*nu./((4*d.)^2)+rho*k./(4*d.),b+rho/(4*d.),exp.interval[1],to)
#   }
#   
#   # Components are here in log space
#   temp <- c(const.component,gauss.component,exp.component)
#   norm.const <- max(temp)
#   # Return result in the probability space
#   neumaierSum(c(log(neumaierSum(exp(temp-norm.const))),norm.const))
# }
# 
# 
# # xx<-seq(0,10000,len=1000)
# # yy<-rep(0,length(xx))
# # for(i in 1:length(xx)) yy[i]<--log10(ExpBoundInt(t.pdf$a.l,t.pdf$b.l,-100,xx[i],log.conc.ineq = log.conc.ineq.l))
# # plot(xx,yy,col="red",type="l")
# #
# # xx<-seq(-1540,-1480,len=1000)
# # yy<-rep(0,length(xx))
# # for(i in 1:length(xx)) yy[i]<--AnalyticBoundInt(t.pdf$a.l,t.pdf$b.l,-1550,xx[i],right = F)/log(10)
# # plot(xx,yy,type="l",col="blue")
# # abline(v=-1501)
# # abline(v=-1526)
# #
# #
# # print(-AnalyticBoundInt(t.pdf$a.r,t.pdf$b.r,2410.521,Inf,right=T,q.=q,L.=L,d.=d,psi.=psi,nu.=nu)/log(10),digits=22)
# #
# # xx<-seq(2410.521,4000,len=1000)
# # yy<-rep(0,length(xx))
# # for(i in 1:length(xx)) yy[i]<--AnalyticBoundInt(t.pdf$a.r,t.pdf$b.r,xx[i],Inf,right = T)/log(10)
# # plot(xx,yy,type="l")
# #
# # xx<-seq(-10000,10000,len=1000)
# # yy<-rep(0,length(xx))
# # for(i in 1:length(xx)) yy[i]<--AnalyticBoundInt(t.pdf$a.l,t.pdf$b.l,-Inf,xx[i],right = F)/log(10)
# # plot(xx,yy,type="l")
# #
# #
# 
# 
# 
# #
# # ExpBoundInt<-function(a,b,from,to,log.conc.ineq,q.=q,L.=L,d.=d,psi.=psi,nu.=nu){
# #
# #   # This performs the integral $ \int_\from^\to exp(a - b*z + log.conc.ineq(z,q,L,d,psi,nu)) dz $ for any sign of a,b.
# #
# #   if(from==to){return(0)}
# #   if(to<from){stop("from must be less than or equal to to")}
# #   if(b==0){stop("b cannot be zero or else this is just an integral of a constant function")}
# #   if(all(is.infinite(c(from,to)))){stop("both from and psi cannot be infinite")}
# #   if(is.infinite(c(from)) & from>0){stop("from cannot be positive infinity")}
# #   if(is.infinite(c(to)) & to<0){stop("psi cannot be negative infinity")}
# #
# #   if(is.infinite(from) & b>0){ stop("b cannot be positive when from is negative infinity --> integral does not converge")}
# #
# #   if(is.infinite(to) & b<0){stop("b cannot be negative when to is infinity --> integral does not converge")}
# #
# #   integrate(f=function(z){
# #     exp(a - b*z + log.conc.ineq(z,q.,L.,d.,psi.,nu.))
# #   }, lower = from, upper = to, rel.tol = .Machine$double.eps)$value
# #
# #
# # }
# #
# 
# # This is h_k(t) meant to be used in obtaining upper bounds
# log.conc.ineq.r <-function(z,q,L,d,psi,nu){
#   k<-q-psi
#   ifelse(z <= k ,0,ifelse(z <= nu/(4*d)+k,-0.5*((k-z)^2)/nu, 0.5*nu/((4*d)^2) + (k-z)/(4*d)))
# }
# 
# # This is h_-k(-t) meant to be used in obtaining lower bounds
# log.conc.ineq.l<-function(z,q,L,d,psi,nu){
#   log.conc.ineq.r(-z,-q,L,d,-psi,nu)
# }
# 
# # Concentration Inequality bounds look right
# # xx<-seq(q-psi-30,q-psi+nu/(4*d)+30,len=1000)
# # plot(xx,-log.conc.ineq.r(xx,q,L,d,psi,nu)/log(10),type="l")
# # abline(v=q-psi)
# # abline(v=q-psi+nu/(4*d))
# 
# #
# # plot(xx,-log.conc.ineq.l(xx,q,L,d,psi,nu)/log(10),type="l")
# # abline(v=(q-psi))
# # abline(v= c(q-psi)-nu/(4*d))
# 
# 
# # # Analytic Exponential Integrals Look Good
# # xx<-seq(-2000,2000,len=1000)
# # yy<-xx
# # for(i in 1:length(xx)) yy[i]<-log10(ExpInt(t.pdf$a.r,t.pdf$b.r,xx[i],Inf))
# # plot(xx,yy)
# 
# 
# 
# 
# ######
# # Example Code Here
# #######
# 
# 
# # Gauss Quadrature Exponential * Bound Integral looks....bad
# #
# # require(QForm)
# # require(PreciseSums)
# #
# # evals<- c(-10:-1,1:10)
# # ncps<-1:20
# #
# # q<--1500
# # L<-1
# # d<-1
# # psi<-1
# # nu<-1000
# #
# # t.cdf<-Tcdf_new(evals,ncps)
# #
# # plot(t.pdf$x,t.pdf$y,type="l")
# #
# #
# # t.pdf<-Tpdf(evals,ncps)
# #
# # QFIntegrate2(-4000,t.pdf,log.conc.ineq.l,log.conc.ineq.r,L = L,d=d,psi=psi,nu=nu)
# 
# #
# #
# # xx<-seq(-2000,2000,len=1000)
# # yy<-matrix(nrow=length(xx),ncol=2)
# # system.time(for(i in 1:length(xx)){
# #   yy[i,]<-QFIntegrate2(xx[i],t.pdf,log.conc.ineq.l,log.conc.ineq.r,L = L,d=d,psi=psi,nu=nu)
# # })
# #
# #
# # plot(xx,yy[,1],type="l",col="red")
# # points(xx,yy[,2],type="l",col="blue")
# #
# #
# # plot(xx,-log1p(-yy[,1])/log(10),type="l",col="red",ylim=c(0,5))
# # lines(xx,-log1p(-yy[,2])/log(10),col="blue")
# #
# # plot(xx,-log(yy[,1])/log(10),type="l",col="red",ylim=c(0,5))
# # lines(xx,-log(yy[,2])/log(10),col="blue")
# #
# #
# #
# # xxx<-seq(-4000,4000,len=100)
# # yyy<-matrix(nrow=length(xxx),ncol=3)
# # system.time(for(i in 1:length(xxx)){
# #   yyy[i,]<-as.numeric(as.vector(QFIntBounds2(obs = xxx[i],evals = evals,ncps=ncps,E_R=psi,nu = nu,resid.op.norm.bd = d)))
# #     })
# #
# # plot(xxx,yyy[,1],type="l",col="red")
# # points(xxx,yyy[,2],type="l",col="blue")
# #
# #
# # plot(xxx,-log1p(-yyy[,1])/log(10),type="l",col="red",ylim=c(0,18))
# # lines(xxx,-log1p(-yyy[,2])/log(10),col="blue")
# #
# # plot(xxx,-log(yyy[,1])/log(10),type="l",col="red",ylim=c(0,18))
# # lines(xxx,-log(yyy[,2])/log(10),col="blue")
# #
# # # Confirm that it's not a bug in how we combine results or a numerical stability problem combining results
# # # Confirm it's not a prob in Analytic integral logic
# # # See if there's a problem in a particular term
# # # If not any of these then we have a problem I think in boole:
# # # check for a problem in boole...where is it hitting 1 or 0 ?
# # # Could rewrite boole function to do integral of exp(g) where we standardize g to keep exp in a better interval
# # # Could go with quad precision
# #
# #
# #
# #
# #
# #
# #
# #
# # xx<-seq(-500,4000,len=1000)
# # yy<-rep(0,length(yy))
# # for(i in 1:length(xx)) yy[i]<--log10(ExpBoundInt(t.pdf$a.r,t.pdf$b.r,xx[i],Inf,log.conc.ineq = log.conc.ineq.l))
# # plot(xx,yy)
# #
# #
# # xx<-seq(-500,4000,len=1000)
# # yy<-rep(0,length(yy))
# # for(i in 1:length(xx)) yy[i]<--log10(AnalyticBoundInt(t.pdf$a.r,t.pdf$b.r,xx[i],Inf,right = F))
# # plot(xx,yy)
# #
# # xx<-seq(-1000,4000,len=1000)
# # yy<-rep(0,length(yy))
# # for(i in 1:length(xx)) yy[i]<--log10(ExpBoundInt(t.pdf$a.l,t.pdf$b.l,-Inf,xx[i],log.conc.ineq = log.conc.ineq.r))
# # plot(xx,yy)
# #
# # ######
# # ### Looks like Analytic Bound does the wrong thing in this case...need to check this out
# # #######
# # xx<-seq(-1000,4000,len=1000)
# # yy<-rep(0,length(yy))
# # for(i in 1:length(xx)) yy[i]<- AnalyticBoundInt(t.pdf$a.l,t.pdf$b.l,-Inf,xx[i])
# # plot(xx,yy)
# #
# # xx<-seq(-2000,2000,len=1000)
# # yy<-rep(0,length(yy))
# # for(i in 1:length(xx)) yy[i]<--log10(ExpBoundInt(t.pdf$a.l,t.pdf$b.l,-Inf,xx[i],log.conc.ineq = log.conc.ineq.l))
# # plot(xx,yy,type="l")
# #
# # ############
# # # This one is also different!!!
# # ##########
# #
# # xx<-seq(-2000,2000,len=1000)
# # yy<-rep(0,length(yy))
# # for(i in 1:length(xx)) yy[i]<--log10(AnalyticBoundInt(t.pdf$a.l,t.pdf$b.l,-Inf,xx[i],right = F))
# # plot(xx,yy,type="l")
# #
# #
# #
# # #######
# #
# #
# #
# #
# # QFIntegrate2(q,t.pdf,log.conc.ineq.l,log.conc.ineq.r,L = L,d=d,psi=psi,nu=nu)
# #
# 
# 
# a<-1e30
# b<-1e10
# la <- log(a)
# lb <- log(b)
# 
# print(la+log(-expm1(lb-la)),digits=22)
# print(log(a)+log1p(-b/a),digits=22)
# 
# print(log(a-b),digits=22)
# 
# 
# 
