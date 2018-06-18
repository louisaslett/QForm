# #
# # Testing Accuracy of CompQuadForm
# #
#
# n.reps<-1000000
# K<-2
# non_square.ncps<-c(-1,20)
# true_samps<-colSums(((matrix(rnorm(n = n.reps*K),nrow=K,ncol=n.reps)+non_square.ncps)^2)*c(50,-2))
#
# xx<-seq(-1500,500,by=1)
#
# truth<-ecdf(true_samps)(xx)
#
# res1<-rep(0,length(xx))
# res2<-rep(0,length(xx))
#
# for(i in 1:length(xx)){
#   res1[i]<-1-davies(q = xx[i],lambda = c(50,-2),delta=(non_square.ncps)^2,lim = 20000, acc = 1e-12)$Qq
#   res2[i]<-1-imhof(q = xx[i],lambda = c(50,-2),delta=(non_square.ncps)^2,limit = 20000, epsrel = 1e-12)$Qq
# }
#
# # Plot Lower Tail
# plot(xx,-log10(res1),type="l",col="red")
# lines(xx,-log10(res2),type="l",col="green")
# lines(xx,-log10(truth),col="blue")
#
# #Plot Upper Tail
# plot(xx,-log10(1-res1),type="l",col="red")
# lines(xx,-log10(1-res2),type="l",col="green")
# lines(xx,-log10(1-truth),col="blue")
#
# # Davies appears to be the most stable method in the package. Imhof seems to have more instability deep in the tails.
#
