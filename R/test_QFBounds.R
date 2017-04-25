# load("~/Documents/lab/thesis/Chapter3/chapter3_R_code/ABO_kenya_min_mat_concentration_bound_demonstration_workspace.RData")
# #
# # # Everything is input on original scale, y^TBy.
#
#
# xx<-seq(-2e5,2e5,by=1000)
# out<-QFBounds(xx,abo_min_mat,2*p-1,2*sqrt(p*(1-p)),25,log=F)
# plot(xx,ecdf(N*true_samps)(xx),ylim=c(0,1),type="l")
# points(out$obs,out$upper,col="red",type="l")
# points(out$obs,out$lower,col="blue",type="l")
#
# xx<-seq(-10e5,10e5,by=10000)
# yy<- -log(ecdf(N*true_samps)(xx))/log(10)
# out<-QFBounds(xx,abo_min_mat,2*p-1,2*sqrt(p*(1-p)),25,log=T)
# plot(xx,yy,type="l",ylim=c(0,13))
# points(out$obs,-out$upper/log(10),col="red",type="l")
# points(out$obs,-out$lower/log(10),col="blue",type="l")
#
# QFBounds(xx,abo_min_mat,2*p-1,2*sqrt(p*(1-p)),25,log=T)
#
#
#
# out
#
# # require(CompQuadForm)
# # require(RSpectra)
# # evals<-eigs(abo_min_mat,k=25)$values
# # evals<-evals[order(abs(evals),decreasing = T)]
# # eps<-seq(0,10,by=0.01)
# # E_R<-0
# # ncps<-rep(0,length(evals))
# # obs<-2000
# # resid_operator_norm_bound<-abs(evals[length(evals)])
# # QFBounds.ineq.lower(eps[100], obs, evals, ncps, E_R, nu2,resid_operator_norm_bound)
# # xx<-seq(0,40000,by=100)
# # plot(xx,RemainderBoundSupportFinder(xx, resid_operator_norm_bound, ncps, E_R, nu2/N))
#
#
# ####
# # This code below corresponds for the browser for checking that the objective/optimization is working
# #####
#
# QFBounds(-810000,abo_min_mat,2*p-1,2*sqrt(p*(1-p)),25,log=T)
# QFBounds(-820000,abo_min_mat,2*p-1,2*sqrt(p*(1-p)),25,log=T)
#
# eps<-seq(0,int_right,length.out = 1000)
#
# potential_bound<-eps
# for(j in 1:length(eps)){
# potential_bound[j]<-exp(QFBounds.ineq.lower(eps[j],obs = obs,
#                                             evals = evals,
#                                             ncps = ncps,
#                                             E_R = E_R,
#                                             nu2 = nu2,
#                                             resid_operator_norm_bound=resid_operator_norm_bound))
# }
#
# plot(eps,potential_bound)
#
#
