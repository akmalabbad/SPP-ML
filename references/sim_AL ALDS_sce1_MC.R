##############################################################
# Simulation study for AL ALDS
##############################################################
library(spatstat)
library(lpSolve)
library(glmnet)
load("~/R Code ALDS/ALDS/bei.covars.rda")
source("~/R Code ALDS/ALDS/sppDantzig_dec19.R")
source("~/R Code ALDS/ALDS/spplasso_dec19.R")
source("~/R Code ALDS/ALDS/find.intercept.R")

###############################
# Initialization
###############################
n.sim<-500; p.sim<-20; beta1<- 2; beta2<- -1
wx=250; wy=125; npts=150; W=owin(c(0,wx),c(0,wy))
kappa=wy/(2*area.owin(W)); scale=5e-2*wx
rho<-exp(beta1*bei.covars$elev+beta2*bei.covars$grad)
beta0<-find.intercept(rho,wx,wy,npts) #see function find intercept
Ifunc<-exp(beta0+beta1*bei.covars$elev+beta2*bei.covars$grad)
trend=as.formula("~elev+grad+Al+B+Ca+Cu+Fe+K+Mg+Mn+P+Zn+N+N.min.+pH+x3x4+x3x5+x3x6+x3x7+x3x8")

fit<- coef<- lambda<-time<-n.lambda<-list()
n.penalty<-2 #AL and ALDS

for (i in 1:n.penalty){
  coef[[i]]=matrix(0,p.sim+1,n.sim)
  lambda[[i]]<-time[[i]]<-n.lambda[[i]]<-rep(0,n.sim)
}

###############################
# Simulation
###############################
for (i in 1:n.sim){
  cat('\r nsim=',i)
  pp <- rThomas(kappa, scale = scale, mu=Ifunc/(kappa*2e-6*area.owin(W)), win=W)
  time[[2]][i]<-system.time(fit[[2]]<- sppDantzig(pp,covariates=bei.covars[1:p.sim],trend=trend,epsilon=1e-5))[3]
  time[[1]][i]<-system.time(fit[[1]]<- spplasso(pp,covariates=bei.covars[1:p.sim],trend=trend,alpha=1,adaptive=T))[3]
  
  for (j in 1:n.penalty){
    coef[[j]][,i]=fit[[j]]$beta
    lambda[[j]][i]=fit[[j]]$lambda.opt
    n.lambda[[j]][i]=fit[[j]]$n.lambda
  }
}
cat('\n')

###################################
# Computing TPR FPR
###################################
n.t.beta<-2; n.f.beta<- p.sim-n.t.beta
TP<-FP<-rep(0,n.penalty)
selection <- list()

for (i in 1:n.penalty)
{
  selection[[i]] <- matrix(0,p.sim,n.sim)
  selection[[i]][which(coef[[i]][-1,]!=0)]<-1
  TP[i]<-sum(selection[[i]][c(1:n.t.beta),])*100/(n.t.beta*n.sim)
  FP[i]<-sum(selection[[i]][c((n.t.beta+1):p.sim),])*100/(n.f.beta*n.sim) 
}

###################################
# Computing Bias SD RMSE
###################################
mean.coef<-list()
bias<-SD<-RMSE<-rep(1e30,n.penalty)
true.beta=c(beta0,beta1,beta2,rep(0,p.sim-2))
for (i in 1:n.penalty){
  mean.coef[[i]]<-apply(coef[[i]],1,mean)
  bias[i]=sqrt(sum((mean.coef[[i]]-true.beta)^2))
  SD[i]=sqrt(sum(apply(coef[[i]],1,var)))
  RMSE[i]=sqrt(sum(apply((coef[[i]]-true.beta)^2,1,mean)))
}

result_tho_1_MC<-rbind(TP,FP,bias,SD,RMSE,c(mean(time[[1]]),mean(time[[2]])),c(mean(lambda[[1]]),mean(lambda[[2]])),c(mean(n.lambda[[1]]),mean(n.lambda[[2]])))