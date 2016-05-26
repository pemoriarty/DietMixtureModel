NLL.mle <- function(par,data){          

  #read in parameters and assign them names
  r_theta=par[1]
  r_thetap1=par[2]
  ms.pres=par[3]
  var.pres=par[4]
  var.abs=par[5]
  betamean=par[6]
  betasd=par[7]
  c_i=par[8]
  
  #calculate beta distribution values from parameter values
  alpha1 <- ((1-betamean)/betasd^2 - 1/betamean)*betamean^2
  alpha2 <- alpha1*(1/betamean-1)
  if(alpha1 <=0){alpha1=0.01}
  if(alpha2 <=0){alpha2=0.01}
  #betamean = Bmean(alpha1,alpha2)
  
  ms.abs<- (ms.pres*r_theta*(betamean + r_thetap1 - c_i - betamean*r_thetap1))/(c_i*(1-r_theta))
  
  #calculate parameters for the gamma distributions from the mean and variance
  shape.pres<-ms.pres^2/var.pres
  scale.pres<-sqrt(var.pres/shape.pres)
  rate.pres<-1/scale.pres
  shape.abs<-ms.abs^2/var.abs
  scale.abs<-sqrt(var.abs/shape.abs)
  rate.abs<-1/scale.abs
  
  
  not.zero.index<-which(data[,1]>0)#subset stomachs that contain the prey
  full.index<-which(data[,1]==1)#subset stomachs that only contain the prey
  beta.use.index<-setdiff(not.zero.index,full.index)#which datapoints have a diet fraction between 0 and 1 
  
  # calculate NLL for empty stomachs
  tmp.data<-data[-not.zero.index,2]
  NLL.m.empty=sum(-dgamma(tmp.data,shape=shape.abs,rate=rate.abs,log=T))
  
  # calculate NLL for stomachs containing the prey
  tmp.data<-data[not.zero.index,2]
  NLL.m.full<-sum(-dgamma(tmp.data,shape=shape.pres,rate=rate.pres,log=T))
  
  # calcualte NLL for diet proportions between 0 and 1
  tmp.data<-data[beta.use.index,1]
  NLL.diet.p<-sum(-dbeta(tmp.data,alpha1,alpha2,log=T))
  
  # calculate NLL of  binomial pdfs
  N.nozero<-max(0,length(not.zero.index))
  N.full<-max(0,length(full.index))
  N.total<-nrow(data)
  NLL.binom=-dbinom(x=N.nozero,size=N.total,prob=r_theta,log=T)
  NLL.full=-dbinom(x=N.full,size=N.nozero,prob=r_thetap1,log=T)
  
  return(NLL.m.empty+NLL.m.full+NLL.diet.p+NLL.binom+NLL.full)#return total negative log likelihood
}