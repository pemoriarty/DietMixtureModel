#Pamela Moriarty
#last edited August 17, 2015
#estimates parameters for the GULPS Model
find.mle <- function(data,p.bin,p.full,mean.pres,var.pres,var.abs,beta.mean,beta.sd,p.pop){
  #data = 2 column dataframe, where first column is proportion of each stomach made up of the prey type and the second column is the total stomach mass
  #p.bin = starting value for probability a stomach contains the prey type
  #p.full = starting value for probability a stomach that contains the prey only contains the prey type
  #mean.pres = starting value for mean stomach mass of those that contain the prey type
  #var.pres = starting value for variance of stomach masses of those that contain the prey type
  #var.abs = starting value for variance of stomach masses of those that do not contain the prey type
  #beta.mean, beta.sd = starting values for the mean and standard deviation of the beta distribution
  #p.pop = starting value for the total contribution of the prey type to the predator's diet at the population level
  
  require(stats4)
  source('Bmean.R')  
  data.fn <- function() data
    
  NLL.mle <- function(par){          

      #pull in the data
      data <- data.fn()
      p.bin=par[1]
      p.full=par[2]
      mean.pres=par[3]
      var.pres=par[4]
      var.abs=par[5]
      beta.mean=par[6]
      beta.sd=par[7]
      p.pop=par[8]

      #calculate beta distribution parameters
      alpha1 <- ((1-beta.mean)/beta.sd^2 - 1/beta.mean)*beta.mean^2
      alpha2 <- alpha1*(1/beta.mean-1)
      if(alpha1 <=0){alpha1=0.01}
      if(alpha2 <=0){alpha2=0.01}
      pdietmean = Bmean(alpha1,alpha2)
      
      mean.abs<- (mean.pres*p.bin*(pdietmean + p.full - p.pop - pdietmean*p.full))/(p.pop*(1-p.bin))
      
      #calculate parameters for the gamma distributions from the mean and variance
      shape.pres<-mean.pres^2/var.pres
      scale.pres<-sqrt(var.pres/shape.pres)
      rate.pres<-1/scale.pres
      shape.abs<-mean.abs^2/var.abs
      scale.abs<-sqrt(var.abs/shape.abs)
      rate.abs<-1/scale.abs
      
      not.zero.index<-which(data[,1]>0)
      full.index<-which(data[,1]==1)
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
      NLL.binom=-dbinom(x=N.nozero,size=N.total,prob=p.bin,log=T)
      NLL.full=-dbinom(x=N.full,size=N.nozero,prob=p.full,log=T)
      
      return(NLL.m.empty+NLL.m.full+NLL.diet.p+NLL.binom+NLL.full)#return total negative log likelihood
  }

#set bounds for the parameters
lower.bnd<-rep(0.0001,8)
upper.bnd<-rep(100000000,8)
upper.bnd[c(1,2,6,8)]<-rep(.999,1)

require(optimx)
scale.warn <- function(w) {
  if(any(grepl("Parameters or bounds appear to have different scalings.",w))){
    invokeRestart("muffleWarning")}
}

optim.result=withCallingHandlers(optimx(par=c(p.bin,p.full,mean.pres,var.pres,var.abs,beta.mean,beta.sd,p.pop),fn=NLL.mle,lower=lower.bnd,upper=upper.bnd,method="nlminb",control=list(maxit=500)),warning=scale.warn)


return(optim.result)
}