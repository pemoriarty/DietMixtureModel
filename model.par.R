#Pamela Moriarty
#last edited May 10, 2016
#estimate contribution of a prey to a predator's diet using the mixture model
#this takes the data, calculates reasonable starting values for the model to run and returns parameter estimates

model.par <- function(prop,total.mass){
  #prey- vector of the proportion of each stomach made up of the prey type
  #total.mass- vector of the total mass of contents in the stomachs
  source('find.mle.R')

  prey.dat <- cbind(prop,total.mass) 
  prey.dat[,1] <- replace(prey.dat[,1],is.nan(prey.dat[,1]),0)#remove NaNs from the data
  prey.dat[,1] <- replace(prey.dat[,1],is.na(prey.dat[,1]),0)#remove NAs from the data
  prey.dat[,1] <- replace(prey.dat[,1],is.infinite(prey.dat[,1]),0)#remove INFs from the data

  prey.dat.no0 <- subset(prey.dat,prey.dat[,2] > 0)#remove empty stomachs 
  prey.dat.no0 <- subset(prey.dat.no0,!is.nan(prey.dat.no0[,2]))
  prey.dat.no0 <- subset(prey.dat.no0,!is.na(prey.dat.no0[,2]))
  prey.dat.no0 <- as.matrix(prey.dat.no0)

  #calculate starting values for the model and check all 3 cases are present in the data (if not, some parameters can't be estimated as they are irrelevant)
  nprey <- length(which(prey.dat.no0[,1]>0))
  nfull <- length(which(prey.dat.no0[,1]==1))
  n <- nrow(prey.dat.no0)
  r_theta.tmp <- nprey/n
  if(r_theta.tmp==0){return("Error. No stomachs contain the prey. Use a reduced model to estimate the relevant parameters.")}
  r_thetap1.tmp <- nfull/nprey
  if(r_thetap1.tmp==0){return("Error. No stomachs contain only the prey. Use a reduced model to estimate the relevant parameters.")}
  notzero.idx <- setdiff(which(prey.dat.no0[,1] > 0), which(prey.dat.no0[,1] ==1))
  betamean.tmp <- mean(prey.dat.no0[notzero.idx,1])
  if(is.nan(betamean.tmp)){return("Error. No stomachs have a diet fraction of 0 < p < 1. Use a reduced model to estimate the relevant parameters.")}
  beta.sd.tmp <- sd(prey.dat.no0[notzero.idx,1])
  ms.pres.tmp <- mean(prey.dat.no0[which(prey.dat.no0[,1]>0),2])
  sd.pres.tmp <- sd(prey.dat.no0[which(prey.dat.no0[,1]>0),2])
  sd.abs.tmp <- sd(prey.dat.no0[which(prey.dat.no0[,1]==0),2])
  wt.mean <- sum(prey.dat.no0[,1]*prey.dat.no0[,2])/sum(prey.dat.no0[,2])

  ms.abs.tmp<-mean(prey.dat.no0[which(prey.dat.no0[,1]==0),2])
  c_i.tmp<-r_theta.tmp*((1-r_thetap1.tmp)*betamean.tmp*ms.pres.tmp+r_thetap1.tmp*ms.pres.tmp)/(r_theta.tmp*ms.pres.tmp+(1-r_theta.tmp)*ms.abs.tmp)
  

  
   #run the model!
  prey.est <- find.mle(prey.dat.no0,p.bin=r_theta.tmp,p.full=r_thetap1.tmp,mean.pres=ms.pres.tmp,var.pres=sd.pres.tmp^2,var.abs=sd.abs.tmp^2,beta.mean=betamean.tmp,beta.sd=beta.sd.tmp,p.pop=c_i.tmp)

  #estimate the maximum likelihood estimate for each parameter individually
  source('par.mle.R')
  bnd <- par.mle(prey.dat.no0,prey.est,0.2)
  
  #check whether individual estimates match the estimates from estimating all parameters at once
  if(bnd[1]<bnd[3] & bnd[2]>bnd[3]){
    return(prey.est)
  }
  #else{print("GRRRRR")}
#  } 
}