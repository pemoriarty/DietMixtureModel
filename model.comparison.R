#Pamela Moriarty
#last edited August 19, 2015
#from data, calculates prey contribution using the mean, weighted mean and GULPS model
#returns either a plot or table with estimates

model.comparison <- function(prop,total.mass,mat=F, CI=F, yaxis=T, ybnd=1.0,...){
  #prop- vector of the proportion of each stomach made up of the prey type
  #total.mass- vector of the total mass of contents in the stomachs
  #mat- if true, returns a matrix of estimates and error, otherwise returns a plot
  #CI- if true, returns 95% confidence intervals for error, otherwise returns standard error
  #yaxis- if plot is returned, should the y axis be displayed
  #ybnd- if plot is returned, what is the upper limit for the y axis on the plot
  #...- other graphical parameters that will be passed to the plotting function
  require(numDeriv)
  require(MASS)
  require(fitdistrplus)
  source('model.par.R')
  source('NLL.mle.R')
  require(plotrix)
  dens.warn <- function(w) {
    if(any(grepl("densfun",w))){
      invokeRestart("muffleWarning")}
  }
    mine <- withCallingHandlers(model.par(prop,total.mass),warning=dens.warn)
    if(!is.list(mine)){return(mine)}
  
    #calculate p.pop estimates from the 3 models
    #mine.est <- mine$par[8]
    mine.est <- as.numeric(mine[8])
    mean.est <- mean(prop,na.rm=T) 
    weighted.est <- sum(prop*total.mass,na.rm=T)/sum(total.mass,na.rm=T)
  
  prop <- replace(prop,is.nan(prop),0)#remove NaNs from the proportion data
  prop <- replace(prop,is.na(prop),0)#remove NAs from the proportion data
  total.mass <- replace(total.mass,is.nan(total.mass),0)#remove NaNs from the total stomach mass data
  total.mass <- replace(total.mass,is.na(total.mass),0)#remova NAs from the total stomach mass data
  prey.dat <- cbind(prop,total.mass) 
  prey.dat[,1] <- replace(prey.dat[,1],is.nan(prey.dat[,1]),0)#remove NaNs from the data
  prey.dat[,1] <- replace(prey.dat[,1],is.na(prey.dat[,1]),0)#remove NAs from the data
  prey.dat[,1] <- replace(prey.dat[,1],is.infinite(prey.dat[,1]),0)#remove INFs from the data
  prey.dat.no0 <- subset(prey.dat,prey.dat[,2] > 0)#remove empty stomachs 
  prey.dat.no0 <- as.matrix(prey.dat.no0)
  data <- prey.dat.no0
  
  #calculate standard error for the 3 models
  hess <- hessian(NLL.mle,c(as.numeric(mine[1:8])),method="Richardson",data=prey.dat.no0)
  mine.se <- try(sqrt(diag(solve(hess))[8]))
    if(class(mine.se)=='try-error'){
      mine.se <- try(sqrt(diag(solve(hess[-2,-2]))[7]))
    }
  #mine.se <- sqrt((diag(solve(mine$hessian)))[8])
  var.se <- (sd(prop,na.rm=T))/sqrt(length(prop))
  weighted.se <- sqrt((length(prop)/((length(prop)-1)*(sum(total.mass)^2)))*(sum((total.mass*prop-mean(total.mass)*weighted.est)^2)-2*weighted.est*sum((total.mass-mean(total.mass))*(total.mass*prop-mean(total.mass)*weighted.est))+weighted.est^2*sum((total.mass-mean(total.mass))^2)))#equation from Gatz and Cochran (1995)
  
  if(CI==TRUE){#calculate confidence intervals, if desired
    mine.var <- 1.96*mine.se
    var.est <- 1.96*var.se
    weighted.var <- 1.96*weighted.se
    error.type = "95% CI"
  }else{
    if(class(mine.se)!='try-error') {mine.var <- mine.se} else{mine.var <- NA}
    var.est <- var.se
    weighted.var <- weighted.se
    error.type = "SE"
  }
  
  if(mat==T){#if matrix output is desired, create the matrix
    est.mat <- data.frame(Mixture=c(mine.est, mine.var), Mean=c(mean.est,var.est), Weighted=c(weighted.est,weighted.var),row.names=c("Estimate",error.type))
  return(est.mat)
  }else{#otherwise, return a plot
  
    par(mar=c(6,6,1,5))#create the plot
    plot(c(mine.est,mean.est,weighted.est) ~ c(1,1.25,1.5),pch=16,cex=4,ylim=c(0,ybnd),xaxt="n",xlab="",ylab="",yaxt='n',col=c("black","blue","red"),...)
    plotCI(x=c(1,1.25,1.5),y=c(mine.est,mean.est,weighted.est),uiw=c(mine.var,var.est,weighted.var),liw=c(mine.var,var.est,weighted.var),err="y",add=T,lwd=4,col=c("black","blue","red"))
    mtext(side=1,at=c(1,1.25,1.5),text=c("Mixture","Mean","Weighted"),cex=1.8,col=c("black","blue","red"),line=1.5)
    mtext(side=1,at=c(1,1.25,1.5),text=c("Model","","Mean"),cex=1.8,col=c("black","blue","red"),line=2.75)
    axis(side=1,at=c(1,1.25,1.5),labels=F)
    mtext(side=1,text="Method",line=4.5,cex=2)
  if(yaxis==T){
    mtext(side=2,text="Prey Contribution to Predator's Diet",line=4,cex=2)
    axis(side=2,at=c(0,ybnd/2,ybnd),las=1,cex.axis=1.8)
  }
    mtext(side=1,text="Error bars = Standard Error",line=4.2,at=1.44,cex=1.5)
  }
}