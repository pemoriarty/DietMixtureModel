#Author: Pamela Moriarty, pmoriart@uw.edu
#Last edited: June 8, 2016
#Purpose: Estimates each parameter individually from the diet mixture model in Moriarty et al. (2016).
#Output: Returns a vector of lower and upper bounds from estimating parameters simultaneously and the c_i estimate from estimating parameters individually


par.mle <- function(prey.dat.no0,prey.est,width){
    r_theta.tmp <- prey.est[1]
    r_thetap1.tmp <- prey.est[2]
    
    notzero.idx <- setdiff(which(prey.dat.no0[,1] > 0), which(prey.dat.no0[,1] ==1))
    beta.mle <- fitdist(prey.dat.no0[notzero.idx],distr="beta",method=c("mle"))[[1]]
      betamean.mle <- beta.mle[1]/sum(beta.mle)
      betasd.mle <- (beta.mle[1]*beta.mle[2])/((beta.mle[1]+beta.mle[2])^2*(beta.mle[1]+beta.mle[2]+1))
      
      ms.abs.mle <- fitdistr(prey.dat.no0[which(prey.dat.no0[,1]==0),2],densfun="gamma")[[1]]  
      ms.abs.mean <- ms.abs.mle[1]/ms.abs.mle[2]
      ms.abs.var <- ms.abs.mle[1]/ms.abs.mle[2]^2
      ms.pres.mle <- fitdistr(prey.dat.no0[which(prey.dat.no0[,1]>0),2],densfun="gamma")[[1]]  
      ms.pres.mean <- ms.pres.mle[1]/ms.pres.mle[2]
      ms.pres.var <- ms.pres.mle[1]/ms.pres.mle[2]^2
      
      c.mle <- r_theta.tmp*((1-r_thetap1.tmp)*betamean.mle*ms.pres.mean+r_thetap1.tmp*ms.pres.mean)/(r_theta.tmp*ms.pres.mean+(1-r_theta.tmp)*ms.abs.mean)

      ubnd<- prey.est[8] + prey.est[8]*width
      lbnd <- prey.est[8] - prey.est[8]*width
      
      return(as.numeric(c(lbnd,ubnd,c.mle)))
}