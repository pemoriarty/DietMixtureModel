par.mle <- function(prey.dat.no0,prey.est,width){
    p.bin.tmp <- prey.est[1]
    p.full.tmp <- prey.est[2]
    
    notzero.idx <- setdiff(which(prey.dat.no0[,1] > 0), which(prey.dat.no0[,1] ==1))
    beta.mlepar <- fitdist(prey.dat.no0[notzero.idx],distr="beta",method=c("mle"))[[1]]
      beta.mle <- beta.mlepar[1]/sum(beta.mlepar)
      beta.mlesd <- (beta.mlepar[1]*beta.mlepar[2])/((beta.mlepar[1]+beta.mlepar[2])^2*(beta.mlepar[1]+beta.mlepar[2]+1))
      
      msa.mle <- fitdistr(prey.dat.no0[which(prey.dat.no0[,1]==0),2],densfun="gamma")[[1]]  
      msa.mean <- msa.mle[1]/msa.mle[2]
      msa.var <- msa.mle[1]/msa.mle[2]^2
      mso.mle <- fitdistr(prey.dat.no0[which(prey.dat.no0[,1]>0),2],densfun="gamma")[[1]]  
      mso.mean <- mso.mle[1]/mso.mle[2]
      mso.var <- mso.mle[1]/mso.mle[2]^2
      
      c.mle <- p.bin.tmp*((1-p.full.tmp)*beta.mle*mso.mean+p.full.tmp*mso.mean)/(p.bin.tmp*mso.mean+(1-p.bin.tmp)*msa.mean)

      ubnd<- prey.est[8] + prey.est[8]*width
      lbnd <- prey.est[8] - prey.est[8]*width
      
      return(as.numeric(c(lbnd,ubnd,c.mle)))
}