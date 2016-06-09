#Author: Pamela Moriarty, pmoriart@uw.edu
#Last edited: June 8, 2016
#Purpose: Calls the diet mixture model described in Moriarty et al. (2016)
#Output: Returns a vector of parameter estimates


run.model <- function(prop,total.mass){
  source('model.par.R')

  out <- model.par(prop,total.mass)
  if(class(out)[1]=="optimx"){
  names(out[1:8])<- c("r_theta","r_thetap=1","ms_pres","var_ms_pres","var_ms_abs","beta_mean","beta_sd","c")
  if(out[13]==0){
    est <- c(out[1:8])
    names(est) <- c("r_{theta}","r_{thetap=1}","m_s|r_theta","var_{rtheta}","var_{not rtheta}","betamean","betasd","c_i")
    
    return(as.matrix(est))
  }else{
    return(print("Model didn't converge"))
  }
  }else{
    print(out)
  }
}