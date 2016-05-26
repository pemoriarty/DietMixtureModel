run.model <- function(prop,total.mass){
  source('model.par.R')
  dens.warn <- function(w) {
    if(any(grepl("densfun",w))){
      invokeRestart("muffleWarning")}
  }
  out <- withCallingHandlers(model.par(prop,total.mass),warning=dens.warn)
  names(out[1:8])<- c("r_theta","r_thetap=1","ms_pres","var_ms_pres","var_ms_abs","beta_mean","beta_sd","c")
  if(out[13]==0){
    est <- c(out[1:8])
    names(est) <- c("r_theta","r_thetap=1","ms_|r_theta","var_{rtheta}","var_{not rtheta}","betamean","betasd","c_i")
    
    return(as.matrix(est))
  }else{
    return(print("Model didn't converge"))
  }
}