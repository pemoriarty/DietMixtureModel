#Author: Pamela Moriarty, pmoriart@uw.edu
#Last edited: June 8, 2016
#Purpose: calculate the expectation of a beta distribution

Bmean <- function(alpha1,alpha2){
  #alpha1, alpha2- parameters of the beta distribution
  return(alpha1/(alpha1 + alpha2))
}