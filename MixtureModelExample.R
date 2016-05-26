#author: Pamela Moriarty
#last edited: May 11, 2016
#Purpose: We demonstrate the use of the mixture model described in:
#Moriarty, P.E., T.E. Essington, and E.J. Ward. 2016. A Novel Method to Estimate Prey Contributions to Predator Diets. Journal of Blah. Submitted. 

#####install and load necessary packages#####
install.packages('plotrix')
install.packages('stats4')
install.packages('optimx')
install.packages('numDeriv')
install.packages('MASS')
install.packages('fitdistrplus')

require(plotrix)
require(stats4)
require(optimx)
require(numDeriv)
require(MASS)
require(fitdistrplus)

#####load example data#####
load('PredatorExample.Rdata')#load example datasets
#sim.data is a simulated dataset
#ling.dat is the lingcod data used in Moriarty et al. (in review)


#####analyze simulated data#####
source('run.model.R')#prints the parameter estimates for all 8 mixture model parameters, 

source('model.comparison.R')#estimates the prey contribution, c_i, using the mixture model, weighted mean and mean and calculates  error for each estimate

run.model(sim.data[,1],sim.data[,2])#returns all parameter estimates from the maximum likelihood estimation procedure
#parameter= true values of parameters
# r_theta = 0.8
# r_thetap=1 = 0.5
# m_s|r_theta = 11
# var_{rtheta} = 4
# var_{not rtheta} = 2
# betamean = 0.7
# betasd = 0.3
# c_i = 0.8

model.comparison(sim.data[,1],sim.data[,2],mat=T)#compare estimate and standard error from our mixture model to a conventional mean and weighted mean

#####analyse lingcod data#####
#ling.data: column1 = total stomach contents mass; column2=proportion of most common prey type; column3=proportion of prey type with covariance

#compare estimates between the mixture model, mean and weighted mean
model.comparison(ling.data[,1],ling.data[,3],mat=T)#common prey type
model.comparison(ling.data[,2],ling.data[,3],mat=T,CI=T)#prey type with covariance

#get estimates for all model parameters
run.model(ling.data[,1],ling.data[,3])
run.model(ling.data[,2],ling.data[,3])
