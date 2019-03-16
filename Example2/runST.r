####################################################################
#article title:  Parallel Tempering via Simulated Tempering Without 
#                Normalizing Constants
#journal name:   Statistics and Computing
#author names:   Biljana Jonoska Stojkova, David A. Campbell
#affiliation 
#and e-mail 
#address of the 
#corresponding 
#author:         Department of Statistics 
#                University Of British Columbia
#                b.stojkova@stat.ubc.ca
####################################################################

source('../Functions/ST_functions.R')
source('ST_Galaxy_functions.R')

library(MASS)
library(MCMCpack)
library(truncnorm)



t1=proc.time()[1]

k1=3
data               = galaxies
data[78]           = 26960 ## there's a typo in the dataset
y                  = data/1000 ## normalize the data 
n                  = length(y)


fmpij      = mpij(y=y,means=rep(20,k1),sigma2=rep(0.5,k1),p=rep(1/3,k1))
pij        = exp(fmpij)/apply(exp(fmpij),1,sum)
Z          = sampleZ(pij)


out_ls=runST(niter =35000,y=y,
	     PriorPars= c(20,100,3,20,rep(1,k1)),
	     IniPar=c(rep(20,k1),rep(0.5,k1),rep(1/3,k1),k1,0.5),
	     tune_pars_init =c(rep(1,k1),rep(0.05,k1),rep(1,k1),1),
	     ttune_pars=c(rep(FALSE,k1),rep(TRUE,k1),rep(FALSE,k1),FALSE),
	     parAdd=Z,ttau=c(2,TRUE))

(t=proc.time()[1]-t1)/60
