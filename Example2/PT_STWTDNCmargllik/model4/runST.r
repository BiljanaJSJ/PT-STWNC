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


source('ST_Galaxy_functions.R')
source('ST_Galaxy_runST.R')

library(MASS)
library(MCMCpack)
library(truncnorm)

set.seed(2555)
t1=proc.time()[1]
k1=3
out_ls=runST(niter=10000,data=galaxies,k=k1,Meanpriorpars= cbind(rep(20,k1),rep(100,k1)),
             Sigmapriorpars=c(3,20),
             Ppriorpars= rep(1,k1),nstops=20,adps_k=0.05,
             kptau=1,kp1=c(1,1),adpm_k= 2,equalVar=1)
(t=proc.time()[1]-t1)/60
