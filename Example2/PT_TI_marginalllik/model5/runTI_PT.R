
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

#run PT-STWNC 20 times to obtain marginal log-likelihood

source('galaxy_functions_PT.R')
source('Galaxy_Mixture_parallel_Temper.R')

library(MASS)
library(MCMCpack)
library(truncnorm)

M=30
sequence=1:M
beta    = (sequence/M)^5

k1=5
for (i in (1:10)){
out_ls=runTI_PT(niter=35000,data=galaxies,k=k1,Meanpriorpars= cbind(rep(20,k1),rep(100,k1)),
                  Sigmapriorpars=c(3,20),
                  Ppriorpars= rep(1,k1),beta=beta,equalVar=NULL,eval_name=i)

}




