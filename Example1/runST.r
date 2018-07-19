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

source("ST_functions_new.r")


set.seed(2555)
out_ls=runST(niter=50000,SigmaPriorPars=c(1,1),n=25,k=1,switchllik='power',q1=0.15,q2=0.1,kp1=c(1,1),kptau=1)



