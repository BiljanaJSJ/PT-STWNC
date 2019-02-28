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


source('../Functions/ST_functions.r')
source('SIR_ST_functions.R')
library(parallel)

data              = read.table("Eyam_time_SIR.csv",sep = ",")
data[data=="NaN"] = NA
data              = as.matrix(data)
colnames(data)    = c("time","S","I","R")
times             = data[,"time"]; 
data              = data[,-1]
y                 = data[,3]
N                 = 261

set.seed(2555)

t1=proc.time()[1]

cl = makeCluster(getOption("cl.cores", 4))

out_ls=runST(niter = 35000,y=y,kp1=1,kptau=1,
	     IniPar=c(.0006,.09,257,4,0.5),
	     tune_pars_init=c(0.02,0.02,5,N-5),
	     ttune_pars=c(rep(FALSE,4)),
             PriorPars=c(1,1,N,5/N),ttau=c(0.1,FALSE),
	     parAdd=list(times=times,N=N), cl=cl)

stopCluster(cl)
(t=proc.time()[1]-t1)/60



