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


source('SIR_ST.R')
source('SIR_ST_functions.R')
#source("SIR-model-set-up_no_fixed_X0.R")
data              = read.table("Eyam_time_SIR.csv",sep = ",")
data[data=="NaN"] = NA
data              = as.matrix(data)
colnames(data)    = c("time","S","I","R")
times             = data[,"time"]; 
data              = data[,-1]
y                 = data[,3]
N=261

set.seed(2555)

t1=proc.time()[1]
out_ls=runST(niter = 35000,y,nstops=20,kp1=c(1,1,1),kptau=1,pickup=NULL,
             x=times,N=N,q1=0.02,q2=0.02,priorpars=c(1,1,N,5/N),tunetau=0.1)
(t=proc.time()[1]-t1)/60



