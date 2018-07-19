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



library(MASS)
library(pryr)
library(doParallel)
source("BNN_PT_fun.R")
source("BNN_PT.R")


data= read.table("gas-furnace.csv",sep = ",",header =  T)
y=data$CO2[1:206]


cluster = makeCluster(detectCores()-1)
registerDoParallel(cluster)

clusterExport(cluster,varlist=list('logpost','loglik','generateInput'))
clusterCall(cluster,function(x) {source('BNN_PT_fun.R');source('BNN_PT.R')})



set.seed(2555)
#niter = 35000
t1=proc.time()[1:3]
y=data$CO2[1:206]
p=4
t1 = as.vector(proc.time())[1:3]
out_ls=runPT(niter = 150000,nChains=15,y=y,p=p,M=5,q1=0.02,q2=0.05, hyperpars=c(rep(10,3),rep(0.05,2)),  x=generateInput(y,p), pickup=NULL)
(t=(as.vector(proc.time())[1:3]-t1)/60)



stopCluster(cluster)


#calculate number of days running
#recent_time="2018-05-04 9:56:00"
#earlier_time="2018-04-07 10:40:00" 

#difftime(recent_time,earlier_time,units="days")
#Time difference of 26.96944 days

