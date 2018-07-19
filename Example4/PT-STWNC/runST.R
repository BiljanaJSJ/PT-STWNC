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
 
source('BNN_ST.R')
source('BNN_ST_functions.R')
library(doParallel)
data         = read.table("gas-furnace.csv",sep = ",",header =  T)
library(doSNOW)
cluster = makeCluster(detectCores()-1)
registerDoParallel(cluster)
 
clusterExport(cluster,list=list('logpost','loglik','generateInput'))
clusterCall(cluster,function(x) {source('BNN_ST.R');source('BNN_ST_functions.R')})
set.seed(2555)
#niter = 35000
t1=proc.time()[1:3]
y=data$CO2[1:206]
p=4
out_ls=runST(niter = 455000,y=y,p=p,M=5,nstops=20,kp1=c(1,1,1),kptau=1,
              q1=0.02,q2=0.05, hyperpars=c(rep(10,3),rep(0.05,2)),tunetau=0.1,  x=generateInput(y,p))


stopCluster(cluster)