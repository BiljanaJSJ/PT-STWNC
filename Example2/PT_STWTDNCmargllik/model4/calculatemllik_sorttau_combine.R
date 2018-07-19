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


#sort the samples [theta, tau | Y] by tau using N=up-low samples
#get M<N unique tau slices 
#sum_i_{1_M} sum_j_{1_ni} log(P(Y | theta_i_j)) [ti+1-ti]

source('ST_Galaxy_functions.R')
source('ST_Galaxy_runST.R')

library(MASS)
library(MCMCpack)
library(truncnorm)


mllik=rep(NA,10)
mean_se=rep(NA,10)


mllik=rep(NA,10)
mean_se=rep(NA,10)
#n=25
low=1000
up=35000

for (i in (1:10)){
  
  out_ls=get(load(paste('mllik',i,'.RData',sep='')))

  tau_samples=c(out_ls$tau[low:up],rep(1,length(out_ls$tau[low:up])))
  mloglik_samples=c(out_ls$mllik[[1]][low:up],out_ls$mllik[[2]][low:up])


  indxSortedtau=sort(tau_samples,index=T)
  
  SortedThetaTau= tau_samples[indxSortedtau$ix]
  #SortedThetaTau[which(is.na(SortedThetaTau))]=1
  
  uniqueTau=c(0,unique(SortedThetaTau))
  tau_int=diff(uniqueTau)
  
  #ind=which(out_ls$PT_chain[[1]][,2]==uniqueTau[1])
  #out_ls$mllik[[1]][ind]*tau_int[1]
  
  indicesUniqueTau=lapply(1:length(uniqueTau),function(x){ which(tau_samples==uniqueTau[x]) })
  
  ThermodynIntergal=lapply(1:length(uniqueTau),function(x){ mean(mloglik_samples[indicesUniqueTau[[x]]],na.rm=T)*tau_int[x] })
  
 # ThermodynIntergal=lapply(2:length(uniqueTau),function(x){ (1/2)*(mean(mloglik_samples[indicesUniqueTau[[x]]],na.rm=T)+mean(mloglik_samples[indicesUniqueTau[[x-1]]],na.rm=T))*tau_int[x] })
  
  unThermDynInt=unlist(ThermodynIntergal)
  unThermDynInt=unThermDynInt[which(!is.na(unThermDynInt))]
  mllik[i]=sum(unThermDynInt)
  
}

Savels=list(mllik=mllik,se=sd(mllik))
save(Savels,file='savelist_allChains.RData')







#combine all the samples together


tau_samples=list()
mloglik_samples=list()

for (i in (1:10)){
out_ls=get(load(paste('mllik',i,'.RData',sep='')))
tau_samples[[i]]=c(out_ls$tau[low:up],rep(1,length(out_ls$tau[low:up])))
mloglik_samples[[i]]=c(out_ls$mllik[[1]][low:up],out_ls$mllik[[2]][low:up])
}

tau_samples       = unlist(tau_samples)
mloglik_samples   = unlist(mloglik_samples)

indxSortedtau     = sort(tau_samples,index=T)

SortedThetaTau    = tau_samples[indxSortedtau$ix]

uniqueTau         = c(0,unique(SortedThetaTau))
tau_int           = diff(uniqueTau)




library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
clusterCall(cl,function(x) {source('ST_Galaxy_functions.R');source('ST_Galaxy_runST.R')})
clusterExport(cl,varlist=ls())

indicesUniqueTau  = parLapply(cl,1:length(uniqueTau),function(x){ which(tau_samples==uniqueTau[x]) })
stopCluster(cl)

ThermodynIntergal = lapply(1:length(uniqueTau),function(x){ mean(mloglik_samples[indicesUniqueTau[[x]]],na.rm=T)*tau_int[x] })

#ThermodynIntergal=lapply(2:length(uniqueTau),function(x){ (1/2)*(mean(mloglik_samples[indicesUniqueTau[[x]]],na.rm=T)+mean(mloglik_samples[indicesUniqueTau[[x-1]]],na.rm=T))*tau_int[x] })

unThermDynInt     = unlist(ThermodynIntergal)
unThermDynInt     = unThermDynInt[which(!is.na(unThermDynInt))]
mllik             = sum(unThermDynInt)

save(mllik,file='AllSamplesCombined_allChains.RData')
