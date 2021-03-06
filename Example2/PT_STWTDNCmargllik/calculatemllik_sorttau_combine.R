
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

source('../../ST_Galaxy_functions.R')
source('../../../Functions/ST_functions.R')

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

k1=2
data               = galaxies
data[78]           = 26960 ## there's a typo in the dataset
y                  = data/1000 ## normalize the data 
n                  = length(y)


fmpij      = mpij(y=y,means=rep(20,k1),sigma2=0.5,p=rep(1/k1,k1))
pij        = exp(fmpij)/apply(exp(fmpij),1,sum)
Z          = sampleZ(pij)


for (i in (1:10)){
  
 # out_ls=get(load(paste('mllik',i,'.RData',sep='')))
 
 out_ls=runST(niter =35000,y=y,
 	     PriorPars= c(20,100,3,20,rep(1,k1)),
	     IniPar=c(rep(20,k1),0.5,rep(1/k1,k1),k1,0.5),
	     tune_pars_init =c(rep(1,k1),0.5,rep(1,k1),1),
	     ttune_pars=c(rep(FALSE,k1),TRUE,rep(FALSE,k1),FALSE),
	     parAdd=Z,ttau=c(6,TRUE))



 
  save( out_ls, file=paste('mllik',i,'.RData',sep=''))

  tau_samples=c(out_ls$PT_chain[[1]][low:up,ncol(out_ls$PT_chain[[1]])],rep(1,length(out_ls$PT_chain[[1]][low:up,ncol(out_ls$PT_chain[[1]])])))
  mloglik_samples=c(out_ls$mllik[[1]][low:up],out_ls$mllik[[2]][low:up])


  indxSortedtau=sort(tau_samples,index=T)
  
  SortedThetaTau= tau_samples[indxSortedtau$ix]
  
  
  uniqueTau=c(0,unique(SortedThetaTau))
  tau_int=diff(uniqueTau)
  
 
  
  indicesUniqueTau=lapply(1:length(uniqueTau),function(x){ which(tau_samples==uniqueTau[x]) })
  
  ThermodynIntergal=lapply(1:length(uniqueTau),function(x){ mean(mloglik_samples[indicesUniqueTau[[x]]],na.rm=T)*tau_int[x] })
  
 
  unThermDynInt=unlist(ThermodynIntergal)
  unThermDynInt=unThermDynInt[which(!is.na(unThermDynInt))]
  mllik[i]=sum(unThermDynInt)
  
}

Savels=list(mllik=mllik,se=sd(mllik))
save(Savels,file='savelist_allChains.RData')




