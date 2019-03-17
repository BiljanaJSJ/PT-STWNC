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



source('../../Functions/ST_functions.R')
source("../ST_pants.r")
library(truncnorm)
library("MASS")
library(coda)
library(MCMCpack)



mllik=rep(NA,20)
mean_se=rep(NA,20)
n=25
low=1000
up=35000
n=25
p_mu=1.5
sigma=1
set.seed(2555)
#set.seed(1555)

y=data(n=n,mu=abs(p_mu),sigma=sigma)


for (i in (1:20)){
##out_ls=get(load(paste('mllik',i,'.RData',sep='')))

out_ls=runST(niter=35000,y=y,PriorPars=c(0,1,1,1),IniPar=c(-0.8,1,0.5),tune_pars_init=c(0.15,0.1),ttune_pars=c(rep(TRUE,2)))

save(out_ls,file=paste('mllik',i,'.RData',sep=''))


tau_samples=c(out_ls$PT_chain[[1]][low:up,3],rep(1,length(out_ls$PT_chain[[1]][low:up,3])))

mloglik_samples=out_ls$mllik[[1]][low:up]



indxSortedtau=sort(tau_samples,index=T)

SortedThetaTau=tau_samples[indxSortedtau$ix]
uniqueTau=unique(SortedThetaTau)
tau_int=diff(uniqueTau)



indicesUniqueTau=lapply(1:length(uniqueTau),function(x){ which(tau_samples==uniqueTau[x]) })

ThermodynIntergal=lapply(1:length(uniqueTau),function(x){ mean(mloglik_samples[indicesUniqueTau[[x]]],na.rm=T)*tau_int[x] })


unThermDynInt=unlist(ThermodynIntergal)
unThermDynInt=unThermDynInt[which(!is.na(unThermDynInt))]
mllik[i]=sum(unThermDynInt)

}

saveall=list(mllik,sd=sd(mllik))
save(mllik, file='savels_AllChains.RData')



#combine all the samples together




#solve the integral analitically


I=log(((2*pi)^(-n/2))*sqrt((n+1)^(-1))*exp( (1/2)*((sum(y)^2)/(n+1)) - (1/2)*(sum(y^2))    ))

save(I,file='analiticalSol.RData')
