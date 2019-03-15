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


####################################################################################################
#The bimodal model - marginal likelihood
#Calculate the thermodinamic integral via Parallel Tempering(PT)
#run PT with optimal temperature schedule as in Calderhead and Girolami (2009) (ti=(i/N)^5) with N=30 chains
####################################################################################################

#rm(list=ls())
source("Bimodal_functions_PT.r")
library(MASS)

runPT=function(niter = 35000,M=30){

y=get(load('ST35000_power.RData'))$y
n     = length(y)

# basic MCMC parameters

 # number of iterations for the algorithm to run


#parallel chains temperature parameters



sequence= 1:M
#beta    = 0.05*(sequence/M)+0.95*(sequence/M)^3
beta    = (sequence/M)^5
npar    = length(beta)

# Set up the matrices
# from all tempered chains:
mu      = list()
mllik   = matrix(NA,niter,npar)
#pn_1=pn=matrix(NA,niter,npar)
#pn_1s=pns=matrix(NA,niter,npar)
pn_1s=matrix(NA,niter,npar)



for(chain in 1:npar){
  mu[[chain]]          = rep(NA,niter)
  mu[[chain]][1]    = rnorm(1, mean=mean(y), sd=sd(y))
}

# fill a matrix to keep track of swap proposal / acceptance
swappers = matrix(0,niter,3)


for(iter in 2:niter){

  # potentially select two chains and propose a swap between parameters
  # we don't really want to swap every time, but almost every time would be nice.  
  # I control this using a lazy algorithm:
  maybeswap = sample(0:(npar+1),2,replace=TRUE)
  
if(all(length(unique(maybeswap))==2,min(maybeswap)>0,max(maybeswap)<=npar,iter>2)){
    #propose a parameter swap
    
    C1like1=loglik(y,mean_mu=mu[[maybeswap[1]]][iter-1],tau=beta[maybeswap[1]])$out
    C2like2=loglik(y,mean_mu=mu[[maybeswap[2]]][iter-1],tau=beta[maybeswap[2]])$out

    C2like1=loglik(y,mean_mu=mu[[maybeswap[1]]][iter-1],tau=beta[maybeswap[2]])$out
    C1like2=loglik(y,mean_mu=mu[[maybeswap[2]]][iter-1],tau=beta[maybeswap[1]])$out

  
    if(runif(1)< exp(C1like2+C2like1  -C1like1 - C2like2)){
      #accept the swap
      #exchange parameters information between the swapping chains
      mu[[maybeswap[1]]][iter]      = mu[[maybeswap[2]]][iter-1]
      mu[[maybeswap[2]]][iter]      = mu[[maybeswap[1]]][iter-1]
      
      swappers[iter,]=c(sort(maybeswap),1)
	  
      
    }else{
      # no swap accepted
      # copy the current values of the parameters of interest
      # from the last accepted value

      mu[[maybeswap[1]]][iter]   = mu[[maybeswap[1]]][iter-1]
      mu[[maybeswap[2]]][iter]   = mu[[maybeswap[2]]][iter-1]
      
      swappers[iter,]=c(sort(maybeswap),0)
    }
    
    #Mutation step
    #The other chains still need to move around
    index=1:npar
    index = index[-maybeswap]
    for(chain in index){

     # if (length(mu[[chain]][1:(iter-1)])==1){
     #  post_sd=1 }else{
     #  post_sd=var(mu[[chain]][1:(iter-1)])}


      output               = STstep_pars(th=c(mu[[chain]][iter-1],beta[chain],1),y=y,n=n,post_sd=post_sd)
      mu[[chain]][iter]    = output$theta[1]
    }
    
    
  }else{
    #skip the parameter swap and update as usual
 
  
  for(chain in 1:npar){
      
#       if (length(mu[[chain]][1:(iter-1)])==1){
#       post_sd=1 }else{
#       post_sd=var(mu[[chain]][1:(iter-1)])}

      output              = STstep_pars(th=c(mu[[chain]][iter-1],beta[chain],1),y=y,n=n,post_sd=post_sd)
      mu[[chain]][iter]   = output$theta[1]
  
    }
   }
 
   for(chain in 1:npar){ 

      mllik[iter,chain]=loglik(y,mean_mu=mu[[chain]][iter],tau=1)$mllik 
      pn_1s[iter,chain] = posterior(y,tau=beta[chain],mu=mu[[chain]][iter])
      #print(pn_1s[iter,chain])
  }


 #   if ((iter %% 1000) == 0) {
 #   out_ls=list(mu=mu,mllik=mllik)
 #   save(out_ls,file=paste('PT_ThermodynInt_',k,'_',iter,'.RData',sep=''))
 #  }


}  

out_ls=list(mu=mu,mllik=mllik,pn_1s=pn_1s,swappers=swappers)
return(out_ls)

}