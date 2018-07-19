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
#BNN model
#Calculate the thermodinamic integral via Parallel Tempering(PT)
#run PT with optimal temperature schedule as in Calderhead and Girolami (2009) (ti=(i/N)^5) with N=30 chains
####################################################################################################

#rm(list=ls())


runPT=function(niter = 35000,nChains=30,y,p,M,q1=0.02,q2=0.05, hyperpars=c(rep(10,3),rep(0.05,2)),  x, pickup=NULL){

  


n     = length(y)

# basic MCMC parameters

 # number of iterations for the algorithm to run


#parallel chains temperature parameters



sequence= 1:nChains
#beta    = 0.05*(sequence/M)+0.95*(sequence/M)^3
beta    = (sequence/nChains)^5
npar    = length(beta)

# Set up the matrices
# from all tempered chains:
mu      = list()
mllik   = matrix(NA,niter,npar)

t1=proc.time()[3]
res=optim(c(rnorm(35,0,10),1/rgamma(1,shape=0.05,scale=1/0.05)),function(l) (-1)*loglik(pars=l,y=y,x=x,tau=0.5,p=p,M=M,hyperpars=hyperpars,optim=1),control=list(maxit=500))
t=(as.vector(proc.time())[3]-t1)/60
print(paste('ini; ',t,sep=''))

for(chain in 1:npar){
  mu[[chain]]       = matrix(NA,niter,ncol=36)
  mu[[chain]][1,]   = res$par
}

index=1

if (!is.null(pickup)){
out_ls_PT=get(load(pickup))

index=length(which(!is.na(out_ls_PT$mu[[chain]][,1])))

for(chain in 1:npar){
  
  mu[[chain]][1:index,]   = out_ls_PT$mu[[chain]][1:index,]
  mllik[1:index,]         = out_ls_PT$mllik[1:index,]
}
mem_change(rm(out_ls_PT))
}



# fill a matrix to keep track of swap proposal / acceptance
swappers = matrix(0,niter,3)


for(iter in ((index+1):niter)){

  print(iter)

  timing=matrix(NA,ncol=3,nrow=niter)
  timing[1,]=as.vector(proc.time())[1:3]
  
  # potentially select two chains and propose a swap between parameters
  # we don't really want to swap every time, but almost every time would be nice.  
  # I control this using a lazy algorithm:
  maybeswap = sample(0:(npar+1),2,replace=TRUE)
  
if(all(length(unique(maybeswap))==2,min(maybeswap)>0,max(maybeswap)<=npar,iter>2)){
    #propose a parameter swap
  
  
    C1like1=loglik(y=y,x=x,tau=beta[maybeswap[1]],p=p,M=M,pars=mu[[maybeswap[1]]][iter-1,], hyperpars=hyperpars)
    C2like2=loglik(y=y,x=x,tau=beta[maybeswap[2]],p=p,M=M,pars=mu[[maybeswap[2]]][iter-1,], hyperpars=hyperpars)
  
    C2like1=loglik(y=y,x=x,tau=beta[maybeswap[2]],p=p,M=M,pars=mu[[maybeswap[1]]][iter-1,], hyperpars=hyperpars)	  
    C1like2=loglik(y=y,x=x,tau=beta[maybeswap[1]],p=p,M=M,pars=mu[[maybeswap[2]]][iter-1,], hyperpars=hyperpars)	  
    
  
    # C1like1=loglik(y,mean_mu=mu[[maybeswap[1]]][iter-1],tau=beta[maybeswap[1]])$out
    # C2like2=loglik(y,mean_mu=mu[[maybeswap[2]]][iter-1],tau=beta[maybeswap[2]])$out
    # 
    # C2like1=loglik(y,mean_mu=mu[[maybeswap[1]]][iter-1],tau=beta[maybeswap[2]])$out
    # C1like2=loglik(y,mean_mu=mu[[maybeswap[2]]][iter-1],tau=beta[maybeswap[1]])$out

  
    if(runif(1)< exp(C1like2+C2like1  -C1like1 - C2like2)){
      #accept the swap
      #exchange parameters information between the swapping chains
      mu[[maybeswap[1]]][iter,]      = mu[[maybeswap[2]]][iter-1,]
      mu[[maybeswap[2]]][iter,]      = mu[[maybeswap[1]]][iter-1,]
      
      swappers[iter,]=c(sort(maybeswap),1)
	  
      
    }else{
      # no swap accepted
      # copy the current values of the parameters of interest
      # from the last accepted value

      mu[[maybeswap[1]]][iter,]   = mu[[maybeswap[1]]][iter-1,]
      mu[[maybeswap[2]]][iter,]   = mu[[maybeswap[2]]][iter-1,]
      
      swappers[iter,]=c(sort(maybeswap),0)
    }
    
    #Mutation step
    #The other chains still need to move around
    # index=1:npar
    # index = index[-maybeswap]
    
    t1=proc.time()[3]
   # for(chain in index){
    res = foreach(chain=1:npar) %dopar%{
      
        if (chain %in% maybeswap){
          output =NULL
        }else{
        output               =  SamplePars(y,x=x,theta=mu[[chain]][iter-1,],tau=beta[chain],
                                         hyperpars=hyperpars, q1=q1,q2=q2,
                                         log_r_bot=NULL,p=p,M=M,
                                         accepts=NULL)                                         
        }
      
        return(output )  
        
      
    }
    t=(as.vector(proc.time())[3]-t1)/60
    print(paste('Update pars exchange-mutation; ',t,sep=''))
    
    for(chain in 1:npar){
    if (!is.null(res[[chain]])){
    mu[[chain]][iter,]    = res[[chain]]$theta
    }
    }
    
  }else{
    #skip the parameter swap and update as usual
 
    t1=proc.time()[3]
    # for(chain in index){
    res = foreach(chain=1:npar) %dopar%{
      
  
      output               =  SamplePars(y,x=x,theta=mu[[chain]][iter-1,],tau=beta[chain],
                                         hyperpars=hyperpars, q1=q1,q2=q2,
                                         log_r_bot=NULL,p=p,M=M,
                                         accepts=NULL)                                         
      
      return(output )  
   
    }
    
    t=(as.vector(proc.time())[3]-t1)/60
    print(paste('Update pars mutation; ',t,sep=''))
    
    for(chain in 1:npar){
      mu[[chain]][iter,]    = res[[chain]]$theta
    }
    
   }
 
   for(chain in 1:npar){ 
     mllik[iter,chain]   = mloglik(y=y,x=x,p=p,M=M,pars=mu[[chain]][iter,], hyperpars=hyperpars)	
   }
 
      timing[iter,]=timing[iter-1,]

      if ((iter %% 1000) == 0) {
      timing[iter,]=(as.vector(proc.time())[1:3]-timing[1,])/60  
      out_ls=list(mu=mu,mllik=mllik,timing=timing)
      save(out_ls,file=paste('PT_ThermodynInt_',iter,'.RData',sep=''))
      }


}  

out_ls=list(mu=mu,mllik=mllik,swappers=swappers,timing=timing)
return(out_ls)

}