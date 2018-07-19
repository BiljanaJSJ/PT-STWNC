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


initialization=function(niter,p,M,nstops,pickup,y,x,kp1,q1,q2,tunetau){
  
  
 
    
  if (is.null(pickup)) {
    
   # theta             = matrix(NA,niter,4)
    theta_1            = matrix(NA,niter,4)
    #some initial parameter values
    #theta[1,]         = c(.0006,.09,257,4)
    tune_q1            = list()
    tune_q1[[1]]       = rep(NA,niter)
    tune_q1[[2]]       = rep(NA,niter)
    tune_q1[[1]][1]    = q1
    tune_q1[[2]][1]    = q1
    
    
    tune_q2           = list()
    tune_q2[[1]]      = rep(NA,niter)
    tune_q2[[2]]      = rep(NA,niter)
    tune_q2[[1]][1]   = q2
    tune_q2[[2]][1]   = q2
    
    
    tau               = rep(NA,niter) 
    tau[1]            = 0.5
    
    tune              = rep(NA,niter) 
    tune[1]           = tunetau
    
    log_r_bot         = list()
    log_r_bot[[1]]    = rep(0,niter)
    log_r_bot[[2]]    = rep(0,niter)
    
    log_r_bot1        = rep(0,niter)
    log_r_top1        = rep(0,niter)
   
    accepts           = list()
    
    
    accepts[[1]]      = matrix(1,nstops,36)
    accepts[[2]]      = matrix(1,nstops,36)
    
    
    
    accepts1          = rep(1,nstops)
    tau_prop          = rep(0,niter)
   
    acc_rate          = list()
    acc_rate[[1]]     = matrix(1,nstops,2)
    acc_rate[[2]]     = matrix(1,nstops,2)
    
    acc_rate1         = rep(1,nstops)
    mllik             = rep(NA,niter)
    failiter          = rep(NA,niter)
    #initialize the values for the first iteration
      
    
    PT_chain           = list()
    PT_chain[[1]]      = matrix(NA,niter,ncol=36)
    PT_chain[[2]]      = matrix(NA,niter,ncol=36)
  
    
    t1=proc.time()[3]
    #cluster = makeCluster(4, type = "SOCK")
    #registerDoParallel(cluster)
    
    #clusterCall(cluster,function(x) {source('BNN_ST.R'); source('BNN_ST_functions.R')})
    
    
    #res = foreach(index=1:2) %dopar%{
      
      
   #   if(index == 1){
        #1. old optimizer
        
   #     return(optim(c(rnorm(35,0,10),1/rgamma(1,shape=0.05,scale=1/0.05)),function(l) (-1)*loglik(pars=l,y=y,x=x,tau=tau[1],p=p,M=M,hyperpars=hyperpars,optim=1),control=list(maxit=500)))
   #   }else if(index==2){
        
   #     return(optim(c(rnorm(35,0,10),1/rgamma(1,shape=0.05,scale=1/0.05)),function(l) (-1)*loglik(pars=l,y=y,x=x,tau=0.99,p=p,M=M,hyperpars=hyperpars,optim=1),control=list(maxit=500)))
   #   }
   # }
    res=optim(c(rnorm(35,0,10),1/rgamma(1,shape=0.05,scale=1/0.05)),function(l) (-1)*loglik(pars=l,y=y,x=x,tau=tau[1],p=p,M=M,hyperpars=hyperpars,optim=1),control=list(maxit=500))
    #stopCluster(cluster)
    
    PT_chain[[1]][1,]  = res$par
    PT_chain[[2]][1,]  = res$par
    
    t=(as.vector(proc.time())[3]-t1)/60
    print(paste('ini:',t,sep=''))   
    
    temp               = matrix(NA,niter,2)
    temp[1,]           = c(tau[1],1)
    
    npar               = 2
    swappers           = matrix(0,niter,3)
    mllik              = list()
    mllik[[1]]         = rep(NA,niter)
    mllik[[2]]         = rep(NA,niter)
    
    kp    = list()
    kp[[1]] = kp1
    kp[[2]] = kp1
   timing=matrix(NA,ncol=3,nrow=niter)
    timing[1,]=as.vector(proc.time())[1:3]
    
    out_ls=list(  #theta             = theta,
      #theta_1           = theta_1,
      PT_chain          = PT_chain,
      tau               = tau,
      timing            = timing,
      log_r_bot         = log_r_bot,
      log_r_top1        = log_r_top1,
      accepts           = accepts,
      accepts1          = accepts1,
      tau_prop          = tau_prop,
      acc_rate          = acc_rate,
      acc_rate1         = acc_rate1,
      tune_q1           = tune_q1,
      tune_q2           = tune_q2,
      mllik             = mllik,
      failiter          = failiter,
      mllik             = mllik,
      temp              = temp,
      tune              = tune,
      tune_q1           = tune_q1,
      tune_q2           = tune_q2,
      kp                = kp,
      swappers          = swappers,
      npar              = npar)
    
  }else{
    
    #if algorithm fails, pick up from where it crushed   
    timing=matrix(NA,ncol=3,nrow=niter)
    
    
    theta_1            = matrix(NA,niter,4)
    #some initial parameter values
    #theta[1,]         = c(.0006,.09,257,4)
    tune_q1            = list()
    tune_q1[[1]]       = rep(NA,niter)
    tune_q1[[2]]       = rep(NA,niter)
    tune_q1[[1]][1]    = q1
    tune_q1[[2]][1]    = q1
    
    
    tune_q2           = list()
    tune_q2[[1]]      = rep(NA,niter)
    tune_q2[[2]]      = rep(NA,niter)
    tune_q2[[1]][1]   = q2
    tune_q2[[2]][1]   = q2
    
    
    tau               = rep(NA,niter) 
    tau[1]            = 0.5
    
    tune              = rep(NA,niter) 
    tune[1]           = tunetau
    
    log_r_bot         = list()
    log_r_bot[[1]]    = rep(0,niter)
    log_r_bot[[2]]    = rep(0,niter)
    
    log_r_bot1        = rep(0,niter)
    log_r_top1        = rep(0,niter)
    
    accepts           = list()
    
    
    accepts[[1]]      = matrix(1,nstops,36)
    accepts[[2]]      = matrix(1,nstops,36)
    
    
    
    accepts1          = rep(1,nstops)
    tau_prop          = rep(0,niter)
    
    acc_rate          = list()
    acc_rate[[1]]     = matrix(1,nstops,2)
    acc_rate[[2]]     = matrix(1,nstops,2)
    
    acc_rate1         = rep(1,nstops)
    mllik             = rep(NA,niter)
    failiter          = rep(NA,niter)
    #initialize the values for the first iteration
    
    
    PT_chain           = list()
    PT_chain[[1]]      = matrix(NA,niter,ncol=36)
    PT_chain[[2]]      = matrix(NA,niter,ncol=36)
    
    
    
    
  #  PT_chain[[1]][1,]  = optim(c(rnorm(35,0,10),1/rgamma(1,shape=0.05,scale=1/0.05)),function(x) (-1)*loglik(pars=x,y=y,tau=tau[1],p=p,M=M,hyperpars=hyperpars,optim=1),control=list(maxit=500))$par
  #  PT_chain[[2]][1,]  = optim(c(rnorm(35,0,10),1/rgamma(1,shape=0.05,scale=1/0.05)),function(x) (-1)*loglik(pars=x,y=y,tau=1,p=p,M=M,hyperpars=hyperpars,optim=1),control=list(maxit=500))$par
    
    
    
    temp               = matrix(NA,niter,2)
    temp[1,]           = c(tau[1],1)
    
    npar               = 2
    swappers           = matrix(0,niter,3)
    mllik              = list()
    mllik[[1]]         = rep(NA,niter)
    mllik[[2]]         = rep(NA,niter)
    
    kp    = list()
    kp[[1]] = kp1
    kp[[2]] = kp1
    
 

    
   #print(out_ls1$PT_chain[[1]][1:10,1:5])
    out_ls1=get(load(pickup))

# newEnv <- function(out_ls1,PT_chain,tau,temp,accepts,accepts1,log_r_bot,tune_q1,tune_q2,kp,tau_prop,log_r_top1,mllik,failiter,swappers, timing,tune ){
#   
#     env <- new.env(parent = globalenv())
#     env$subset <- subset 
#   
#     out_ls <- with(env, {
#   
      
          
    PT_chain[[1]][1:dim(out_ls1$PT_chain[[1]])[1],]=out_ls1$PT_chain[[1]]
    PT_chain[[2]][1:dim(out_ls1$PT_chain[[2]])[1],]=out_ls1$PT_chain[[2]]
    tau[1:length(out_ls1$tau)] = out_ls1$tau
    temp[1:length(out_ls1$tau),] = cbind(out_ls1$tau,1)
    accepts  = out_ls1$accepts 
    accepts1 = out_ls1$accepts1 
    log_r_bot[[1]][1:length(out_ls1$tau)] = out_ls1$log_r_bot[[1]]
    log_r_bot[[2]][1:length(out_ls1$tau)] = out_ls1$log_r_bot[[2]] 
    tune_q1[1:length(out_ls1$tune_q1)] = out_ls1$tune_q1  
    tune_q2[1:length(out_ls1$tune_q2)] = out_ls1$tune_q2  
    kp  = out_ls1$kp    
    tau_prop[1:length(out_ls1$tau_prop )] = out_ls1$tau_prop 
    log_r_top1[1:length(out_ls1$log_r_top1)] = out_ls1$log_r_top1
    mllik[1:length(out_ls1$mllik)]  = out_ls1$mllik    
    failiter[1:length(out_ls1$failiter)] = out_ls1$failiter 
    swappers[1:dim(out_ls1$swappers )[1],] = out_ls1$swappers 
    timing[1:dim(out_ls1$timing )[1],] = out_ls1$timing 
    tune[1:length(out_ls1$tune)] = out_ls1$tune      
    
    
    out_ls=list(PT_chain=PT_chain,tau=tau,accepts=accepts,accepts1=accepts1,log_r_bot=log_r_bot,tune_q1= tune_q1,tune_q2= tune_q2,temp=temp,
                kp=kp,tau_prop=tau_prop,log_r_top1=log_r_top1,mllik=mllik,failiter=failiter,swappers=swappers,tune=tune,timing=timing)
    
      
    # out_ls=list(PT_chain=PT_chain,tau=tau,accepts=accepts,accepts1=accepts1,log_r_bot=log_r_bot,tune_q1= tune_q1,tune_q2= tune_q2,
    #             kp=kp,tau_prop=tau_prop,log_r_top1=log_r_top1,mllik=mllik,failiter=failiter,swappers=swappers,tune=tune,timing=timing)
    # 
#      return(out_ls)
# #     
#    })
# 
#   return(out_ls)
# }
# 
# 
# 
# out_ls=newEnv(PT_chain=PT_chain,tau=tau,temp=temp,accepts=accepts,accepts1=accepts1,log_r_bot=log_r_bot,tune_q1=tune_q1,tune_q2=tune_q2,kp=kp,tau_prop=tau_prop,log_r_top1=log_r_top1,mllik=mllik,failiter=failiter,swappers=swappers, timing=timing,tune=tune)
# 
#print(out_ls$PT_chain[[1]][1:10,1:5])

out_ls$npar                = 2
mem_change(rm(out_ls1))
gc()



  }
  return(out_ls)
}
######################################################################################
#runST runs ST algorithm.
#
# input:       niter            - number of iterations
#              n                - sample size
#              p_mu             - true mean value for the data
#              k                - stdev of the prior of mu (=1), this is adapted
#              sigma            - stdev of the data (=1)
# output:      theta            - bivariate chain with length niter
#              acc_rate         - aceptance rate
#              accepts          - number of accepted proposals
#              log_r_bot        - vector of length niter -(un-normalized log posterior 
#                                 at the last accepted parameter value)
#              kp               - counter for accepted values
#              k_adp            - adapted transition probabilities matrix in order
#                                 to improve acceptance rate
#              y                - data vector
#              tau_profile      - two chains for the profile likelihood of tau
#                                 one obtained from the proposed tau and the other from
#                                 the current iteration (they should be the same)
######################################################################################
runST=function(niter = 100,y,x,p,M,nstops=20,kp1=c(1,1,1),kptau=1,pickup=NULL,
                q1=0.02,q2=0.02,hyperpars=c(rep(10,3),rep(0.05,2)),tunetau){
  
    #timing=matrix(NA,ncol=3,nrow=niter)
    #timing[1,]=as.vector(proc.time())[1:3]
   
    init=initialization(niter=niter,p=p,M=M,nstops=nstops,pickup=pickup,y=y,x=x,kp1=kp1,q1=q1,q2=q2,tunetau=tunetau)
  
    print('init finished')
   
    theta_1           = init$PT_chain[[1]]
    tau               = init$tau
    tau_prop          = init$tau_prop
    log_r_bot         = init$log_r_bot
    log_r_bot1        = init$log_r_bot1
    log_r_top1        = init$log_r_top1
    accepts           = init$accepts
    accepts1          = init$accepts1
    tau_prop          = init$tau_prop
    acc_rate          = init$acc_rate
    acc_rate1         = init$acc_rate1
    mllik             = init$mllik
    failiter          = init$failiter
  
    PT_chain          = init$PT_chain
    temp              = init$temp
    swappers          = init$swappers
    npar              = init$npar
    tune_q1           = init$tune_q1
    tune_q2           = init$tune_q2
    kp                = init$kp
    tune              = init$tune
    timing            = init$timing
    npar              = init$npar


    log_r_bot[[1]]    = logpost(pars=PT_chain[[1]][1,],y=y,x=x,tau=tau[1],p=p,M=M,  hyperpars=hyperpars)
    log_r_bot[[2]]    = log_r_bot[[1]]
    
  
    
    if (!(is.null(pickup))) {
      index            = which(is.na(init$tau))[1]
    }else{
      index            = 2
    }
    
    temp[index,]       = c(tau[index-1],1)
    

          
   out = list()
   out1 = list()
  # library(parallel)
  # cl = makeCluster(getOption("cl.cores", 4))
  
  for(iter in index:niter){
    
    print(iter)
    print(PT_chain[[1]][iter-1,])
    print(PT_chain[[2]][iter-1,])
    print(c(tau[iter-1],tau_prop[iter-1],tune[iter-1]))
    
    maybeswap = sample(0:(npar+1),2,replace=TRUE)
    
    if(all(length(unique(maybeswap))==2,min(maybeswap)>0,max(maybeswap)<=npar,iter>2)){
      #propose a parameter swap
      t1=proc.time()[3] 
   
	  l11    = loglik(y=y,x=x,tau=temp[iter-1,maybeswap[1]],p=p,M=M,pars=PT_chain[[maybeswap[1]]][iter-1,], hyperpars=hyperpars)	  
	  l22    = loglik(y=y,x=x,tau=temp[iter-1,maybeswap[2]],p=p,M=M,pars=PT_chain[[maybeswap[2]]][iter-1,], hyperpars=hyperpars)	  
		  
	  l21    = loglik(y=y,x=x,tau=temp[iter-1,maybeswap[2]],p=p,M=M,pars=PT_chain[[maybeswap[1]]][iter-1,], hyperpars=hyperpars)	  
	  l12    = loglik(y=y,x=x,tau=temp[iter-1,maybeswap[1]],p=p,M=M,pars=PT_chain[[maybeswap[2]]][iter-1,], hyperpars=hyperpars)	  

	  t=(as.vector(proc.time())[3]-t1)/60
	  print(paste('swap chains',t,sep=''))        
     
    if (runif(1)< exp(l12 + l21 -l11 - l22)){
        #accept the swap
        
        PT_chain[[maybeswap[1]]][iter,] =PT_chain[[maybeswap[2]]][iter-1,]
        PT_chain[[maybeswap[2]]][iter,] =PT_chain[[maybeswap[1]]][iter-1,]
        
        swappers[iter,]=c(maybeswap,1)
        
      }else{
        # no swap accepted
        PT_chain[[maybeswap[1]]][iter,] =PT_chain[[maybeswap[1]]][iter-1,]
        PT_chain[[maybeswap[2]]][iter,] =PT_chain[[maybeswap[2]]][iter-1,]
        
        
        swappers[iter,]=c(maybeswap,0)
      }
      
      tau[iter]   = tau[iter-1]
      temp[iter,] = c(tau[iter-1],1)
    
      
    }else{
    
    chain=1
    theta_1[iter,]            =  PT_chain[[chain]][iter-1,]
    
    t1=proc.time()[3]
    
    out                       =  SamplePars(y,x=x,theta=PT_chain[[chain]][iter-1,],tau=tau[iter-1],
                                            hyperpars=hyperpars,
                                            log_r_bot=log_r_bot[[chain]][iter-1],p=p,M=M,
                                            accepts=accepts[[chain]][kp[[chain]],],
                                            q1=tune_q1[[chain]][iter-1],q2=tune_q2[[chain]][iter-1],
                                            adaptvars=apply(PT_chain[[chain]][,1:35],2,var,na.rm=T ))
    
    t=(as.vector(proc.time())[3]-t1)/60
    print(paste('Update pars chain 1;',t,sep=''))  
 
    PT_chain[[chain]][iter,]                 = out$theta
    for (i in (1:36)){
    accepts[[chain]][kp[[chain]][1],i]       = out$accepts[i]
    # accepts[[chain]][kp[[chain]][2],2]       = out$accepts[2]
    # accepts[[chain]][kp[[chain]][3],3]       = out$accepts[3]
    }
    log_r_bot[[chain]][iter]                 = out$log_r_bot
    
	
	
    t1=proc.time()[3]
    
    out1                = tryCatch({   
	
	

      
                                    STstep_tau(te=tau[iter-1],theta=PT_chain[[chain]][iter,],y=y,x=x,p=p,M=M, hyperpars=hyperpars,acc=accepts1[kp[[1]]],tune=tune[iter-1] ,iter=iter,samples=cbind(PT_chain[[chain]],tau))  
                                     
                                    },
                                    error=function(e){
                                       message(paste("er:",e,sep=''))
                                       return(NA)
                                       },                                
                                    
                                    warning=function(wr) {
                                       message(paste("wr:",wr,sep=''))
                                       return(NULL)
                                    }
                                    )   
   
  t=(as.vector(proc.time())[3]-t1)/60
  print(paste('Update tau; ',t,sep=''))
    

 if (is.list(out1)) {
      
       theta_1[iter,]             = out$theta
 
       tau[iter]                  = out1$te
       tau_prop[iter]             = out1$tau_prop
       temp[iter,]                = c(tau[iter],1)
       accepts1[kp[[1]]]          = out1$accepts1
       log_r_bot1[iter]           = out1$log_r_bot
       tau_prop[iter]             = out1$tau_prop
       log_r_top1[iter]           = out1$log_r_top1
  
     
      
   
     }else{
      
       mllik[[chain]][iter]                = mllik[[chain]][iter-1] 
       tau[iter]                           = tau[iter-1]
       tau_prop[iter]                      = tau_prop[iter-1]
       temp[iter,]                         = c(tau[iter],1)
       PT_chain[[chain]][iter,]            = theta_1[iter,]
   
       failiter[iter]                      = iter
                      
    }

   
    chain=2
  
    t1=proc.time()[3]
    
    out                       =  SamplePars(y,x=x,theta=PT_chain[[chain]][iter-1,],tau=tau[iter-1],
                                            hyperpars=hyperpars,
                                            log_r_bot=log_r_bot[[chain]][iter-1],p=p,M=M,
                                            accepts=accepts[[chain]][kp[[chain]],],
                                            q1=tune_q1[[chain]][iter-1],q2=tune_q2[[chain]][iter-1],
                                            adaptvars=apply(PT_chain[[chain]][,1:35],2,var,na.rm=T ) )
  
    t=(as.vector(proc.time())[3]-t1)/60
    print(paste('Update pars chain 2;',t,sep=''))  
    
    #  print(out$theta)
    PT_chain[[chain]][iter,]                 = out$theta
    accepts[[chain]][kp[[chain]][1],1]       = out$accepts[1]
    accepts[[chain]][kp[[chain]][2],2]       = out$accepts[2]
    accepts[[chain]][kp[[chain]][3],3]       = out$accepts[3]
    log_r_bot[[chain]][iter]                 = out$log_r_bot
  
    }
  
  chain=1
  
  mllik[[chain]][iter]                 = mloglik(y=y,x=x,p=p,M=M,pars=PT_chain[[chain]][iter,], hyperpars=hyperpars)	
  
  chain=2
  
  mllik[[chain]][iter]                 = mloglik(y=y,x=x,p=p,M=M,pars=PT_chain[[chain]][iter,], hyperpars=hyperpars)	  
  
  
  tune[iter]                           = tune[iter-1]  
  

for (chain in (1:npar)){
  
  tune_q1[[chain]][iter]=tune_q1[[chain]][iter-1]
  tune_q2[[chain]][iter]=tune_q2[[chain]][iter-1]
  
  
}
  
  # for (chain in (1:npar)){
  # if(iter==kp*niter/nstops){  
  #   # regularly save and make updates
  #   rate[[chain]][,kp[[chain]]] = accepts[[chain]][,kp[[chain]]]/(niter/nstops);
  #   kp[[chain]]       = kp[[chain]]+1;
  #   #save.image(filename)
  #   
  #   if((rate[[chain]][1,kp[[chain]]-1]>.49 ||  rate[[chain]][1,kp[[chain]]-1]<.39)   && iter< niter/2){
  #     # adjust the transition density if we are in the first half of the iterations. and the acceptance rate is too high or too low.
  #     stepvar[1,1] = stepvar[1,1]*rate[1,kp-1]/.44; 
  #   }
  #   
  #   if((rate[2,kp-1]>.49 ||  rate[2,kp-1]<.39)   && iter< niter/2){
  #     # adjust the transition density if we are in the first half of the iterations. and the acceptance rate is too high or too low.
  #     stepvar[2,1] = stepvar[2,1]*rate[2,kp-1]/.44; 
  #   }
  #   
  # }
  # }
  # 
  
    timing[iter,]=timing[iter-1,]
   
    if (((iter %% 1000) == 0) | iter==7001){
      
      
      timing[iter,]=(as.vector(proc.time())[1:3]-timing[1,])/60  
      
      out_ls=list(PT_chain=PT_chain,tau=tau,accepts=accepts,accepts1=accepts1,log_r_bot=log_r_bot,tune_q1= tune_q1,tune_q2= tune_q2,
                  kp=kp,tau_prop=tau_prop,log_r_top1=log_r_top1,mllik=mllik,failiter=failiter,swappers=swappers,tune=tune,timing=timing)
      
     
      
      save(out_ls, file=paste("ST",iter,".RData",sep=""))
      save(".Random.seed",file=paste("random_state_",iter,".RData",sep=''))
    }
  }
  
#stopCluster(cl)
  
  return(list(PT_chain=PT_chain,tau=tau,accepts=accepts,accepts1=accepts1,log_r_bot=log_r_bot,tune_q1= tune_q1,tune_q2= tune_q2,
              kp=kp,tau_prop=tau_prop,log_r_top1=log_r_top1,mllik=mllik,failiter=failiter,swappers=swappers,tune=tune,timing=timing))
}
######################################################################################
