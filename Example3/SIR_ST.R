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


initialization=function(niter,nstops,pickup,y,kp1,q1,q2,tunetau){
  
  
 
    
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
    accepts[[1]]      = matrix(1,nstops,3)
    accepts[[2]]      = matrix(1,nstops,3)
    
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
    PT_chain[[1]]      = matrix(NA,niter,ncol=4)
    PT_chain[[2]]      = matrix(NA,niter,ncol=4)
    PT_chain[[1]][1,]  = c(.0006,.09,257,4)
    PT_chain[[2]][1,]  = c(.0006,.09,257,4)
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
  }else{
    
    #if algorithm fails, pick up from where it crushed   
    
    load(pickup)
    
    
    
    #theta             = out_ls$theta
    theta_1           = matrix(NA,niter,4)
    #some initial parameter values
      
    tau               = out_ls$tau
    tau_prop          = out_ls$tau_prop
    log_r_bot         = out_ls$log_r_bot
    log_r_bot1        = out_ls$log_r_bot1
    log_r_top1        = out_ls$log_r_top1
    accepts           = out_ls$accepts
    accepts1          = out_ls$accepts1
    tau_prop          = out_ls$tau_prop
    acc_rate          = out_ls$acc_rate
    acc_rate1         = out_ls$acc_rate1
    mllik             = out_ls$mllik
    failiter          = out_ls$failiter
    tune              = out_ls$tune
    tune_q1           = out_ls$tune_q1
    tune_q2           = out_ls$tune_q2
    kp                = out_ls$kp
   # mllik            = out_ls$mllik
    PT_chain         = out_ls$PT_chain
   # temp             = out_ls$temp
    temp             = matrix(NA,niter,2)
    swappers         = out_ls$swappers
    npar             = 2
    
  }
  return(list(  #theta             = theta,
                theta_1           = theta_1,
                tau               = tau,
                log_r_bot         = log_r_bot,
                log_r_bot1        = log_r_bot1,
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
                PT_chain          = PT_chain,
                temp              = temp,
                tune              = tune,
                tune_q1           = tune_q1,
                tune_q2           = tune_q2,
                kp                = kp,
                swappers          = swappers,
                npar              = npar))
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
runST=function(niter = 100,y,nstops=20,kp1=c(1,1,1),pickup=NULL,
               x=times,N=261,priorpars=c(1,1,N,5/N),q1=0.1,q2=0.1,kptau=1,tunetau){
  
    init=initialization(niter=niter,nstops=nstops,pickup=pickup,y,kp1,q1,q2,tunetau)
  
   
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
    
    R                 = Rhat(PT_chain[[1]][1,])
    log_r_bot[[1]]    = SIRlogpost(pars=PT_chain[[1]][1,],y=y,R=R[,2],I=R[,1],tau=tau[1],priorpars=priorpars,N=N)
    log_r_bot[[2]]    = log_r_bot[[1]]
      
    
    if (!(is.null(pickup))) {
      index            = which(is.na(init$tau))[1]
    }else{
      index            = 2
    }
    
    temp[index,]       = c(tau[index-1],1)
    

          
   out = list()
   out1 = list()
   library(parallel)
   cl = makeCluster(getOption("cl.cores", 4))
  
  for(iter in index:niter){
    
    print(iter)
    print(PT_chain[[1]][iter-1,])
    print(PT_chain[[2]][iter-1,])
    print(c(tau[iter-1],tau_prop[iter-1],tune[iter-1]))
    
    maybeswap = sample(0:(npar+1),2,replace=TRUE)
    
    if(all(length(unique(maybeswap))==2,min(maybeswap)>0,max(maybeswap)<=npar,iter>2)){
      #propose a parameter swap
      
   
      R1     = Rhat(PT_chain[[maybeswap[1]]][iter-1,])
      R2     = Rhat(PT_chain[[maybeswap[2]]][iter-1,])                                
          
      l11    = SIRloglik(y=y,R=R1[,2],I=R1[,1],tau=temp[iter-1,maybeswap[1]])[1]
      l22    = SIRloglik(y=y,R=R2[,2],I=R2[,1],tau=temp[iter-1,maybeswap[2]])[1]
      
      
      l21    = SIRloglik(y=y,R=R1[,2],I=R1[,1],tau=temp[iter-1,maybeswap[2]])[1]
      l12    = SIRloglik(y=y,R=R2[,2],I=R2[,1],tau=temp[iter-1,maybeswap[1]])[1]
      
     
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
 
    out                       =  SamplePars(y,theta=PT_chain[[chain]][iter-1,],tau=tau[iter-1],
                                            priorpars=priorpars,
                                            log_r_bot=log_r_bot[iter-1],N=N,
                                            accepts=c(accepts[[chain]][kp[[chain]],1],accepts[[chain]][kp[[chain]],2],accepts[[chain]][kp[[chain]],3]),
                                            q1=tune_q1[[chain]][iter-1],q2=tune_q2[[chain]][iter-1])
    
      
 
    PT_chain[[chain]][iter,]                 = out$theta
    accepts[[chain]][kp[[chain]][1],1]       = out$accepts[1]
    accepts[[chain]][kp[[chain]][2],2]       = out$accepts[2]
    accepts[[chain]][kp[[chain]][3],3]       = out$accepts[3]
    log_r_bot[[chain]][iter]                 = out$log_r_bot
    
    
 
    
    out1                = tryCatch({   
      
                                    STstep_tau(te=tau[iter-1],theta=PT_chain[[chain]][iter,],y, priorpars=priorpars,N=N,acc=accepts1[kp[[1]]],tune=tune[iter-1] ,iter=iter,cl)  
                                     
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
  
    out   = SamplePars(y,theta=PT_chain[[chain]][iter-1,],tau=1,priorpars=priorpars,
                       log_r_bot=log_r_bot[iter-1],N=N,
                       accepts=c(accepts[[chain]][kp[[chain]],1],accepts[[chain]][kp[[chain]],2],accepts[[chain]][kp[[chain]],3]),
                       q1=tune_q1[[chain]][iter-1],q2=tune_q2[[chain]][iter-1])
    
  
    #  print(out$theta)
    PT_chain[[chain]][iter,]                 = out$theta
    accepts[[chain]][kp[[chain]][1],1]       = out$accepts[1]
    accepts[[chain]][kp[[chain]][2],2]       = out$accepts[2]
    accepts[[chain]][kp[[chain]][3],3]       = out$accepts[3]
    log_r_bot[[chain]][iter]                 = out$log_r_bot
  
    }
  
  chain=1
  R                                    = Rhat(PT_chain[[chain]][iter,])
  mllik[[chain]][iter]                 = SIRloglik(y,R=R[,2],I=R[,1],tau[iter])[2]
  chain=2
  R                                    = Rhat(PT_chain[[chain]][iter,])
  mllik[[chain]][iter]                 = SIRloglik(y,R=R[,2],I=R[,1],tau[iter])[2]
  
  
  tune[iter]                           = tune[iter-1]  
  

for (chain in (1:npar)){
  
  tune_q1[[chain]][iter]=tune_q1[[chain]][iter-1]
  tune_q2[[chain]][iter]=tune_q2[[chain]][iter-1]
  
  
}

   
    if (((iter %% 1000) == 0) || (iter==10487) || (iter==10486)) {
      out_ls=list(PT_chain=PT_chain,tau=tau,accepts=accepts,accepts1=accepts1,log_r_bot=log_r_bot,tune_q1= tune_q1,tune_q2= tune_q2,
                  kp=kp,tau_prop=tau_prop,log_r_top1=log_r_top1,mllik=mllik,failiter=failiter,swappers=swappers,tune=tune)
      
        
      save(out_ls, file=paste("ST",iter,".RData",sep=""))
      save(".Random.seed",file=paste("random_state_",iter,".RData",sep=''))
    }
  }
  
stopCluster(cl)
  
  return(list(PT_chain=PT_chain,tau=tau,accepts=accepts,accepts1=accepts1,log_r_bot=log_r_bot,tune_q1= tune_q1,tune_q2= tune_q2,
              kp=kp,tau_prop=tau_prop,log_r_top1=log_r_top1,mllik=mllik,failiter=failiter,swappers=swappers,tune=tune))
}
######################################################################################
