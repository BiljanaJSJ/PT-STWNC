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


######################################################################################
#sample_tau
#input:         
#               y                 - data
#               means             - vector with last accepted value of means for each component
#               sigma2            - vector with last accepted value of sigma2 for each component
#               p                 - vector with last accepted value of mixture weights for each component 
#               Meanpriorpars     - prior for the means
#               Sigmapriorpars    - prior for sigmas
#               k                 - number of components
#               tau               - inverse temperature parameter
#               Ppriorpars        - prior for weights
#               acc               - counts accepted values of tau
#               Z                 - last sampled indicator (latent) variable 
#               tunetau           - tune the transitional step of tau
#               equalVar          - indicator that switches to the model with equal variances 
#                                   or to the model with unequal variance
# output:      
#               te                - last accepted value of tau
#               log_r_bot         - log of denominator of acceptance ratio when sampling tau
#               accepts1          - count of accepted tau
#               tau_prop          - proposed value of tau for the current iteration
#               log_r_top1        - log of numerator of acceptance ratio when sampling tau 
#               max_sigma2_prop   - optimized sigma2 value at tau_prop
#               max_sigma2_it     - optimized sigma2 value at current tau
######################################################################################
sample_tau=function(y,means,sigma2,p,Meanpriorpars,Sigmapriorpars,k,tau,Ppriorpars,acc,Z,tunetau,equalVar){
  
  #propose tau from the normal truncated on (0,1)
  tau_prop            = rtruncnorm(1,a=0,b=1,tau,tunetau)
  #obtain maximum values for mu and sigma2 for the proposed value of tau
  debugThis=list(tau=tau,tau_prop=tau_prop,means=means,sigma2=sigma2,p=p,Z=Z)
  save(debugThis,file='debug.RData')
  N=10000
  #initialize the objects
  max_means_prop=max_means_it=matrix(NA,N,length(means))
  max_means_prop[1,]=max_means_it[1,]=means
  max_sigma2_prop=max_sigma2_it=matrix(NA,N,length(sigma2))
  max_sigma2_prop[1,]=max_sigma2_it[1,]=sigma2
  p_max_prop=p_max_it=matrix(NA,N,length(p))
  p_max_prop[1,]=p_max_it[1,]=p
  Z_prop_max=Z_it_max=list()
  Z_prop_max[[1]]=Z_it_max[[1]]=Z
  Nk_prop=Nk_it=matrix(NA,N,length(means))
  Nk_prop[1,]=Nk_it[1,]=apply(Z,1,sum)
  post_prop=post_it=rep(NA,N)
  post_prop[1]=post_it[1]=0
  #mu is optimized using the posterior mean of the means formula 
  eps=10^-2
  for (i in (2:N)){
 
    
       #resample Z for the maximized means, sigma2 and p at tau_prop
           fmpij                  = mpij(y,tau=tau_prop,means=max_means_prop[i-1,],sigma2=max_sigma2_prop[i-1,],p=p_max_prop[i-1,])
           pij                    = exp(fmpij)/apply(exp(fmpij),1,sum)
           Z_prop_max[[i]]        = sampleZ(pij)
           Nk_prop[i,]            = apply(Z_prop_max[[i]],1,sum)

    
       #resample Z for the maximized means, sigma2 and p at current tau
           fmpij                  = mpij(y,tau=tau,means=max_means_it[i-1,],sigma2=max_sigma2_it[i-1,],p=p_max_it[i-1,])
           pij                    = exp(fmpij)/apply(exp(fmpij),1,sum)
           Z_it_max[[i]]          = sampleZ(pij)
           Nk_it[i,]              = apply(Z_it_max[[i]],1,sum)
       
       
       #maximize p for both proposed and current tau
           p_max_prop[i,]         = Nk_prop[i,]/sum(Nk_prop[i,])
           p_max_it[i,]           = Nk_it[i,]/sum(Nk_it[i,])
    
       #maximize means for tau_prop 
          sumY                    = Z_prop_max[[i]]%*%y
          denom_prop              = max_sigma2_prop[i-1,]+tau_prop*Nk_prop[i,]*Meanpriorpars[,2]
          max_means_prop[i,]      = (tau_prop*sumY*Meanpriorpars[,2]+max_sigma2_prop[i-1,]*Meanpriorpars[,1])/denom_prop
  
       #maximize sigma2 at proposed tau
       #sigma2 is optimized using the posterior mode
       SSE                        = SSE_f(y,Z_prop_max[[i]],means=max_means_prop[i,],k=k) 
  
       if (is.null(equalVar)){
           d0_prop                = (Nk_prop[i-1,]*tau_prop/2+Sigmapriorpars[1] )
           v0_prop                = (tau_prop*SSE/2+Sigmapriorpars[2])
           max_sigma2_prop[i,]    = v0_prop/(d0_prop+1)
       }else{
           SSE                    = sum(SSE)
           d0_prop                = (sum(Nk_prop[i-1,])*tau_prop/2+Sigmapriorpars[1] )
           v0_prop                = tau_prop*SSE/2+Sigmapriorpars[2]  
           max_sigma2_prop[i,]    = v0_prop/(d0_prop+1)
        }
    
  
        #maximize means and sigma2 for the current value of tau
          sumY                     = Z_it_max[[i]]%*%y
          denom_it                 = max_sigma2_it[i-1,]+tau*Nk_it[i,]*Meanpriorpars[,2]
          max_means_it[i,]         = (tau*sumY*Meanpriorpars[,2]+max_sigma2_it[i-1,]*Meanpriorpars[,1])/denom_it
      
        
        #optimize sigmas from the posterior mode at current tau
          SSE                      = SSE_f(y,Z_it_max[[i]],means=max_means_it[i,],k=k) 
      
        if (is.null(equalVar)){
           #if the model runs with unequal variances
           # the maximized sigma2 is a vector of length k  
            d0                     = (Nk_it[i-1,]*tau/2+Sigmapriorpars[1])
            v0                     = (tau*SSE/2+Sigmapriorpars[2])
            max_sigma2_it[i,]      = v0/(d0+1)
        }else{
           #if the model runs with equal variances 
           # the maximized sigma2 is a one value
            SSE                    = sum(SSE)
            d0                     = (sum(Nk_it[i-1,])*tau/2+Sigmapriorpars[1] )
            v0                     = (tau*SSE/2+Sigmapriorpars[2])
            max_sigma2_it[i,]      = v0/(d0+1)
        }

            #evaluate the posterior P(mu, sigma2,p,Z | y)  at the maximums for proposed and current tau
                      post_prop[i]  =  posterior_optim(y=y,Z=Z_prop_max[[i]],p=p_max_prop[i,],means=max_means_prop[i,],sigma2=max_sigma2_prop[i,],tau=tau_prop,Meanpriorpars,Sigmapriorpars,Ppriorpars,k=k)
                        post_it[i]  =  posterior_optim(y=y,Z=Z_it_max[[i]],p=p_max_it[i,],means=max_means_it[i,],sigma2=max_sigma2_it[i,],tau=tau,Meanpriorpars,Sigmapriorpars,Ppriorpars,k=k)
       #check if the relative error is smaller than the tolerance eps
        if (!(i==1)){
        if ((is.nan(post_prop[i])) || (is.nan(post_it[i])) || (is.nan(post_prop[i-1])) || (is.nan(post_it[i-1]))) break
        if (( (post_prop[i]-post_prop[i-1])/post_prop[i]<eps) & ( (post_it[i]-post_it[i-1])/post_it[i]<eps)) {
              break}
      
        }
  }

   #evaluate the joint posterior  P(means,sigma2,Z,p,tau | Y)  at proposed tau:
  log_r_top = posterior(y=y,Z=Z,p=p,
                        means=means,sigma2=sigma2,tau=tau_prop,
                        max_means=max_means_prop[i,],max_sigma2=max_sigma2_prop[i,],max_p=p_max_prop[i,],
                        max_Z=Z_prop_max[[i]],Meanpriorpars=Meanpriorpars,
                        Sigmapriorpars=Sigmapriorpars,Ppriorpars=Ppriorpars,k=k)
  
  #evaluate the joint posterior  P(means,sigma2,Z,p,tau | Y)  at current tau:
  log_r_bot = posterior(y=y,Z=Z,p=p,
                        means=means,sigma2=sigma2,tau=tau,
                        max_means=max_means_it[i,],max_sigma2=max_sigma2_it[i,],
                        max_p=p_max_it[i,],max_Z=Z_it_max[[i]],Meanpriorpars=Meanpriorpars,
                        Sigmapriorpars=Sigmapriorpars,Ppriorpars=Ppriorpars,k=k)
 
  #calculate the acceptance probability
  alpha=log_r_top$output-log_r_bot$output+log(dtruncnorm(tau,a=0,b=1,tau_prop,tunetau))-log(dtruncnorm(tau_prop,a=0,b=1,tau,tunetau))
  
  if (all(!is.na(alpha) , runif(1) < exp(alpha))){
    tau = tau_prop
    acc=acc+1
  }
 
  ptau=c(log_r_top$ptau,log_r_top$ptau)
  
  return(list(tau=tau,acc=acc,tau_prop=tau_prop,ptau=ptau))
}

######################################################################################

######################################################################################
#runST runs ST algorithm This is the main function
#
# input:       niter            - number of iterations
#              data             - Galaxy data from the MASS package
#              k                - number of components
#              Meanpriorpars    - prior for means
#              Sigmapriorpars   - prior for sigmas
#              Ppriorpars       - prior for p, the mixtures weigths
#              nstops           - number of stops in when calculating acceptance rate
#              kptau            - initialize the counter of the acceptance of tau
#              kp1              - initialize the counter of the acceptance of means and sigma2
#              adps_k           - initialize the tuning parameter of the transition variance of sigma2
#              adpm_k           - initialize the tuning parameter of the transition variance of means
#              equalVar         - indicator that switches to the model with equal variances 
#                                  or to the model with unequal variance 
# output:      means            - matrix of sampled means
#              sigma2           - matrix of sampled sigmas
#              p                - matrix of sampled weights for the mixture
#              Z                - list of sampled latent variable
#              accm             - counts of the accepted means
#              accs             - counts of the accepted sigma2
#              acc              - counts of the accepted tau 
#              tau              - vector of sampled inverse temperature
#              tau_prop         - for debugging proposed value of tau
#              mllik            - marginal likelihood
#              adps             - tuning parameter of the transitional variance of sigma2
#              tunetau          - tuning parameter for the transition step of tau
#              acc_rate         - acceptance rate of the means components
#              acc_rate1        - acceptance rate of the tau
#              acc_rate2        - acceptance rate of the sigma2 components
#              adpm             - tuning parameter of the transitional variance of means
#              swappers         - a matrix that keeps track of the swaps
######################################################################################

runST=function(niter,data,k,Meanpriorpars,Sigmapriorpars,Ppriorpars,nstops=20,kptau,kp1,adps_k,adpm_k,equalVar=NULL){
  
  #correct the typo in the Galaxy data set
  data[78]           = 26960 ## there's a typo in the dataset
  y                  = data/1000 ## normalize the data 
  n                  = length(y)
  #run PT-STWTDNC two chains
  #initialize the algorithm
 
  means                = list()
  means[[1]]           = matrix(NA,niter,k)
  means[[2]]           = matrix(NA,niter,k)
  means[[1]][1,]       = rep(20,k)
  means[[2]][1,]       = rep(20,k)
  sigma2               = list()  
  if (is.null(equalVar)){
  
  sigma2[[1]]             = matrix(NA,niter,k)
  sigma2[[2]]             = matrix(NA,niter,k)
  sigma2[[1]][1,]         = rep(0.5,k)
  sigma2[[2]][1,]         = rep(0.5,k)
  }else{
    
  sigma2[[1]]             = matrix(NA,niter,1)
  sigma2[[2]]             = matrix(NA,niter,1)
  sigma2[[1]][1]          = 0.5
  sigma2[[2]][1]          = 0.5
  }
  
  p                  = list()
  p[[1]]             = matrix(NA,niter,k)
  p[[2]]             = matrix(NA,niter,k)
  p[[1]][1,]         = rep(1/3,k)
  p[[2]][1,]         = rep(1/3,k)
  Z                  = list()
  
  tau                = rep(NA,niter)
  tau[1]             = 0.5
  mllik              = list()
  mllik[[1]]         = rep(NA,niter)
  mllik[[2]]         = rep(NA,niter)
  tau_prop           = rep(NA,niter)
  ptau               = matrix(NA,niter,2)
  
  kp                 = list()
  kp[[1]]            = kp1
  kp[[2]]            = kp1
  
  accm               = list()
  accm[[1]]          = matrix(1,nstops,k)
  accm[[2]]          = matrix(1,nstops,k)
  
  
  accs               = list()
  acc_rate           = list()
  adps               = list()
  
  if (is.null(equalVar)){
  accs[[1]]          = matrix(1,nstops,k)
  accs[[2]]          = matrix(1,nstops,k)
  acc_rate[[1]]      = matrix(1,nstops,k)
  acc_rate[[2]]      = matrix(1,nstops,k)
  adps[[1]]          = matrix(NA,niter,k)
  adps[[2]]          = matrix(NA,niter,k)
  adps[[1]][1,]      = adps_k
  adps[[2]][1,]      = adps_k
  }else{
    
  accs[[1]]          = rep(1,nstops)
  accs[[2]]          = rep(1,nstops)
  acc_rate[[1]]      = rep(1,nstops)
  acc_rate[[2]]      = rep(1,nstops)
  
  adps[[1]]          = rep(NA,niter)
  adps[[2]]          = rep(NA,niter)
  adps[[1]][1]       = adps_k
  adps[[2]][1]       = adps_k
  }
  
  acc                = rep(1,nstops)
  acc_rate1          = rep(1,nstops)
 
  acc_rate2          = list()
  acc_rate2[[1]]     = matrix(1,nstops,k)
  acc_rate2[[2]]     = matrix(1,nstops,k)
  
  
  adpm               = list()
  adpm[[1]]          = matrix(NA,niter,k)
  adpm[[2]]          = matrix(NA,niter,k)
  adpm[[1]][1,]      = adpm_k
  adpm[[2]][1,]      = adpm_k
  tunetau            = rep(NA,niter)
  tunetau[1]         = 2
    
 
  beta               = matrix(NA,niter,2)
  beta[1,]           = c(tau[1],1)
  npar               = 2
  swappers           = matrix(0,niter,3)
 
 
 
 
  var_means =rep(NA,k)
  
  for (iter in (2:niter)){
  
   # print(paste('iter=',iter,sep=''))
  
    #propose chains for a swap
    maybeswap = sample(0:(npar+1),2,replace=TRUE)
    
    if(all(length(unique(maybeswap))==2,min(maybeswap)>0,max(maybeswap)<=npar,iter>2)){
     
      
      #if swap starts calculate the exchange probability
      #and decide if the swap is accepted
      if (is.null(equalVar)){
      sigma2_1=sigma2[[maybeswap[1]]][iter-1,]
      sigma2_2=sigma2[[maybeswap[2]]][iter-1,]
      }else{
      sigma2_1=sigma2[[maybeswap[1]]][iter-1]
      sigma2_2=sigma2[[maybeswap[2]]][iter-1]
        
      }
      
      
      l11 = loglike(y,Z[[maybeswap[1]]],p[[maybeswap[1]]][iter-1,],means[[maybeswap[1]]][iter-1,],sigma2_1,beta[maybeswap[1]],k)$llik
      l22 = loglike(y,Z[[maybeswap[2]]],p[[maybeswap[2]]][iter-1,],means[[maybeswap[2]]][iter-1,],sigma2_2,beta[maybeswap[2]],k)$llik
     
      
      l21 = loglike(y,Z[[maybeswap[1]]],p[[maybeswap[1]]][iter-1,],means[[maybeswap[1]]][iter-1,],sigma2_1,beta[maybeswap[2]],k)$llik
      l12 = loglike(y,Z[[maybeswap[2]]],p[[maybeswap[2]]][iter-1,],means[[maybeswap[2]]][iter-1,],sigma2_2,beta[maybeswap[1]],k)$llik
      
      
      
      if(runif(1)< exp(l12 + l21   -l11 - l22)){
        #if swap is accepted
        #exchange the parameter of interests between two chains
        means[[maybeswap[1]]][iter,]    = means[[maybeswap[2]]][iter-1,]
            p[[maybeswap[1]]][iter,]    =     p[[maybeswap[2]]][iter-1,]
        
        means[[maybeswap[2]]][iter,]    = means[[maybeswap[1]]][iter-1,]
            p[[maybeswap[2]]][iter,]    =     p[[maybeswap[1]]][iter-1,]
        
        
        if (is.null(equalVar)){
         sigma2[[maybeswap[1]]][iter,]  = sigma2[[maybeswap[2]]][iter-1,]
         sigma2[[maybeswap[2]]][iter,]  = sigma2[[maybeswap[1]]][iter-1,]
        }else{
         sigma2[[maybeswap[1]]][iter]   = sigma2[[maybeswap[2]]][iter-1]
         sigma2[[maybeswap[2]]][iter]   = sigma2[[maybeswap[1]]][iter-1]
        }

        
        
       
        
        tmp    = Z[[1]]
    	  Z[[1]] = Z[[2]]
    	  Z[[2]] = tmp
        
        swappers[iter,]=c(sort(maybeswap),1)
        
      }else{
        # if no swap accepted
        # update the current parameters of interest by copying their 
        # last accepted values
        means[[maybeswap[1]]][iter,]   = means[[maybeswap[1]]][iter-1,]
            p[[maybeswap[1]]][iter,]   =     p[[maybeswap[1]]][iter-1,]
        
        means[[maybeswap[2]]][iter,]   = means[[maybeswap[2]]][iter-1,]
            p[[maybeswap[2]]][iter,]   =     p[[maybeswap[2]]][iter-1,]
        
        
        if (is.null(equalVar)){
          sigma2[[maybeswap[1]]][iter,]  = sigma2[[maybeswap[1]]][iter-1,]
          sigma2[[maybeswap[2]]][iter,]  = sigma2[[maybeswap[2]]][iter-1,]
        }else{
          sigma2[[maybeswap[1]]][iter]   = sigma2[[maybeswap[1]]][iter-1]
          sigma2[[maybeswap[2]]][iter]   = sigma2[[maybeswap[2]]][iter-1]
        }
   
        swappers[iter,]=c(sort(maybeswap),0)
      }
      
      tau[iter]                = tau[iter-1]
      beta[iter,]              = beta[iter-1,]
   
      
    }else{
      # Mutation step
      # Each of the two chains undergo the mutation step 
      
      #Mutation step for the 'tempered' chain via STWTDNC   
  chain=1  

  if (is.null(equalVar)){
      sigma2_it   =  sigma2[[chain]][iter-1,]
      accs_it     =  accs[[chain]][kp[[chain]][1],]
      adps_it     =  adps[[chain]][iter-1,]
  }else{
      sigma2_it   =  sigma2[[chain]][iter-1]
      accs_it     =  accs[[chain]][kp[[chain]][1]]
      adps_it     =  adps[[chain]][iter-1]
  }

  #Transition step 1: update the parameters while tau is held fixed
  #via Metropolis within Gibbs   
  out=Gibbs(y=y,means=means[[chain]][iter-1,],sigma2=sigma2_it,p=p[[chain]][iter-1,],Meanpriorpars=Meanpriorpars,
            tau=tau[iter-1],Sigmapriorpars=Sigmapriorpars,Ppriorpars=Ppriorpars,accm=accm[[chain]][kp[[chain]][2],],
            accs=accs_it,k=k,adps=adps_it,adpm=adpm[[chain]][iter-1,],equalVar=equalVar)
    
  means[[chain]][iter,]               = out$means
  
  if (is.null(equalVar)){
  sigma2[[chain]][iter,]              = out$sigma2
  accs[[chain]][kp[[chain]][1],]      = out$accs
  }else{
  sigma2[[chain]][iter]               = out$sigma2  
  accs[[chain]][kp[[chain]][1]]       = out$accs
  }
  
  p[[chain]][iter,]                   = out$p
  Z[[chain]]                          = out$Z
  accm[[chain]][kp[[chain]][2],]      = out$accm
  

  if (is.null(equalVar)){
    sigma2_it1                        = sigma2[[chain]][iter,]
  }else{
    sigma2_it1                        = sigma2[[chain]][iter]
  }

  #Transition step 2: update the tau while parameters of interest are being held fixed
  #via standard Metropolis-Hastings
  out1       = sample_tau(y=y,means=means[[chain]][iter,],sigma2=sigma2_it1,p=p[[chain]][iter,],Meanpriorpars=Meanpriorpars,
                          Sigmapriorpars=Sigmapriorpars,
                          k=k,tau=tau[iter-1],Ppriorpars=Ppriorpars,acc=acc[kptau],Z=Z[[chain]],
                          tunetau=tunetau[iter-1],equalVar=equalVar)

  tau[iter]                          = out1$tau
  beta[iter,]                        = c(tau[iter],1)
  acc[kptau]                         = out1$acc
  tau_prop                           = out1$tau_prop
  ptau[iter,]                        = out1$ptau
  

  
#mutation step for the 'target' chain:
#update the parameters of interest at tau=1 via standard Metrpolis-Hastings
chain=2
 
  if (is.null(equalVar)){
    sigma2_it                        = sigma2[[chain]][iter-1,]
    accs_it                          = accs[[chain]][kp[[chain]][1],]
    adps_it                          = adps[[chain]][iter-1,]

  }else{
    sigma2_it                        = sigma2[[chain]][iter-1]
    accs_it                          = accs[[chain]][kp[[chain]][1]]
    adps_it                          = adps[[chain]][iter-1]

  }

  out = Gibbs(y=y,means=means[[chain]][iter-1,],sigma2=sigma2_it,p=p[[chain]][iter-1,],Meanpriorpars=Meanpriorpars,
              tau=1,Sigmapriorpars=Sigmapriorpars,Ppriorpars=Ppriorpars,accm=accm[[chain]][kp[[chain]][2],],
              accs=accs_it,k=k,adps=adps_it,adpm=adpm[[chain]][iter-1,],equalVar=equalVar)

 
 means[[chain]][iter,]                 = out$means
 
 if (is.null(equalVar)){
   sigma2[[chain]][iter,]              = out$sigma2
   accs[[chain]][kp[[chain]][1],]      = out$accs
 }else{
   sigma2[[chain]][iter]               = out$sigma2  
   accs[[chain]][kp[[chain]][1]]       = out$accs
 }
 p[[chain]][iter,]                     = out$p
 Z[[chain]]                            = out$Z
 accm[[chain]][kp[[chain]][2],]        = out$accm
 }


#tune transitional variance of means and sigma2 for each of the two chains
for (chain in (1:npar)){
  
  
  if (is.null(equalVar)){
    adps[[chain]][iter,]               = adps[[chain]][iter-1,]
  }else{
    adps[[chain]][iter]                = adps[[chain]][iter-1]
  }
  
  adpm[[chain]][iter,]                 = adpm[[chain]][iter-1,]
  
  if (iter== kp[[chain]][1]*niter/nstops){  
    
  
    
    if (is.null(equalVar)){
      acc_rate[[chain]][kp[[chain]][1],] = (accs[[chain]][kp[[chain]][1],]+1)/(2+(niter/nstops));
    }else{
      acc_rate[[chain]][kp[[chain]][1]]  = (accs[[chain]][kp[[chain]][1]]+1)/(2+(niter/nstops));  
    }
    kp[[chain]][1]                       =  kp[[chain]][1]+1;
    if(iter< niter/2){
      if (is.null(equalVar)){
        if(acc_rate[[chain]][kp[[chain]][1]-1,]>.29 || acc_rate[[chain]][kp[[chain]][1]-1,]<.19 ){
          # adjust the transition density if we are in the first half of the iterations. 
          adps[[chain]][iter,]           = adps[[chain]][iter,]*acc_rate[[chain]][kp[[chain]][1]-1,]/.24
        }}else{
          
          if(acc_rate[[chain]][kp[[chain]][1]-1]>.29 || acc_rate[[chain]][kp[[chain]][1]-1]<.19 ){
            # adjust the transition density if we are in the first half of the iterations. 
            adps[[chain]][iter]          = adps[[chain]][iter]*acc_rate[[chain]][kp[[chain]][1]-1]/.24
          }
        }
        
      }
    }
  
  
  if (iter == kp[[chain]][2]*niter/nstops){  
    acc_rate2[[chain]][kp[[chain]][2],]  = (accm[[chain]][kp[[chain]][2],]+1)/(2+(niter/nstops));
    kp[[chain]][2]                       =  kp[[chain]][2]+1;
    # adjust the transition density if we are in the first half of the iterations. 
    if(iter< niter/2){
      if(acc_rate2[[chain]][kp[[chain]][2]-1,]>.29 || acc_rate2[[chain]][kp[[chain]][2]-1,]<.19 ){
         adpm[[chain]][iter,]            = adpm[[chain]][iter,]*acc_rate2[[chain]][kp[[chain]][2]-1,]/.24
      }
    }
    
  }
 }
#tune transitional step for tau
tunetau[iter]                            = tunetau[iter-1]

if (iter == kptau*niter/nstops){  
        acc_rate1[kptau]                 = (acc[kptau]+1)/(2+(niter/nstops));
  kptau                                  =  kptau+1;
  # adjust the transition density if we are in the first half of the iterations. 
  if(iter< niter/2){
    if(acc_rate1[kptau-1]>.19 || acc_rate1[kptau-1]<.9 ){
          tunetau[iter]                  = tunetau[iter]*acc_rate1[kptau-1]/.14
    }
  }
}

#calculate marginal likelihood
for (chain in (1:npar)){
   if (is.null(equalVar)){
       sigma2_mllik                     = sigma2[[chain]][iter,]
   }else{
       sigma2_mllik                     = sigma2[[chain]][iter]
   }
  mllik[[chain]][iter]                  = loglike(y=y,Z=Z[[chain]],p=p[[chain]][iter,],means=means[[chain]][iter,],sigma2=sigma2_mllik,tau=tau[iter],k=k)$mllik
 
 }
 
 

 
  if ((iter %% 1000) == 0) {
    out_ls  = list(means=means,sigma2=sigma2,p=p,Z=Z,accm=accm,accs=accs,tau=tau,acc=acc,tau_prop=tau_prop,mllik=mllik,adps=adps,tunetau=tunetau,acc_rate=acc_rate,acc_rate1=acc_rate1,acc_rate2=acc_rate2,adpm=adpm,swappers=swappers)

    save(out_ls, file=paste("ST",iter,".RData",sep=""))
  }
  }


  return(list(means=means,sigma2=sigma2,p=p,Z=Z,accm=accm,accs=accs,tau=tau,acc=acc,tau_prop=tau_prop,mllik=mllik,adps=adps,tunetau=tunetau,acc_rate=acc_rate,acc_rate1=acc_rate1,acc_rate2=acc_rate2,adpm=adpm,swappers=swappers))
  
  }
  
  
  
  
  
  
  
  
