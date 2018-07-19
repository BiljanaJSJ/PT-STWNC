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


library("MASS")
library(coda)

surf.colors <- function(x, col = terrain.colors(20)) {
  # First we drop the 'borders' and average the facet corners
  # we need (nx - 1)(ny - 1) facet colours!
  x.avg <- (x[-1, -1] + x[-1, -(ncol(x) - 1)] +
              x[-(nrow(x) -1), -1] + x[-(nrow(x) -1), -(ncol(x) - 1)]) / 4
  # Now we construct the actual colours matrix
  colors = col[cut(x.avg, breaks = length(col), include.lowest = T)]
  return(colors)
}
#######################################################
#prior_mu
# mu is N(0,k^2) distributed
#input:         k- standard deviation of mu
#output:        data point for mu from its prior
#######################################################

dprior_mu=function(x,mean_mu=0,k=1,log=T){
  return(sum(dnorm(x,mean=mean_mu,sd=k,log=log)))
}

############################################################
# logpriorSigma:Prior of Sigmas 
# sigmas are distributed IG(priorpars[1],priorpars[2])
# input:     x              - values at which prior of sigma2 is evaluated
#            SigmaPriorPars - vector of shape and scale of the prior sigma2
#            log            - whether log of the prior is evaluated, 
#                             default is true
#output:     evaluated prior of sigma2 for each of the components
############################################################
dprior_sig=function(x,SigmaPriorPars,log=T){
   return(sum(dgamma(1/x,shape=SigmaPriorPars[1],scale=1/SigmaPriorPars[2],log=log)))
}

#############################################################
# loglik:  Tempered Likelihood is Gaussian
# here the tempered likelihood is evaluated, as well as
# the untempered log likelihood which later
# is used to calculate the log marginal likelihood
# input:  x           - data
#         mean_mu     - mean
#         sigma       - sigma2
#         tau         - value of tau
#         log         - whether log of the full likelihood is evaluated, 
#                       default is true
#         switchllik  - switch the tempered likelihood from 
#                       power tempered into tau in sigma
# output:   a list with two elements: llik - tempered likelihood
#                                     mllik - untempered likelihood
#############################################################
loglik=function(x,mean_mu,sigma=1,tau,log=T, switchllik){
  if (log==T)
  {
    if (switchllik=='power'){
    out=tau*sum(dnorm(x,mean=abs(mean_mu),sd=sqrt(sigma),log=log))
	  }else{
	  out=sum(dnorm(x,mean=abs(mean_mu),sd=sqrt(sigma/tau),log=log))
	  }
	 
    }else{
    if (switchllik=='power'){
    out=(prod(dnorm(x,mean=abs(mean_mu),sd=(sigma),log=log))^tau)
	  }else{
	  out=(prod(dnorm(x,mean=abs(mean_mu),sd=(sigma/tau),log=log)))
	  }
	  }
    mllik=sum(dnorm(x,mean=abs(mean_mu),sd=sqrt(sigma),log=log)) 
    return(list(out=out,mllik=mllik))

}


#######################################################
#data - simulate data y from N(abs(mu),sigma^2)
#input:         n - sample size
#            mean - mean of the data
#           sigma - standard deviation of the data
#output:       data from N(abs(mu),sigma^2) where mu
# is simulated with prior_mu function
#######################################################
data=function(n=5,mu,sigma=1){
  
  return(rnorm(n=n,mean=abs(mu),sd=sigma))
}
#######################################################


######################################################################################
#STStep_mu
#input:       th        - bivariate data point (mu,tau)
#             y         - data vector
#             n         - sample size
#             p_mu      - true mu value
#             k         - stdev of prior on mu (=1)
#             sigma     - stdev of the data vector y (=1)
#             mean_mu   - prior mean of mu (=0)
#             rate_exp  - rate of prior of tau (=0.1)
#             log_r_bot - this is passed parameter so that we calculate it only if 
#                         there are changes
#             accepts   - keeps track of the number of accepted values
#             kp        - counter for accepted values
#             rexp      - prior(tau)=exp(rexp)
#output:      theta     - data point estimated from MH
#             accepts   - number of accepted proposals
#             log_r_bot - un-normalized log posterior at the last accepted 
#                         parameter value
#             kp=kp     - a counter for acceptance rate

######################################################################################
STstep_pars = function(th,y=y,n=n,p_mu=p_mu,k=k,mean_mu=mean_mu,
                       log_r_bot=NA,acc=accepts,niter=niter,SigmaPriorPars,tune=tune,switchllik,tune_q1,post_sd){
  

  # th is the last accepted parameter value
   #tune=2.38
   #calculate the optimal posterior variance of mu
   #post_sd    = (((th[3])*k^2)/(n*th[2]*k^2+th[3]))
   #post_sd    = sqrt(post_sd*tune/2)
   var_target   = th[3]/(n*th[2]+th[3]) 
  
   #th is the last accepted parameter value
   post_sd    = 2.4*sqrt(var_target)
  
   #propose mu
   mu_prop=rnorm(n=1,mean=mean_mu,sd=post_sd)

   #the top part of alpha
   log_r_top_mu=posterior_mu(n=n,y=y,tau=th[2],mu=mu_prop,k=k,sigma=th[3],log=T,switchllik=switchllik)
  
   #the last accepted parameter value
   log_r_bot_mu=posterior_mu(n=n,y=y,tau=th[2],mu=th[1],k=k,sigma=th[3],log=T,switchllik=switchllik)
  
  
   #calculate alpha
   alpha_mu = log_r_top_mu - log_r_bot_mu
   # make a decision
  
   if (all(!is.na(alpha_mu) , runif(1) < exp(alpha_mu))){
    # accept the move
    acc[1]=acc[1]+1;
    th[1] = mu_prop;
    log_r_bot_mu = log_r_top_mu;
    
   }
    
#   #propose sigma2
#   #sig_prop=abs(th[3]+runif(1,-0.05,0.05))
#   # q1=0.05
#    delta           = rnorm(1,0,tune_q1)
#    logsig_prop     = log(abs(th[3]))+delta
##   logsig_prop     = rlnorm(1,th[3],q1)
#    sig_prop        = exp(logsig_prop)

## print('th[3],logsig_prop,sig_prop')
## print(c(th[3],logsig_prop,sig_prop))
#   #sig_prop=abs(rnorm(n=1,mean=th[3],sd=post_sig_sd))
#   #the top part of alpha
#   log_r_top_sig=posterior_sig(sigma=sig_prop,n=n,y=y,tau=th[2],mu=th[1],SigmaPriorPars=SigmaPriorPars,log=T,switchllik=switchllik)
   
#   #the last accepted parameter value
#   log_r_bot_sig=posterior_sig(sigma=th[3],n=n,y=y,tau=th[2],mu=th[1],SigmaPriorPars=SigmaPriorPars,log=T,switchllik=switchllik)
   
   log_r_bot_sig=0
   
#   #calculate alpha
#   alpha_sig = log_r_top_sig - log_r_bot_sig +dlnorm(sig_prop,log(th[3]),tune_q1,log=T)-dlnorm(th[3],logsig_prop,tune_q1,log=T)
#   # make a decision
   
#   if (all(!is.na(alpha_sig) , runif(1) < exp(alpha_sig))){
#     # accept the move
#     acc[2]=acc[2]+1;
#     th[3] = sig_prop;
#     #log_r_bot_sig = log_r_top_sig;
#     
#   }
#  # do not do anything since we maintain the last value
  return(list(theta=th,accepts=acc,log_r_bot=c(log_r_bot_mu,log_r_bot_sig),mu1prop=mu_prop))
}
#####################################################################################
#posterior_mu function 
#draws samples from the posterior 
#        distribution of (mu/Y,tau)
#input:  n              - sample size
#        y              - univariate data
#        tau            - temperature parameter
#        mu             - abs(mu) is mean of the data
#        k              - standard deviation of the mean
#        sigma          - variance of the data
#        log            - whether log of the full likelihood is evaluated, 
#                         default is true
#        switchllik     - switch the tempered likelihood from 
#                         power tempered into tau in sigma
#output: draws samples from the conditional posterior of mu 
#        distribution (mu/Y,tau), given the input parameters
#####################################################################################

posterior_mu=function(n=5,y,tau,mu,k=1,sigma=1,log=T,switchllik){

    llik  = loglik(y,mean_mu=mu,sigma=sigma,tau=tau,log=log,switchllik=switchllik)$out
    pr_mu = dprior_mu(x=mu,k=k,log=log)
    
    
    if (log==T){ 
      out_mu=llik+pr_mu
    }else{
      out_mu=llik*pr_mu
    }
  
  return(out_mu)
}

#####################################################################################
#posterior_sig function 
#draws samples from the posterior 
#        distribution of (mu/Y,tau)
#input:  
#        sigma          - variance of the data
#        n              - sample size
#        y              - univariate data
#        tau            - temperature parameter
#        mu             - abs(mu) is mean of the data
#        log            - whether log of the full likelihood is evaluated, 
#                         default is true
#        SigmaPriorPars - shape and scale of the prior of sigma2
#        switchllik     - switch the tempered likelihood from 
#                         power tempered into tau in sigma
#output: draws samples from the conditional posterior of mu 
#        distribution (mu/Y,tau), given the input parameters
#####################################################################################
posterior_sig=function(sigma,n=5,y,tau,mu,log=T,SigmaPriorPars,switchllik){
  
  llik  = loglik(y,mean_mu=mu,sigma=sigma,tau=tau,log=log,switchllik=switchllik)$out
  pr_sig = dprior_sig(x=sigma,SigmaPriorPars=SigmaPriorPars,log=log)
  
  if (log==T){ 
    out_sig=llik+pr_sig
  }else{
    out_sig=llik*pr_sig
  }
  
  return(out_sig)
}
#######################################################################
# prior_tau: calculate the prior of tau
# input:  x              - data
#         mean_mu        - mean
#         sigma          - sigma2
#         Sigmapriorpars - prior of sigma2
#         tau            - value of tau
#         k              - prior variance for the mean
#         log            - whether log of the full likelihood is evaluated, 
#                          default is true
#         switchllik     - switch the tempered likelihood from 
#                          power tempered into tau in sigma
# output: a value of calculated prior of tau
########################################################################
prior_tau=function(x,mean_mu,sigma=1,SigmaPriorPars,tau,k=1,log=T,switchllik){
  
  llik=loglik(x,mean_mu=mean_mu,sigma=sigma,tau=tau,log=log,switchllik=switchllik)$out
  pr_mu=dprior_mu(x=mean_mu,k=k,log=log)
  #pr_sig = dprior_sig(x=sigma,SigmaPriorPars=SigmaPriorPars,log=log)
  
  if (log==T){
    prior_t=-llik-pr_mu#-pr_sig
  }else{
    prior_t=1/(llik*pr_mu*pr_sig)
  }
  
  return(sum(prior_t ))
}


#######################################################################
# posterior: calculate the prior of tau
# input:  n              - sample size
#         y              - data
#         tau            - value of tau
#         mu             - mean
#         k              - prior variance for the mean
#         sigma          - sigma2
#         log            - whether log of the full likelihood is evaluated, 
#                          default is true
#         Sigmapriorpars - prior of sigma2
#         max_mu         - maximized mu
#         max_sig        - maximized sigma
#         switchllik     - switch the tempered likelihood from 
#                          power tempered into tau in sigma
# output: a value of evaluated posterior
########################################################################
posterior=function(n=5,y,tau,mu,k=1,sigma=1,log=T,SigmaPriorPars,max_mu,max_sig,switchllik){
  
  ptau  = prior_tau(y,mean_mu=max_mu,sigma=sigma,SigmaPriorPars=SigmaPriorPars,
                    tau=tau,k=k,log=log,switchllik=switchllik)
  llik  = loglik(y,mean_mu=mu,sigma=sigma,tau=tau,log=log,switchllik=switchllik)$out
  pr_mu = dprior_mu(x=mu,k=k,log=log)
#  pr_sig= dprior_sig(x=sigma,SigmaPriorPars=SigmaPriorPars,log=log)
  if (log==T){
    output=llik+pr_mu+ptau#+pr_sig
  }else{
    output=llik*pr_mu*ptau
  }
  return(output)
}

#############################################################
# Compute the SSE between the data and the mean
#############################################################
SSEfun = function(y,mu){  
  sum((y-abs(mu))^2)
}

######################################################################################
#STStep_tau: Transition two, update tau with fixed parameters
#input:       th              - bivariate data point (mu,tau)
#             y               - data vector
#             n               - sample size
#             k               - stdev of prior on mu (=1)
#             acc             - count the accepted tau
#             SigmaPriorPars  - shape and scale for the prior of sigma2
#             switchllik      - switch the tempered likelihood from 
#                               power tempered into tau in sigma
#          
#output:      theta           - sampled tau
#             log_r_bot       - denominator of the acceptance ratio
#             accepts1        - count of accepted tau
#             tau_prop        - proposed tau
#             log_r_top1      - numerator of the acceptance ratio
######################################################################################
STstep_tau = function(th,y=y,n=n,k=k,
                      acc=accepts1,SigmaPriorPars,switchllik){
  
  tau_prop        = runif(n=1,0,1)

  y_bar           = mean(y)
  
  
  N=10000
  max_sigma2_prop=max_sigma2_it= rep(NA,N)
  max_sigma2_prop[1]=max_sigma2_it[1]=th[3]
  Optimx_prop=Optimx_it= rep(NA,N)
  Optimx_prop=Optimx_it=th[1]
  #max_sigma2_prop_old = max_sigma2_it_old=max_sigma2_prop
  
  eps=10^-3
  for (i in (2:N)){
    
    #maximize mean from the posterior mean at tau_prop
    Optimx_prop[i]     = n*y_bar*((tau_prop*k^2)/(n*tau_prop*k^2+max_sigma2_prop[i-1]))
    #maximize mean from the posterior mean at current tau
    Optimx_it[i]       = n*y_bar*((th[2]*k^2)/(n*th[2]*k^2+max_sigma2_it[i-1]))
    
    
    #maximize the sigma2 from the posterior mode at tau_prop
    SSE             = SSEfun(y,Optimx_prop[i] ) 
    d0_prop         = ifelse(switchllik=='power',n*tau_prop/2+SigmaPriorPars[1],n/2+SigmaPriorPars[1])
    v0_prop         = 1/(tau_prop*SSE/2+1/SigmaPriorPars[2])
    max_sigma2_prop[i] = 1/((d0_prop+1)*v0_prop)
    
    SSE             = SSEfun(y,Optimx_it[i] ) 
    #maximize the sigma2 from the posterior mode at current tau
    d0              = ifelse(switchllik=='power',n*th[2]/2+SigmaPriorPars[1],n/2+SigmaPriorPars[1])
    v0              = 1/(th[2]*SSE/2+1/SigmaPriorPars[2])
    max_sigma2_it[i]   = 1/((d0+1)*v0)
    
    if (!(i==1)){
      if (( (abs(Optimx_prop[i]-Optimx_prop[i-1])/Optimx_prop[i])<eps) & ( abs( (Optimx_it[i]-Optimx_it[i-1])/Optimx_it[i])< eps) & (abs( (max_sigma2_prop[i]-max_sigma2_prop[i-1])/max_sigma2_prop[i])< eps) & (abs( (max_sigma2_it[i]-max_sigma2_it[i-1])/max_sigma2_it[i])< eps)){
        break}
    }
    #     Optimx_prop_old       = Optimx_prop
    #     Optimx_it_old         = Optimx_it
    #     max_sigma2_prop_old   = max_sigma2_prop
    #     max_sigma2_it_old     = max_sigma2_it
    
    
  }
  print(i)  
  
  #calculate the MH acceptance prob
  log_r_top       = posterior(n=n,y=y,tau=tau_prop,mu=th[1],k=k,sigma=th[3],SigmaPriorPars=SigmaPriorPars,log=T,max_mu=Optimx_prop[i],max_sig=max_sigma2_prop[i],switchllik=switchllik)
  log_r_bot       = posterior(n=n,y=y,tau=th[2],mu=th[1],k=k,sigma=th[3],SigmaPriorPars=SigmaPriorPars,log=T,max_mu=Optimx_it[i],max_sig=max_sigma2_it[i],switchllik=switchllik)
  
  #calculate alpha
  alpha           = log_r_top - log_r_bot
  # make a decision
  
  if (all(!is.na(alpha) , runif(1) < exp(alpha))){
    # accept the move
    acc           = acc+1;
    th[2]         = tau_prop
  
  }
  
  
  return(list(theta=th,
              log_r_bot = log_r_bot,accepts1=acc,
              log_r_top1=log_r_top,tau_prop=tau_prop)) #RHS=RHS
}
######################################################################################
initialization=function(niter=niter,nstops=nstops,pickup=pickup,kp1,q1,q2,k){
  
  if (is.null(pickup)) {
    
    #intialize the algorithm for the first run
    log_r_bot        = matrix(0,niter,2)
    log_r_bot1       = rep(0,niter)
    log_r_top1       = rep(0,niter)
    accepts          = list()
    accepts[[1]]     = matrix(1,nstops,2)
    accepts[[2]]     = matrix(1,nstops,2)
    
        
    acc_rate          = list()
    acc_rate[[1]]     = matrix(1,nstops,2)
    acc_rate[[2]]     = matrix(1,nstops,2)
    
    kp    = list()
    kp[[1]] = kp1
    kp[[2]] = kp1
    
    accepts1         = rep(1,nstops)
   # acc_rate         = rep(1,nstops)
    acc_rate1        = rep(1,nstops)
    #theta            = matrix(NA,niter,ncol=3)
 
#    theta[1,]        = c(-0.8,0.5,0.2)
    tau_prop         = rep(0,niter) 
    
    #mllik            = rep(NA,niter)
tune_q1           = list()
tune_q1[[1]]      = rep(NA,niter)
tune_q1[[2]]      = rep(NA,niter)
tune_q1[[1]][1]   = q2
tune_q1[[2]][1]   = q2

tune           = list()
tune[[1]]      = rep(NA,niter)
tune[[2]]      = rep(NA,niter)
tune[[1]][1]   = q1
tune[[2]][1]   = q1


    
    
    PT_chain           = list()
    PT_chain[[1]]      = matrix(NA,niter,ncol=3)
    PT_chain[[2]]      = matrix(NA,niter,ncol=3)
    PT_chain[[1]][1,]  = c(-0.8,0.5,1)
    PT_chain[[2]][1,]  = c(-0.8,1,1)
    PT_chain[[2]][,2]  = 1
 #   beta               = matrix(NA,niter,2)
  #  beta[1,]           = c(PT_chain[[1]][1,2],1)
   
    post_sd             = list()
    post_sd[[1]]        = rep(NA,niter) 
    post_sd[[2]]        = rep(NA,niter) 
    post_sd[[1]][1]     = (((PT_chain[[1]][1,3])*k^2)/(n*PT_chain[[1]][1,2]*k^2+PT_chain[[1]][1,3]))
    post_sd[[2]][1]     = (((PT_chain[[2]][1,3])*k^2)/(n*PT_chain[[2]][1,2]*k^2+PT_chain[[2]][1,3]))
    
    npar               = 2
    swappers           = matrix(0,niter,3)
    mllik              = list()
    mllik[[1]]         = rep(NA,niter)
    mllik[[2]]         = rep(NA,niter)
    
  }else{
    
    #if algorithm fails, pick up from where it crushed   
    
    load(pickup)
    
    
    log_r_bot        = out_ls$log_r_bot
    log_r_bot1       = out_ls$log_r_bot1
    log_r_top1       = out_ls$log_r_top1
    accepts          = out_ls$accepts
    accepts1         = out_ls$accepts1
    acc_rate         = out_ls$acc_rate
    acc_rate1        = out_ls$acc_rate1
  #  theta            = out_ls$theta
    tau_prop         = out_ls$tau_prop
    mllik            = out_ls$mllik
    PT_chain         = out_ls$PT_chain
  #  beta             = out_ls$beta
    swappers         = out_ls$swappers
    npar             = out_ls$npar
    kp               = out_ls$kp
    tune             = out_ls$tune
    tune_q1          = out_ls$tune_q1
  post_sd             = out_ls$post_sd
  
  }
  return(list(log_r_bot        = log_r_bot,
              log_r_bot1       = log_r_bot1,
              log_r_top1       = log_r_top1,
              accepts          = accepts,
              accepts1         = accepts1,
              acc_rate         = acc_rate,
              acc_rate1        = acc_rate1,
           #   theta            = theta,
              tau_prop         = tau_prop,
              mllik            = mllik,
              PT_chain         = PT_chain,
           #   beta             = beta,
              swappers         = swappers,
              npar             = npar,
              kp               = kp,
              tune             = tune,
              tune_q1          = tune_q1,
              post_sd=post_sd))
}
######################################################################################
#runST runs ST algorithm.
#
# input:       niter            - number of iterations
#              n                - sample size
#              p_mu             - true mean value for the data
#              k                - stdev of the prior of mu (=1), this is adapted
#              sigma            - stdev of the data (=1)
#              nstops           - number of stops in case we are tuning
#              kp               - count for the number of stops 
#              pickup           - default NULL, refering to that we are starting the algorithm from the beggining
#                                 if we pickup from where it failed we just need to supply the last successfully written file
#              SigmaPriorPars   - prior for sigmas
#              switchllik       - switch the tempered likelihood from 
#                                 power tempered into tau in sigma



# output:      theta            - bivariate chain with length niter
#              acc_rate         - aceptance rate
#              accepts          - number of accepted proposals
#              log_r_bot        - vector of length niter -(un-normalized log posterior 
#                                 at the last accepted parameter value)
#              kp               - counter for accepted values
#              y                - data vector
#              tau_profile      - two chains for the profile likelihood of tau
#                                 one obtained from the proposed tau and the other from
#                                 the current iteration (they should be the same)
######################################################################################
runST=function(niter = 500,n=25,p_mu=1.5,k=1,sigma=1,nstops=20,kp1=c(1,1),pickup=NULL,kptau=1,SigmaPriorPars,switchllik,q1=2.4,q2=0.05){
  
  #y=data(n=n,mu=abs(p_mu),sigma=sigma)

  #y=get(load('data.RData'))
  #y=get(load('/Users/bjonoska/Documents/runPants/ST50000_power.RData'))$y
  y=get(load('ST35000_power.RData'))$y


  init=initialization(niter=niter,nstops=nstops,pickup=pickup,kp1=kp1,q1=q1,q2=q2,k=k)
  
  log_r_bot          = init$log_r_bot
  log_r_bot1         = init$log_r_bot1
  log_r_top1         = init$log_r_top1
  accepts            = init$accepts
  accepts1           = init$accepts1
  acc_rate           = init$acc_rate
  acc_rate1          = init$acc_rate1
 # theta            = init$theta
  tau_prop           = init$tau_prop  
  mllik              = init$mllik

  
  PT_chain           = init$PT_chain
 # beta               = init$beta
  swappers           = init$swappers
  npar               = init$npar
  kp                 = init$kp
  tune               = init$tune
  tune_q1            = init$tune_q1
  post_sd            = init$post_sd
  
  if (!(is.null(pickup))) {
    index=which(is.na(init$PT_chain[[1]][,2]))[1]
  }else{
    index            = 2
  }
  
   
  #tune      = 2.38
  
  #initialize log_r_bot
  log_r_bot[1,1]  =  posterior_mu(n=n,y=y,tau=PT_chain[[1]][1,2],mu=PT_chain[[1]][1,1],
                                  k=PT_chain[[1]][1,3],log=T,switchllik=switchllik)
  
  log_r_bot[1,2]  =  posterior_sig(n=n,y=y,tau=PT_chain[[1]][1,2],mu=PT_chain[[1]][1,1],
                                   sigma=PT_chain[[1]][1,3],log=T,SigmaPriorPars=SigmaPriorPars,
                                   switchllik=switchllik)
  
  out = list()
 
  

  for(iter in index:niter){
    
    
    
    maybeswap = sample(0:(npar+1),2,replace=TRUE)
    
    if(all(length(unique(maybeswap))==2,min(maybeswap)>0,max(maybeswap)<=npar,iter>2)){
      #propose a parameter swap
      l11=loglik(y,mean_mu=PT_chain[[maybeswap[1]]][iter-1,1],sigma=PT_chain[[maybeswap[1]]][iter-1,3],tau=PT_chain[[maybeswap[1]]][iter-1,2],log=T, switchllik=switchllik)$out
      l22=loglik(y,mean_mu=PT_chain[[maybeswap[2]]][iter-1,1],sigma=PT_chain[[maybeswap[2]]][iter-1,3],tau=PT_chain[[maybeswap[2]]][iter-1,2],log=T, switchllik=switchllik)$out
        
      l21=loglik(y,mean_mu=PT_chain[[maybeswap[1]]][iter-1,1],sigma=PT_chain[[maybeswap[1]]][iter-1,3],tau=PT_chain[[maybeswap[2]]][iter-1,2],log=T, switchllik=switchllik)$out
      l12=loglik(y,mean_mu=PT_chain[[maybeswap[2]]][iter-1,1],sigma=PT_chain[[maybeswap[2]]][iter-1,3],tau=PT_chain[[maybeswap[1]]][iter-1,2],log=T, switchllik=switchllik)$out

      #       print(c(maybeswap[1],PT_chain[[maybeswap[1]]][iter-1,],beta[iter-1,maybeswap[1]]))
      #       print(c(l11,l22,l21,l12))
      #       print(iter)
      if(runif(1)< exp(l12 + l21 -l11 - l22)){
        #accept the swap
        
        PT_chain[[maybeswap[1]]][iter,c(1,3)] =PT_chain[[maybeswap[2]]][iter-1,c(1,3)]
        PT_chain[[maybeswap[2]]][iter,c(1,3)] =PT_chain[[maybeswap[1]]][iter-1,c(1,3)]
      
        swappers[iter,]=c(maybeswap,1)
        
      }else{
        # no swap accepted
        PT_chain[[maybeswap[1]]][iter,c(1,3)] =PT_chain[[maybeswap[1]]][iter-1,c(1,3)]
        PT_chain[[maybeswap[2]]][iter,c(1,3)] =PT_chain[[maybeswap[2]]][iter-1,c(1,3)]
        
        
        swappers[iter,]=c(maybeswap,0)
      }
   
      PT_chain[[1]][iter,2] = PT_chain[[1]][iter-1,2]
      PT_chain[[2]][iter,2] = 1
      
    }else{
   
   
   chain=1
   out          =      STstep_pars(th=PT_chain[[chain]][iter-1,],n=n,p_mu=p_mu,k=k,
                                   mean_mu=PT_chain[[chain]][iter-1,1],
                                   y=y,log_r_bot=log_r_bot[iter-1,],
                                   acc=c(accepts[[chain]][kp[[chain]][1],1],accepts[[chain]][kp[[chain]][2],2]),
                                   tune=tune[[chain]][iter-1],
                                   SigmaPriorPars=SigmaPriorPars,
                                   switchllik=switchllik,
                                   tune_q1=tune_q1[[chain]][iter-1],post_sd=post_sd[[chain]][iter-1])
   
   
   #update mu from the STStep_mu
   (PT_chain[[chain]][iter,]                    = out$theta)
  # beta[iter,]                 = c(out1$theta[2],1)
   mu_last                                      = PT_chain[[chain]][iter-1,1]
   accepts[[chain]][kp[[chain]][1],1]           = out$accepts[1]
   accepts[[chain]][kp[[chain]][2],2]           = out$accepts[2]
   
   
  (out1        = STstep_tau(th=PT_chain[[chain]][iter,],n=n,k=k,
                             y=y,
                             acc=accepts1[kptau],
                             SigmaPriorPars=SigmaPriorPars,
                             switchllik=switchllik))
   
   (PT_chain[[chain]][iter,2]  = out1$theta[2])
   
   
  # beta[iter,]             = c(out1$theta[2],1)
   accepts1[kptau]            = out1$accepts
   log_r_bot1[iter]        = out1$log_r_bot
   tau_prop[iter]          = out1$tau_prop
   log_r_top1[iter]        = out1$log_r_top1

 
   
   
   
   
   chain=2
   out          =      STstep_pars(th=PT_chain[[chain]][iter-1,],n=n,p_mu=p_mu,k=k,
                                   mean_mu=PT_chain[[chain]][iter-1,1],
                                   y=y,log_r_bot=log_r_bot[iter-1,],
                                   acc=c(accepts[[chain]][kp[[chain]][1],1],accepts[[chain]][kp[[chain]][2],2]),
                                   tune=tune[[chain]][iter-1],
                                   SigmaPriorPars=SigmaPriorPars,
                                   switchllik=switchllik,
                                   tune_q1=tune_q1[[chain]][iter-1],post_sd=post_sd[[chain]][iter-1])
   
   
   #update mu from the STStep_mu
   (PT_chain[[chain]][iter,]                    = out$theta)
   mu_last                                      = PT_chain[[chain]][iter-1,1]
   accepts[[chain]][kp[[chain]][1],1]           = out$accepts[1]
   accepts[[chain]][kp[[chain]][2],2]           = out$accepts[2]
 
   
}
chain=1   
mllik[[chain]][iter]    = loglik(x=y,mean_mu=PT_chain[[chain]][iter,1],sigma=PT_chain[[chain]][iter,3],tau=PT_chain[[chain]][iter,2],log=T,switchllik=switchllik)$mllik
chain=2
mllik[[chain]][iter]    = loglik(x=y,mean_mu=PT_chain[[chain]][iter,1],sigma=PT_chain[[chain]][iter,3],tau=PT_chain[[chain]][iter,2],log=T,switchllik=switchllik)$mllik



for (chain in (1:npar)){
  #post_sd[[chain]][iter]  = post_sd[[chain]][iter-1]
  
  #tune[[chain]][iter]     = tune[[chain]][iter-1]
 
  
  ##tune transition of alpha 
  #if (iter== kp[[chain]][1]*niter/nstops){  
    
       
  #  post_sd[[chain]][iter]    = var(PT_chain[[chain]][(1:iter),1])
  #  acc_rate[[chain]][kp[[chain]][1],1] = (accepts[[chain]][kp[[chain]][1],1]+1)/(2+(niter/nstops));
  #  kp[[chain]][1]                            =  kp[[chain]][1]+1;
  #  if(iter< niter/2){
      
  #    if(acc_rate[[chain]][kp[[chain]][1]-1,1]>.49 || acc_rate[[chain]][kp[[chain]][1]-1,1]<.39 ){
  #      # adjust the transition density if we are in the first half of the iterations. 
  #      tune[[chain]][iter] = tune[[chain]][iter]*acc_rate[[chain]][kp[[chain]][1]-1,1]/.44
  #    }
  #    
  #  }
  #}

  tune_q1[[chain]][iter]  = tune_q1[[chain]][iter-1]
  
  if (iter== kp[[chain]][2]*niter/nstops){  
    
    acc_rate[[chain]][kp[[chain]][2],2] = (accepts[[chain]][kp[[chain]][2],2]+1)/(2+(niter/nstops));
    kp[[chain]][2]                            =  kp[[chain]][2]+1;
    if(iter< niter/2){
      
      if(acc_rate[[chain]][kp[[chain]][2]-1,2]>.49 || acc_rate[[chain]][kp[[chain]][2]-1,2]<.39 ){
        # adjust the transition density if we are in the first half of the iterations. 
        tune_q1[[chain]][iter] = tune_q1[[chain]][iter]*acc_rate[[chain]][kp[[chain]][2]-1,2]/.44
      }
      
    }
  }

}

   
     
     if ((iter %% 1000) == 0) {
        out_ls=list(PT_chain=PT_chain,rate=acc_rate,accepts=accepts,accepts1=accepts1,log_r_bot=log_r_bot,log_r_bot1=log_r_bot1,
                   kp=kp,y=y,tau_prop=tau_prop,log_r_top=log_r_top1,mllik=mllik)
       
  	     save(out_ls, file=paste("ST",iter,"_",switchllik,".RData",sep=""))
        
        png(paste("DensityPlots_",iter,'_',switchllik,".png",sep=''))
        par(mfrow=c(3,2))
        
        plot(out_ls$PT_chain[[1]][,c(1,2)],main=paste("Joint posterior ",iter," ",switchllik,sep=''))
        plot(out_ls$PT_chain[[1]][,c(1,2)], type='l',main=paste("Joint posterior ",iter," ",switchllik,sep=''))
        
        plot(density(out_ls$PT_chain[[1]][which(!is.na(out_ls$PT_chain[[1]][,1])),2]),main=paste("Density marginal tau "," ",switchllik,sep=''))
        plot(density(out_ls$PT_chain[[1]][which(!is.na(out_ls$PT_chain[[1]][,1])),1]),main=paste("Density marginal mu "," ",switchllik,sep=''))
        plot(density(out_ls$PT_chain[[1]][which(!is.na(out_ls$PT_chain[[1]][,1])),3]),main=paste("Density marginal sigma2 "," ",switchllik,sep=''))
        
        #plot(density(out_ls$mu1[which(!out_ls$mu1==0)]),main="Marginal mu1")
        if (length(which(out_ls$PT_chain[[1]][,2]>0.96))>2){
          plot(density(out_ls$PT_chain[[1]][,1][which(out_ls$PT_chain[[1]][,2]>0.96)]), main=paste("Density marginal mu for tau>0.95 "," ",switchllik,sep=''))
        }
        
        dev.off()
        
        
        png(paste("TracePlots_",iter,'_',switchllik,".png",sep=''))
        par(mfrow=c(3,2))
        
        library(coda)
        par(mfrow=c(3,2))
        traceplot(as.mcmc(out_ls$PT_chain[[1]][which(!is.na(out_ls$PT_chain[[1]][,1])),1]),main=paste("Trace plot mu "," ",switchllik,sep=''))
        traceplot(as.mcmc(out_ls$PT_chain[[1]][which(!is.na(out_ls$PT_chain[[1]][,1])),2]),main=paste("Trace plot tau "," ",switchllik,sep=''))
        traceplot(as.mcmc(out_ls$log_r_bot[which(!(out_ls$log_r_bot==0))]),main=paste("Trace plot log_r_bot mu "," ",switchllik,sep=''))
        traceplot(as.mcmc(out_ls$log_r_bot1[which(!(out_ls$log_r_bot1==0))]),main=paste("Trace plot log_r_bot1 (mu,tau) "," ",switchllik,sep=''))
        traceplot(as.mcmc(out_ls$PT_chain[[1]][which(!is.na(out_ls$PT_chain[[1]][,1])),3]),main=paste("Trace plot sigma2 "," ",switchllik,sep=''))
         
        dev.off()
        
  
     }
    
  }
  
  return(list(PT_chain=PT_chain,rate=acc_rate,accepts=accepts,accepts1=accepts1,log_r_bot=log_r_bot,log_r_bot1=log_r_bot1,
              kp=kp,y=y,tau_prop=tau_prop,log_r_top=log_r_top1,mllik=mllik,swappers=swappers,tune=tune,tune_q1=tune_q1))
}
######################################################################################






#########################################################################################
#runRaftery_Lewis function performs Raftery_Lewis diagnostic test for thre 
#chains. This function returns mcmc objects returned from Raftery_Lewis.
#input: out_list - output set obtained from the samplers
#       q   - raftery.diag argument: the quantile to be estimated
#       t   - raftery.diag argument: the desired margin of error of the estimate
#       eps - raftery.diag argument: precision required for estimate of time to 
#             convergence
#       CI  - raftery.diag argument: the probability of obtaining an estimate in the 
#             interval (q-t,q+t).
#       p   - is switching paramter that alows treating separatelly
#             extraction of the data from different
#             outputs from JAGS sampler and Our sampler  
#       l   - lower bound for number of iterations(allows for performing the test on 
#             burned chain)
#       u   - upper bound for number of iterations(allows for performing the test on 
#             burned chain)
#output: a list of: 
#       mc1 - an mcmc object 
#       mc2 - an mcmc object 
#       rl1$resmatrix - raftery.diag 3-d array result containing M the length of 
#                       "burn in", N the required sample size, Nmin the minimum 
#                       sample size based on zero autocorrelation and I = (M+N)/Nmin 
#                       the "dependence factor" for the house effect1
#       rl2$resmatrix 
#       rl3$resmatrix  

#########################################################################################
runRaftery_Lewis = function(out_list,q = c(0.025),t = 0.01,eps = .01,CI = .95,
                            l=1,u=u)
{
  
  
  #extract chains the sampler output  
  #and turn them into mcmc objects
  mc1=as.mcmc(out_list[[1]][,1][l:u])
  mc2=as.mcmc(out_list[[1]][,2][l:u])
  
  
  #run raftery.diag for each chain
  rl1=raftery.diag(mc1, q=q, r=t, 
                   s=CI, converge.eps=eps)
  
  
  rl2=raftery.diag(mc2, q=q, r=t, 
                   s=CI, converge.eps=eps)
  
  
  #prepare a return list
  return(list(mc1,mc2,rl1$resmatrix,rl2$resmatrix))
}
#########################################################################################



#########################################################################################
#Gewake_diagnostic runs Gewake diagnostic 
#returns matrix of mcmc objects of all chains as well as results from Gewake
#in separate list objects
#input:          out_list  - output set from samplers  
#                l         - lower bound for number of iterations(allows for performing  
#                            the test on burned chain)
#                u         - upper bound for number of iterations(allows for performing  
#                            the test on burned chain)
#                Party     - enables run for separate parties (1-partyI,2-partyII)
#output:         he        - matrix of mcmc objects 
#                he.gew    - list of results from geweke.diag
#########################################################################################

Gewake_diagnostic= function(out_list,l=1,u=u){
  
  
  chains=cbind(out_list[[1]][,1][l:u],out_list[[1]][,2][l:u])
  
  #give names to H.E. matrix columns
  
  colnames(chains)=c("mu","tau")
  
  
  chains.geweke <- apply(chains,2,geweke.diag,0.2,0.5)
  chains.mcmc <- as.mcmc(chains)
  
  return(list(mat=chains.mcmc,gew=chains.geweke))
  
  
}


surf.colors <- function(x, col = terrain.colors(20)) {
  # First we drop the 'borders' and average the facet corners
  # we need (nx - 1)(ny - 1) facet colours!
  x.avg <- (x[-1, -1] + x[-1, -(ncol(x) - 1)] +
              x[-(nrow(x) -1), -1] + x[-(nrow(x) -1), -(ncol(x) - 1)]) / 4
  # Now we construct the actual colours matrix
  colors = col[cut(x.avg, breaks = length(col), include.lowest = T)]
  return(colors)
}

