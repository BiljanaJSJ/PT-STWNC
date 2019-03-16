####################################################################
#article title:  Parallel Tempering via Simulated Tempering Without 
#                Normalizing Constants
#journal name:   Statistics and Computing
#author names:   Biljana Jonoska Stojkova, 
#                David A. Campbell
#affiliation 
#and e-mail 
#address of the 
#corresponding 
#author:         Department of Statistics 
#                University Of British Columbia
#                b.stojkova@stat.ubc.ca
####################################################################







#######################################################
#prior_mu
# mu is N(mean_mu,k^2) distributed
#input:         x         - point at which the prior is evaluated
#               mean_mu   - mean of mu,
#               k         - standard deviation of mu
#               log       - whether log of the prior is evaluated, 
#                           default is true
#output:        prior of mu evaluated at x
#######################################################

dprior_mu=function(x,mean_mu=0,k=1,log=T){
  return(sum(dnorm(x,mean=mean_mu,sd=k,log=log)))
}



#############################################################
# loglik:  Tempered Likelihood is Gaussian
# here the tempered likelihood is evaluated, as well as
# the untempered log likelihood which later
# is used to calculate the log marginal likelihood
# input:  x           - data
#         pars        - vector of parameters of interest
#         tau         - a scalar value for tau
#         log         - TRUE/FALSE should likelihood be evaluated on a log scale or on the original scale
#         parAdd      - a list of additional parameters, can be either sampled 
#                       parameteres or additional parameters
#         mllik_eval  - TRUE/FALSE whether to evaluate untempered likelihood
#                       needed for marginal likelihood calculation
#                       by default is FALSE and only lok likelihood is evaluated
# output:   a list with two elements: out   - tempered likelihood
#                                     mllik - untempered likelihood
#############################################################
loglik=function(x,pars,tau,log=T,parAdd=NULL, mllik_eval=FALSE){
  
    if (log==T)
    {
      out=tau*sum(dnorm(x,mean=abs(pars[1]),sd=sqrt(pars[2]),log=log))

    }else{
      out=(prod(dnorm(x,mean=abs(pars[1]),sd=sqrt(pars[2]),log=log))^tau)
    }
	  mllik    = NULL
	  if (mllik_eval){
    mllik=sum(dnorm(x,mean=abs(pars[1]),sd=sqrt(pars[2]),log=log)) 
	  }
    return(list(out=out,mllik=mllik))

}


#######################################################
#data - simulate data y from N(abs(mu),sigma^2)
#input:         n - sample size
#              mu - mean of the data
#           sigma - standard deviation of the data
#output:    generated data from N(abs(mu),sigma^2) 
#######################################################
data=function(n=5,mu,sigma=1){
  
  return(rnorm(n=n,mean=abs(mu),sd=sigma))
}
#######################################################


#############################################################
#STstep_pars: Transition one of STWTDNC, sample parameters with fixed tau
# input:  
#           pars      - vector of the sampled parameters 
#           tau       - scalar with value of current tau
#           y         - data
#           log_r_bot - a chain of posterior evaluations for updating the parameters
#           acc       - a vector of counts for the accepted updates for each parameter
#           PriorPars - priors for all of the parameters  
#           tune_q1   - tunning chain for the variance of the transition kernerl of the first parameter
#           tune_q2   - tunning chain for the variance of the transition kernerl of the second parameter
#           parAdd    - a list of additional parameters, can be either sampled 
#                       parameteres or additional parameters

#output:  a list 
#           theta     - a vector of updated parameters, the last element is tau
#           parAdd    - sampled indicator matrix denoting which data point was assigned to which mixture component
#           accepts   - a vector of updated counts for the accepted updates for each parameter
#           log_r_bot - posterior evaluation for updated the parameters
#############################################################
STstep_pars = function(pars,tau,y=y,log_r_bot=NA,acc=accepts,
											 PriorPars,tune_pars,ttune_pars=NULL,parAdd=NULL){
  
   
	 n            = length(y)
   #using the Gelman optimal step for the transition variance
   #the transition step is 2.4*Var_of_the_target_poserior
   var_target   = pars[2]/(n*tau+pars[2]) 
  
   #th is the last accepted parameter value
   post_sd    = 2.4*sqrt(var_target)
  
   #propose mu
   mu_prop=rnorm(n=1,mean=pars[1],sd=post_sd)

   #the top part of alpha
   log_r_top_mu=posterior_mu(y=y,tau=tau,pars=c(mu_prop,pars[2]),PriorPars=PriorPars,log=T)
  
   #the last accepted parameter value
   log_r_bot_mu=posterior_mu(y=y,tau=tau,pars=pars,PriorPars=PriorPars,log=T)
  
  
   #calculate alpha
   alpha_mu = log_r_top_mu - log_r_bot_mu
   # make a decision
  
   if (all(!is.na(alpha_mu) , runif(1) < exp(alpha_mu))){
    # accept the move
    acc[1]  = acc[1]+1;
    pars[1] = mu_prop;
    #log_r_bot_mu = log_r_top_mu;
    
   }
    

   log_r_bot=posterior_notau(y=y,pars=pars,tau=tau,log=T,PriorPars=PriorPars)$output
  # do not do anything since we maintain the last value
  return(list(theta=c(pars,tau),accepts=acc,log_r_bot=log_r_bot))
}
#####################################################################################
#posterior_mu function 
#evaluates from the posterior 
#        distribution of (mu/Y,tau)
#input:  
#         y           - univariate data
#         tau         - value of tau
#         pars        - vector of the sampled parameters, last element contains current tau 
#         PriorPars   - priors for all of the parameters  
#         parAdd      - a list of additional parameters, can be either sampled 
#                       parameteres or additional parameters
#output: evaluated conditional posterior of mu 
#        distribution P(mu / Y,tau), given the input parameters
#####################################################################################

posterior_mu=function(y,tau,pars,PriorPars,log=T,parAdd=NULL){

    llik  = loglik(y,pars=pars,tau=tau,log=log)$out
    pr_mu = dprior_mu(x=pars[1],k=PriorPars[2],log=log)
    
    
    if (log==T){ 
      out_mu=llik+pr_mu
    }else{
      out_mu=llik*pr_mu
    }
  
  return(out_mu)
}


#######################################################################
# prior_tau: calculate the prior of tau
# input:  
#         x           - data
#         pars        - vector of parameters of interest
#         tau         - a scalar value for tau
#         PriorPars   - priors for all of the parameters  
#         log         - whether log of the prior is evaluated, 
#                       default is true
#         parAdd      - a list of additional parameters, can be either sampled 
#                       parameteres or additional parameters
# output: a value of calculated prior of tau
########################################################################
prior_tau=function(x,pars,PriorPars,tau,log=T,parAdd=NULL){
  
  llik=loglik(x,pars=pars,tau=tau,log=log)$out
  pr_mu=dprior_mu(x=pars[1],k=PriorPars[2],log=log)
  
  if (log==T){
    prior_t=-llik-pr_mu
  }else{
    prior_t=1/(llik*pr_mu)
  }
  
  return(sum(prior_t ))
}


#######################################################################
# posterior: evaluate the joint posterior distribution P(mu, sigma2,tau / Y)
# input:  
#         x           - data
#         pars        - vector of parameters of interest
#         tau         - a scalar value for tau
#         log         - whether log of the prior is evaluated, 
#         PriorPars   - priors for all of the parameters  
#                       default is true
#         max         - a vector of maximized parameters
#         parAdd      - a list of additional parameters, can be either sampled 
#                       parameteres or additional parameters
# output: a value of evaluated posterior
########################################################################
posterior=function(y,pars,tau,log=T,PriorPars,par_max,parAdd=NULL){
  #max[2]=pars[2]
  ptau  = prior_tau(y,pars=par_max,tau=tau,PriorPars=PriorPars,log=log)
  llik  = loglik(y,pars=pars,tau=tau,log=log)$out
  pr_mu = dprior_mu(x=pars[1],k=PriorPars[2],log=log)
  if (log==T){
    output=llik+pr_mu+ptau
  }else{
    output=llik*pr_mu*ptau
  }
  return(output)
}
#######################################################################
# posterior_notau: evaluate the joint posterior distribution P(mu, sigma2|y, tau)
# input:  
#         y           - data
#         pars        - vector of parameters of interest
#         tau         - a scalar value for tau
#         log         - whether log of the prior is evaluated, 
#         PriorPars   - priors for all of the parameters  
#                       default is true
#         parAdd      - a list of additional parameters, can be either sampled 
#                       parameteres or additional parameters
# output: a value of evaluated posterior
########################################################################
posterior_notau=function(y,pars,tau,log=T,PriorPars,parAdd=NULL){
	
	llik  = loglik(y,pars=pars,tau=tau,log=log)$out
	pr_mu = dprior_mu(x=pars[1],k=PriorPars[2],log=log)
	if (log==T){
		output=llik+pr_mu
	}else{
		output=llik*pr_mu
	}
	return(list(output=output))
}
#############################################################
# Compute the SSE between the data and the mean
#############################################################
SSEfun = function(y,mu){  
  sum((y-abs(mu))^2)
}
##################################################
# OptimizePars  - optimize parameters of interest 
#                 called from Tstep_tau, returns vectors of
#                 optimized parameters for proposed tau and for current tau
# Input:  y           - data
#         tau_prop    - value of proposed tau
#         pars        - vector of the sampled parameters, last element contains current tau 
#         PriorPars   - priors for all of the parameters  
#         parAdd      - a list of additional parameters, can be either sampled 
#                       parameteres or additional parameters
# Output: list
#         max_prop         - a vector of maximized parameters at proposed tau
#         max_it           - a vector of maximized parameters at current tau
#         parAdd_max_prop  - a matrix of maximized Z (indicator matrix of which data point belongs to which mixture component at proposed tau)
#         parAdd_max_it    - a matrix of maximized Z (indicator matrix of which data point belongs to which mixture component at current tau)
##################################################
OptimizePars <- function(y,tau_prop,pars,PriorPars,parAdd=NULL,cl=NULL){
	
	n=length(y)
	y_bar           = mean(y)
	
	N=10000
	mu_prop=mu_it=rep(NA,N)
	mu_prop[1]=mu_it[1]=pars[1]
	post_prop=post_it=rep(NA,N)
	post_prop[1]=post_it[1]=0
	
	
	
	eps=10^-2
	for (i in (2:N)){
		
	mu_prop[i]      = n*y_bar*((tau_prop*PriorPars[2]^2)/(n*tau_prop*PriorPars[2]^2+pars[2]))
	#maximize mean from the posterior mean at current tau
	mu_it[i]        = n*y_bar*((pars[3]*PriorPars[2]^2)/(n*pars[3]*PriorPars[2]^2+pars[2]))
	
	
	post_prop[i]  =  posterior_notau(y=y,pars=c(mu_prop[i],pars[2]),tau=tau_prop,log=T,PriorPars=PriorPars)$output
	post_it[i]    =  posterior_notau(y=y,pars=c(mu_it[i],pars[2]),tau=pars[3],log=T,PriorPars=PriorPars)$output
	#check if the relative error is smaller than the tolerance eps
	if (!(i==1)){
		if ((is.nan(post_prop[i])) || (is.nan(post_it[i])) || (is.nan(post_prop[i-1])) || (is.nan(post_it[i-1]))) break
		if (( (post_prop[i]-post_prop[i-1])/post_prop[i]<eps) & ( (post_it[i]-post_it[i-1])/post_it[i]<eps)) {
			break}
		
	}
	}
	
	return(list(max_prop=c(mu_prop[i],pars[2]), 
							max_it=c(mu_it[i],pars[2])))
	
}


# ######################################################################################
# #END!!





