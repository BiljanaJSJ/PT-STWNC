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
loglik=function(x,mean_mu,sigma=1,tau,log=T, switchllik='power'){
  if (log==T)
  {
    if (switchllik=='power'){
    out=tau*sum(dnorm(x,mean=abs(mean_mu),sd=sqrt(sigma),log=log))
	  }else{
	  out=sum(dnorm(x,mean=abs(mean_mu),sd=sqrt(sigma/tau),log=log))
	  }
	 
    }else{
    if (switchllik=='power'){
    out=(prod(dnorm(x,mean=abs(mean_mu),sd=sqrt(sigma),log=log))^tau)
	  }else{
	  out=(prod(dnorm(x,mean=abs(mean_mu),sd=sqrt(sigma/tau),log=log)))
	  }
	  }
    mllik=sum(dnorm(x,mean=abs(mean_mu),sd=sqrt(sigma),log=T)) 
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
#sample_pars
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
STstep_pars = function(th,y=y,n=n,k=1,switchllik='power',post_sd){
  

  # th is the last accepted parameter value
  # tune=2.38
   #calculate the optimal posterior variance of mu
  # post_sd    = (((th[3])*k^2)/(n*th[2]*k^2+th[3]))
  # post_sd    = sqrt(post_sd*tune/2)
   #post_sd    = sqrt(((tune)^2)*post_sd)
   var_target   = th[3]/(n*th[2]+th[3]) 
   #th is the last accepted parameter value
   post_sd    = 2.4*sqrt(var_target)
   
   #propose mu
   mu_prop=rnorm(n=1,mean=th[1],sd=post_sd)

   #the top part of alpha
   log_r_top_mu=posterior_mu(n=n,y=y,tau=th[2],mu=mu_prop,k=k,sigma=th[3],log=T,switchllik=switchllik)
  
   #the last accepted parameter value
   log_r_bot_mu=posterior_mu(n=n,y=y,tau=th[2],mu=th[1],k=k,sigma=th[3],log=T,switchllik=switchllik)
  
  
   #calculate alpha
   alpha_mu = log_r_top_mu - log_r_bot_mu
   # make a decision
  
   if (all(!is.na(alpha_mu) , runif(1) < exp(alpha_mu))){
    # accept the move
   # acc[1]=acc[1]+1;
    th[1] = mu_prop;
    log_r_bot_mu = log_r_top_mu;
    
   }
    
  # do not do anything since we maintain the last value
  return(list(theta=th,log_r_bot=log_r_bot_mu,mu1prop=mu_prop))
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
#         switchllik     - switch the tempered likelihood from 
#                          power tempered into tau in sigma
# output: a value of evaluated posterior
########################################################################
posterior=function(y,tau,mu,k=1,sigma=1,log=T,switchllik='power'){
  
  llik  = loglik(y,mean_mu=mu,sigma=sigma,tau=tau,log=log,switchllik=switchllik)$out
  pr_mu = dprior_mu(x=mu,k=k,log=log)
 
  if (log==T){
    output=llik+pr_mu
  }else{
    output=llik*pr_mu
  }
  return(output)
}


