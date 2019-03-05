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


############################################################
# surf.colors : Graphical, used to color the surface of the 2D function with 
# terrain colors
# input :    x  - matrix of function surface estimates
#Source: https://stat.ethz.ch/pipermail/r-help/2003-September/039104.html
############################################################
surf.colors <- function(x, col = terrain.colors(20)) {
  # First we drop the 'borders' and average the facet corners
  # we need (nx - 1)(ny - 1) facet colours!
  x.avg <- (x[-1, -1] + x[-1, -(ncol(x) - 1)] +
              x[-(nrow(x) -1), -1] + x[-(nrow(x) -1), -(ncol(x) - 1)]) / 4
  # Now we construct the actual colours matrix
  colors = col[cut(x.avg, breaks = length(col), include.lowest = T)]
  return(colors)
}

#############################################################
# loglik:  Tempered Likelihood is Gaussian
# here the tempered likelihood is evaluated, as well as
# the untempered log likelihood which later
# is used to calculate the log marginal likelihood
# input:  x           - Galaxy data from MASS package
#         pars        - vector of parameters of interest
#         tau         - a scalar value for tau
#         log         - TRUE/FALSE should likelihood be evaluated on a log scale or on the original scale
#         parAdd      - a list of additional parameters, can be either sampled 
#                       parameteres or additional parameters
# output:   a list 
#             out    - tempered likelihood
#             mllik  - untempered likelihood
#             sigma2 - extracted parameters from pars vector
#             p      - extracted parameters from pars vector  
#############################################################
loglik =function(x,pars,tau,log=T,parAdd){
 
	k=pars[length(pars)]
	if (((length(pars)-1) %% k)==0) {
		equalVar=NULL
		means=pars[1:k]
		sigma2=pars[(k+1):(2*k)]
		p=pars[(2*k+1):(2*k+k)]
	}else{
		equalVar=1
    	        means=pars[1:k]
	        sigma2=pars[((k+1):(k+1))]
	        p=pars[(k+2):(2*k+1)]
	}
	
  llik     = dnorm(x,apply(t(parAdd)*matrix(means,length(x),k,byrow=T),1,sum),
                     apply(t(parAdd)*matrix(sqrt(sigma2),length(x),k,byrow=T),1,sum),log=TRUE)
  
  tllik=tau*sum(llik)
 
 
   output = matrix(NA,length(x),length(p))
   for(kp in (1:length(p))){
   
        if (length(sigma2)==length(p)) 
           {sigma2_t=sigma2[kp]
            }else{
             sigma2_t=sigma2
            }
      
        output[,kp]=(p[kp])*dnorm(x,means[kp],sqrt(sigma2_t),log=F)
   
    }

    mllik=sum(log(apply(output,1,sum)))



  return(list(out=tllik,mllik=mllik,sigma2=sigma2,p=p))
}

############################################################


############################################################
#prior_means  - evaluate prior of means
#Gaussian prior
#input:         
#            means         - value at which prior is evaluated
#            Meanpriorpars - mean and st.dev of the prior
#                            defined by the user 
#            log           - whether log of the prior is evaluated, 
#                            default is true
#output:     evaluated prior of means for each of the components
############################################################
prior_means=function(means,Meanpriorpars,log=T){
  sum(dnorm(means,Meanpriorpars[1],Meanpriorpars[2],log=log))
}
############################################################


############################################################
# prior_sigma2:Prior of Sigmas 
# sigmas are distributed IG(priorpars[1],priorpars[2])
# input:     Sigma         - values at which prior of sigma2 is evaluated
#            priorpars     - vector of shape and scale of the prior sigma2
#            log           - whether log of the prior is evaluated, 
#                            default is true
#output:     evaluated prior of sigma2 for each of the components
############################################################
prior_sigma2=function(sigma2,Sigmapriorpars, log=T){
  sum((dgamma(1/sigma2,shape=Sigmapriorpars[1],scale=1/Sigmapriorpars[2],log=log)))
}
#######################################################################
# posterior_mean: evaluate the conditional posterior distribution of sigma2
# input:  
#         y           - data
#         pars        - vector of parameters of interest
#         tau         - a scalar value for tau
#         PriorPars   - priors for all of the parameters  
#         parAdd      - a list of additional parameters, can be either sampled 
#                       parameteres or additional parameters
# output: evaluated conditional posterior of means
########################################################################

posterior_mean=function(y,pars,tau,PriorPars,parAdd){

	k        = pars[length(pars)]
	
	llik     = loglik(x=y,tau=tau,pars=pars,parAdd=parAdd)$out
  pmeans   = prior_means(means=pars[1:k],Meanpriorpars=PriorPars)

  return(llik+pmeans)#+sum(p_z))
}
############################################################


############################################################
# priorp:  Prior of weights (p) of the mixing distributions
# p are distributed Dirichlet(1,..,1)
# input:   p              - values at which prior of p is evaluated
#          Ppriorpars     - vector of concentration parameters 
#                           for the prior of p
############################################################
priorp = function(p,Ppriorpars){
  d=ddirichlet(p,Ppriorpars)
  return(log(d))
}
############################################################




#######################################################################
# posterior_sig: evaluate the conditional posterior distribution of sigma2
# input:  
#         y           - data
#         pars        - vector of parameters of interest
#         tau         - a scalar value for tau
#         PriorPars   - priors for all of the parameters  
#         parAdd      - a list of additional parameters, can be either sampled 
#                       parameteres or additional parameters
# output: evaluated conditional posterior of sigma2
########################################################################
posterior_sig=function(y,pars,tau,PriorPars,parAdd){
		
		llike     = loglik(x=y,pars=pars,tau=tau,parAdd=parAdd)
		llik      = llike$out
	    psig      = prior_sigma2(sigma2=llike$sigma2,Sigmapriorpars=PriorPars)

    return(llik+psig)
}
########################################################################


#######################################################################
# prior_tau: calculate the prior of tau
# input:  
#         y           - data
#         pars        - vector of parameters of interest
#         tau         - a scalar value for tau
#         PriorPars   - priors for all of the parameters  
#         log         - whether log of the prior is evaluated, 
#                       default is true
#         parAdd      - a list of additional parameters, can be either sampled 
#                       parameteres or additional parameters
# output: a value of calculated prior of tau
########################################################################
prior_tau=function(y,pars,tau,PriorPars,log=T,parAdd=NULL){


	k         = pars[length(pars)]	
	llike     = loglik(x=y,pars=pars,tau=tau,parAdd=parAdd)
	llik      = llike$out
	
  pmeans  = prior_means(means=pars[1:k],Meanpriorpars=PriorPars[1:2])
  psig    = prior_sigma2(sigma2=llike$sigma2,Sigmapriorpars=PriorPars[3:4])
  pp      = priorp(p=llike$p,Ppriorpars=PriorPars[(5:(4+k))])
  
  p_z     = apply(t(parAdd)*matrix(log(llike$p),length(y),k,byrow=T),1,sum)
 

  output  = -llik-pmeans-psig-pp-sum(p_z)
  return(output)
}

#######################################################################
# posterior_optim: optimization function called from OptimizePars
# input:  
#         y           - data
#         Z           - matrix of indicator variable
#         p           - vector of mixture weights
#         means       - vector of means
#         sigma2      - vector of sigma2
#         tau         - value of tau
#         PriorPars   - priors for all of the parameters  
#         k           - number of components in the mixture of Gaussians model
# output: evaluated posterior_optim
########################################################################

posterior_optim=function(y,Z,p,means,sigma2,tau,PriorPars,k){
	
	  	
	 pars      = c(means,sigma2,p,k)
	 llike     = loglik(x=y,pars=pars,tau=tau,parAdd=Z)
	 llik      = llike$out


	pmeans  = prior_means(means=means,Meanpriorpars=PriorPars[1:2])
	psig    = prior_sigma2(sigma2=sigma2,Sigmapriorpars=PriorPars[3:4])
	pp      = priorp(p=p,Ppriorpars=PriorPars[5:(4+k)])
	p_z     = sum(apply(t(Z)*matrix(log(p),length(y),k,byrow=T),1,sum))
	output  = llik+pmeans+psig+pp+p_z
	return(output)
}

#######################################################################
# posterior_notau: evaluate the posterior distribution
#                  this function calls the prior of tau function
# input:  y              - Galaxy data from MASS package
#         pars        - vector of parameters of interest
#         tau         - a scalar value for tau
#         PriorPars   - priors for all of the parameters  
#         log         - whether log of the prior is evaluated, 
#                       default is true
#         parAdd      - a list of additional parameters, can be either sampled 
#                       parameteres or additional parameters
# output: a value of the posterior evaluation
########################################################################
posterior_notau=function(y,pars,tau,log=T,PriorPars,parAdd=NULL){
	
	#k         = pars[length(pars)]	
	llike     = loglik(x=y,pars=pars,tau=tau,parAdd=parAdd)
	llik      = llike$out	
	
	output=llik
	return(list(output=output))#+sum(p_z))
}

#######################################################################
# posterior: evaluate the posterior distribution
# this function calls the prior of tau function
# so it needs current values of the means, sigma2, p and Z as well
# as maximized values of the parameters to supply the prior_t function
# input:  y              - Galaxy data from MASS package
#         pars        - vector of parameters of interest
#         tau         - a scalar value for tau
#         log         - whether log of the prior is evaluated, 
#                       default is true
#         PriorPars   - priors for all of the parameters  
#         max         - a vector of optimized paramters of interest
#         parAdd      - a list of additional parameters, can be either sampled 
#                       parameteres or additional parameters
# output: a value of the posterior evaluation
########################################################################
posterior=function(y,pars,tau,log=T,PriorPars,par_max,parAdd=NULL,parAdd_max=NULL){
	

		#k         = pars[length(pars)]	
		llike     = loglik(x=y,pars=pars,tau=tau,parAdd=parAdd)
		llik      = llike$out	
			
		ptau    = prior_tau(y=y,pars=par_max,tau=tau,PriorPars=PriorPars,parAdd=parAdd_max)
        output  = llik+ptau
 
  return(output)
}
#############################################################
# Compute the SSE between the data and the group means
#############################################################
SSE_f = function(y,Z,means,k){  
  apply((t(Z)*y-t(Z)*(matrix(means,length(y),k,byrow=T)))^2,2,sum)
}



#############################################################
#mpij:   calculate full data log Likelihood 
# This is only used to sample Z, the indicator variable
# which is a matrix of (k X niter), indicating each of the
# data points to which mixture component belongs
# input:  y            - Galaxy data from MASS package
#         means        - vector of means
#         sigma2       - vector of sigma2
#         p            - vector of mixture weights
#output: matrix (n X k) containiing probabilities
#        for each of the data point belong to each of the components          
#        which are used as probabilities for the multinomial distribution
#        when sampling the latent variable Z
#############################################################
mpij = function(y,means,sigma2,p,tau){
  output = matrix(NA,length(y),length(p))
  for(kp in (1:length(p))){
     if (length(sigma2)==length(p)) 
         {sigma2_t=sigma2[kp]
         }else{
          sigma2_t=sigma2
         }
      output[,kp]=log(p[kp])+ dnorm(y,means[kp],sqrt(sigma2_t),log=T)
  }
  return((output-apply(output,1,max)))
}
#############################################################

#############################################################
# sampleZ - samples Z as multinomial using probabilities obtained from 
# mpij
# input:        pij - probabilities obrtained from 
#                     mpij
# output:       sampled indicator variable
#############################################################
sampleZ = function(pij){
  apply(pij,1,function(x){rmultinom(n=1,size= 1, x)})
}

#############################################################
# sampleMeans - samples means (this function is used by Gibbs sampler)
# input:        y                    - Galaxy data
#               Z                    - the indicator variable that matches data points to the
#                                      mixture components              
#               Nk                   - a vector with number of data points in each of the components
#               sigma2               - sigma2 last accepted value
#               Meanpriorpars        - prior mean and st.dev of means
#               tau                  - temperature
# output:       sampled means 
#############################################################
# sampleMeans = function(y,Z,Nk,sigma2,Meanpriorpars,tau){
#   sumY  = Z%*%y
#   denom = sigma2+tau*Nk*Meanpriorpars[,2]
#   mu    = (tau*sumY*Meanpriorpars[,2]+sigma2*Meanpriorpars[,1])/denom
#   sd    = sqrt(sigma2*Meanpriorpars[,2]/denom)
#   rnorm(n=length(Nk),mu,sd)
# }


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
STstep_pars =	function(pars,tau,y,log_r_bot=NA,acc,
            					 PriorPars,tune_pars,ttune_pars,parAdd){  
  

	k          = pars[length(pars)]
	if (((length(pars)-1) %% k)==0) {
	 	equalVar = NULL
		means    = pars[1:k]
		sigma2   = pars[(k+1):(2*k)]
		p        = pars[(2*k+1):(2*k+k)]
	}else{
		equalVar = 1
		means    = pars[1:k]
		sigma2   = pars[((k+1):(k+1))]
		p        = pars[(k+2):(2*k+1)]
	}
	pars_tune  = which(ttune_pars)
	if (length(pars_tune)==0){
	if (is.null(equalVar)){
	pars_tune=(k+1):(2*k)
	}else{
	pars_tune=((k+1):(k+1))
    }	
	}
	
  fmpij      = mpij(y=y,means=means,sigma2=sigma2,p=p,tau=tau)
  pij        = exp(fmpij)/apply(exp(fmpij),1,sum)
  Z          = sampleZ(pij)

  Nk         = apply(Z,1,sum)

  
  p          = rdirichlet(1,PriorPars[5:(4+k)]+Nk)[1,]
  
  if (((length(pars)-1) %% k)==0) {
  pars[(2*k+1):(2*k+k)] = p
  }else{
  pars[(k+2):(2*k+1)]   = p	
  }

  var_target=sigma2/(Nk*tau +sigma2) 
  sd_mean=2.4*sqrt(var_target)

  
  for (l in (1:k)){
  
  means_prop    = means
  means_prop[l] = rtruncnorm(1,a=5,b=40,means[l],sd_mean[l])

  pars_prop     = pars
  pars_prop[l]  = means_prop[l]
  
  alpha=posterior_mean(y=y,pars=pars_prop,tau=tau,PriorPars=PriorPars[1:2],parAdd=Z)-
        posterior_mean(y=y,pars=pars,tau=tau,PriorPars=PriorPars[1:2],parAdd=Z)+
        log(dtruncnorm(means,a=5,b=40,means_prop,sd_mean))-log(dtruncnorm(means_prop,a=5,b=40,means,sd_mean))
  
  if (all(!is.na(alpha) , runif(1) < exp(alpha))){
    means      = means_prop;
    acc[l]     = acc[l] +1
    pars       = pars_prop
  }
  }
 

if (is.null(equalVar)){


for (l in (1:k)){

  sigma2_prop            = sigma2
  logsig_prop            = log(sigma2)
  logsig_prop[l]         = rnorm(1,log(sigma2[l]),tune_pars[pars_tune[l]])
  sigma2_prop[l]         = exp(logsig_prop[l])
  
  pars_prop                  = pars
  pars_prop[(((k+1):(2*k)))] = sigma2_prop
  
  alpha=posterior_sig(y=y,pars=pars_prop,tau=tau,PriorPars=PriorPars[3:4],parAdd=Z)-
        posterior_sig(y=y,pars=pars,tau=tau,PriorPars=PriorPars[3:4],parAdd=Z)+
        dlnorm(sigma2,logsig_prop,tune_pars[pars_tune],log=T)-dlnorm(sigma2_prop,log(sigma2),tune_pars[pars_tune],log=T)
  
  
  if (all(!is.na(alpha) , runif(1) < exp(alpha))){
    sigma2              = sigma2_prop;
    acc[k+l]            = acc[k+l]+1
    pars                = pars_prop
  }
}
}else{
  
  sigma2_prop            = sigma2
  logsig_prop            = log(sigma2)
  logsig_prop            = rnorm(1,log(sigma2),tune_pars[pars_tune])
  sigma2_prop            = exp(logsig_prop)
  
  pars_prop              = pars
  pars_prop[k+1]         = sigma2_prop
  
  
  
  alpha=posterior_sig(y=y,pars=pars_prop,tau=tau,PriorPars=PriorPars[3:4],parAdd=Z)-
      	posterior_sig(y=y,pars=pars,tau=tau,PriorPars=PriorPars[3:4],parAdd=Z)+
        dlnorm(sigma2,logsig_prop,tune_pars[pars_tune],log=T)-dlnorm(sigma2_prop,log(sigma2),tune_pars[pars_tune],log=T)
  
  
  
  if (all(!is.na(alpha) , runif(1) < exp(alpha))){
    sigma2                   = sigma2_prop;
    acc[k+1]                 = acc[k+1] +1
    pars                     = pars_prop
  }
}
   # th=c(means,sigma2,p,k)
    log_r_bot=posterior_notau(y=y,pars=pars,tau=tau,log=T,PriorPars=PriorPars,parAdd=Z)$output

    return(list(theta=c(pars,tau),parAdd=Z,accepts=acc,log_r_bot=log_r_bot))
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

OptimizePars <- function(y,tau_prop,pars,PriorPars,parAdd=NULL, cl=NULL){
	
	tau=pars[length(pars)]
	pars=pars[(1:(length(pars)-1))]
	
	
	k=pars[length(pars)]
	if (((length(pars)-1) %% k)==0) {
		equalVar=NULL
		means=pars[1:k]
		sigma2=pars[(k+1):(2*k)]
		p=pars[(2*k+1):(2*k+k)]
	}else{
		equalVar=1
		means=pars[1:k]
		sigma2=pars[((k+1):(k+1))]
		p=pars[(k+2):(2*k+1)]
	}
	Z=parAdd
	
	N=10000
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
		denom_prop              = max_sigma2_prop[i-1,]+tau_prop*Nk_prop[i,]*rep(PriorPars[2],k)
		max_means_prop[i,]      = (tau_prop*sumY*rep(PriorPars[2],k)+max_sigma2_prop[i-1,]*rep(PriorPars[1],k))/denom_prop
		
		#maximize sigma2 at proposed tau
		#sigma2 is optimized using the posterior mode
		SSE                        = SSE_f(y=y,Z_prop_max[[i]],means=max_means_prop[i,],k=k) 
		
		if (is.null(equalVar)){
			d0_prop                = (Nk_prop[i,]*tau_prop/2+PriorPars[3] )
			v0_prop                = (tau_prop*SSE/2+PriorPars[4])
			max_sigma2_prop[i,]    = v0_prop/(d0_prop+1)
		}else{
			SSE                    = sum(SSE)
			d0_prop                = (sum(Nk_prop[i,])*tau_prop/2+PriorPars[3] )
			v0_prop                = tau_prop*SSE/2+PriorPars[4]  
			max_sigma2_prop[i,]    = v0_prop/(d0_prop+1)
		}
		
		
		#maximize means and sigma2 for the current value of tau
		sumY                     = Z_it_max[[i]]%*%y
		denom_it                 = max_sigma2_it[i-1,]+tau*Nk_it[i,]*rep(PriorPars[2],k)
		max_means_it[i,]         = (tau*sumY*rep(PriorPars[2],k)+max_sigma2_it[i-1,]*rep(PriorPars[1],k))/denom_it
		
		
		#optimize sigmas from the posterior mode at current tau
		SSE                      = SSE_f(y,Z_it_max[[i]],means=max_means_it[i,],k=k) 
		
		if (is.null(equalVar)){
			#if the model runs with unequal variances
			# the maximized sigma2 is a vector of length k  
			d0                     = (Nk_it[i,]*tau/2+PriorPars[3])
			v0                     = (tau*SSE/2+PriorPars[4])
			max_sigma2_it[i,]      = v0/(d0+1)
		}else{
			#if the model runs with equal variances 
			# the maximized sigma2 is a one value
			SSE                    = sum(SSE)
			d0                     = (sum(Nk_it[i,])*tau/2+PriorPars[3] )
			v0                     = (tau*SSE/2+PriorPars[4])
			max_sigma2_it[i,]      = v0/(d0+1)
		}
		
		#evaluate the posterior P(mu, sigma2,p,Z | y)  at the maximums for proposed and current tau
		post_prop[i]  =  posterior_optim(y=y,Z=Z_prop_max[[i]],p=p_max_prop[i,],means=max_means_prop[i,],sigma2=max_sigma2_prop[i,],tau=tau_prop,PriorPars=PriorPars,k=k)
		post_it[i]    =  posterior_optim(y=y,Z=Z_it_max[[i]],p=p_max_it[i,],means=max_means_it[i,],sigma2=max_sigma2_it[i,],tau=tau,PriorPars=PriorPars,k=k)
		#check if the relative error is smaller than the tolerance eps
		if (!(i==1)){
			if ((is.nan(post_prop[i])) || (is.nan(post_it[i])) || (is.nan(post_prop[i-1])) || (is.nan(post_it[i-1]))) break
			if (( (post_prop[i]-post_prop[i-1])/post_prop[i]<eps) & ( (post_it[i]-post_it[i-1])/post_it[i]<eps)) {
				break}
			
		}
	}
	
	return(list(max_prop=c(max_means_prop[i,],max_sigma2_prop[i,],p_max_prop[i,],k), 
							max_it=c(max_means_it[i,],max_sigma2_it[i,],p_max_it[i,],k),parAdd_max_prop=Z_prop_max[[i]],parAdd_max_it=Z_it_max[[i]]))
	
}	


