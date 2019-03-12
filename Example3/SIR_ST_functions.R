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


###################################################################################
# SIR_ST_functions.R
# contains all the functions needed for the algorithm to run on the SIR example
###################################################################################




###################################################################################
#Rhat - solution to the ODE 
#We modeled the data which corresponds to number of death people until time t
#we take the column of the ODE solution that corresponds to R (removed)
#input:     par - a vector of values (beta,alpha,S0,I0)
#output:    a matrix with 2 columns that contains solution for (R,I)
###################################################################################
Rhat = function(par,times){
  pars = par[1:2]
  I0=par[4];
  names(pars) = c("beta","alpha");
  N = 261
  x0 = c(N-I0,I0,0)
  names(x0) = varnames = c("S","I","R")
  SIR = make.SIR()

  sol  = rk4(x0,times,SIR$fn.ode,pars)

  sol  = sol[,2:4]
  fits = sol[,2:3]
  return(fits)
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
#         mllik_eval  - TRUE/FALSE whether to evaluate untempered likelihood
#                       needed for marginal likelihood calculation
#                       by default is FALSE and only lok likelihood is evaluated
# output:   a list 
#             out    - tempered likelihood
#             mllik  - untempered likelihood
#             R1     - solution of DE at pars 
#############################################################
loglik = function(x,pars,tau,parAdd,log=T, mllik_eval=FALSE) {  

	R1      = Rhat(pars,parAdd$times)
	if (length((which(is.nan(R1[,1]))))>0) return(NA)
	
	I       = R1[,1]
	R       = R1[,2]
	N       = parAdd$N
	
	n       = length(I)
        llik    = tau*(sum(dbinom(x,N,R/N,log=log))+dbinom(0,N,I[n]/N,log=log)+dbinom(1,N,I[n-1]/N,log=log))
        mllik    = NULL
        if (mllik_eval){
             mllik   = sum(dbinom(x,N,R/N,log=log))+dbinom(0,N,I[n]/N,log=log)+dbinom(1,N,I[n-1]/N,log=log)
        }
  return(list(out=llik,mllik=mllik,R1=R1))
}
###################################################################################


###################################################################################
#logprior  evaluate log prior
#input:        pars      - a vector of (beta,alpha, S0,I0)
#              PriorPars - prior parameters for alpha, beta and I0
#output:       evaluated log prior     
###################################################################################
logprior=function(pars,PriorPars){

  term12   = sum(dgamma(pars[1:2],PriorPars[1],PriorPars[2],log=T))
  term3    = dbinom(pars[4],PriorPars[3],PriorPars[4],log=T)
  out      = term12+term3
  return(out)
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
#         discrete    - used for optimization, passes the value of I0 
# output: a value of the posterior evaluation
########################################################################
posterior_notau = function(y,pars,tau,PriorPars,parAdd,discrete=NULL,log=T) {
	
		if (!(is.null(discrete))){
                  N         = parAdd$N
  		  pars      = c(pars,N-discrete,discrete)
  		  R=Rhat(pars, parAdd$times)
  		  if (length((which(is.nan(R[,1]))))>0) return(NA)
		}

	      loglik     = loglik(x=y,pars=pars,tau=tau,parAdd=parAdd)$out

        logpriors  = logprior(pars=pars,PriorPars=PriorPars)
        out        = loglik+logpriors
 

  return(list(output=out))
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
posterior =function(pars,y,tau,PriorPars,par_max,parAdd){
  
  llik      = loglik(x=y,pars=pars,tau=tau,parAdd=parAdd)$out
  logpriors = logprior(pars=pars,PriorPars=PriorPars)
  pr_tau    = prior_tau(par_max=par_max,y=y,tau=tau,PriorPars=PriorPars,parAdd=parAdd)
  out       = llik+logpriors+pr_tau
  return(out)
}
#######################################################################
# prior_tau: calculate the prior of tau
# input:  
#         y           - data
#         max         - vector of parameters of interest
#         tau         - a scalar value for tau
#         PriorPars   - priors for all of the parameters  
#         log         - whether log of the prior is evaluated, 
#                       default is true
#         parAdd      - a list of additional parameters, can be either sampled 
#                       parameteres or additional parameters
# output: a value of calculated prior of tau
########################################################################
prior_tau=function(par_max,y,tau,PriorPars,parAdd){
	loglik    = loglik(x=y,pars=par_max,tau=tau,parAdd = parAdd)$out
	logpriors = logprior(pars=par_max,PriorPars=PriorPars)
  return(-loglik-logpriors)
}


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
STstep_pars = function(pars,y,tau,PriorPars,log_r_bot=NA,
                       acc,tune_pars,ttune_pars=NULL,parAdd){

    N                 = parAdd$N
    #evaluate  log_r_bot at current theta
 
    log_r_bot         = posterior_notau(pars=pars,y=y,tau=tau,PriorPars=PriorPars,parAdd = parAdd)$output

    #sample I0
    I0_prop           = rbinom(1,N,pars[4]/N)
 
    pars_prop      = pars
    pars_prop[4]   = I0_prop
    pars_prop[3]   = N-I0_prop
    R              = Rhat(pars_prop,parAdd$times)
    if (length(R[R<0])>0| length(R[R>N]>0)) log_r_top=-Inf else log_r_top=posterior_notau(pars=pars_prop,y=y,tau=tau,PriorPars=PriorPars,parAdd = parAdd)$output
    ratio    = log_r_top-log_r_bot +dbinom(pars[4],N,I0_prop/N,log=T)-dbinom(I0_prop,N,pars[4]/N,log=T)
    prob     = min(1,exp(ratio))
    u        = runif(1)
    (u<=prob)
    if (u<=prob) {
      pars        = pars_prop
      log_r_bot   = log_r_top
      acc[3]      = acc[3]+1
      acc[4]      = acc[4]+1
    }
  
    #sample beta with log normal proposal
    delta           = rnorm(1,0,tune_pars[1])
    logbeta_prop    = log(pars[1])+delta
    (beta_prop      = exp(logbeta_prop))
    pars_prop       = pars
    pars_prop[1]    = beta_prop
    R               = Rhat(pars_prop,parAdd$times)
    if (length(R[R<0])>0 | length(R[R>N])>0) log_r_top = -Inf else log_r_top = posterior_notau(pars=pars_prop,y=y,tau=tau,PriorPars,parAdd=parAdd)$output
    ratio = log_r_top -log_r_bot+dlnorm(pars[1],logbeta_prop,tune_pars[1],log=T)-dlnorm(beta_prop,log(pars[1]),tune_pars[1],log=T)
    prob  = min(1,exp(ratio))
    u     = runif(1)
    (u<=prob)
    if (u<=prob) {
      pars        = pars_prop
      log_r_bot   = log_r_top
      acc[1]      = acc[1]+1
    }
  
    #sample alpha with log normal proposal
    delta           = rnorm(1,0,tune_pars[2])
    logalpha_prop   = log(pars[2])+delta
    (alpha_prop     = exp(logalpha_prop))
    pars_prop       = pars
    pars_prop[2]    = alpha_prop

    R               = Rhat(pars_prop,parAdd$times)
    if (length(R[R<0])>0 | length(R[R>N])>0) log_r_top = -Inf else log_r_top = posterior_notau(pars=pars_prop,y=y,tau=tau,PriorPars=PriorPars,parAdd=parAdd)$output
    ratio = log_r_top-log_r_bot+dlnorm(pars[2],logalpha_prop,tune_pars[2],log=T)-dlnorm(alpha_prop,log(pars[2]),tune_pars[2],log=T)
    prob  = min(1,exp(ratio))
    u     = runif(1)
 
    if (u<=prob) {
      pars            = pars_prop
      log_r_bot       = log_r_top
      acc[2]          = acc[2]+1
    }

  return(list(theta=c(pars,tau),log_r_bot=log_r_bot,accepts=acc,parAdd=parAdd))

}

####################################################################################################
# A wrapper function used for optimization called from the STstep_tau
#         optimizes alpha and beta for 8 different fixed values for I0
#         at proposed and current value of tau
#         this function is called 8 times in parallel
# input:  discrete     - when not NULL used in the optimization routine
#                        to find the solution at each function evaluation of SIRlogpost 
#         pars         - a vector of current values (beta,alpha, S0,I0)
#         y            - data
#         tau_prop     - proposed value of tau
#         te           - current value of tau
#         PriorPars    - prior parameters for alpha, beta and I0
#         parAdd       - a list of additional parameters, can be either sampled 
#                        parameteres or additional parameters
# output: max_prop         - a vector of maximized parameters at proposed tau
#         max_it           - a vector of maximized parameters at current tau
#         eval_fun_p       - maximum value of joint posterior (beta,alpha, I0) at proposed tau
#         eval_fun_i       - maximum value of joint posterior (beta,alpha, I0) at current tau
####################################################################################################
wrapp_fun=function(discrete,pars,y,tau_prop,te,PriorPars,parAdd){
 
 
 if ((is.na(posterior_notau(pars=pars[1:2],y=y,tau=tau_prop,PriorPars=PriorPars,parAdd=parAdd,discrete=discrete)$output)) || (is.na(posterior_notau(pars=pars[1:2],y=y,tau=te,PriorPars=PriorPars,parAdd=parAdd,discrete=discrete)$output)) )   return(c(-19999999999,-19999999999,-19999999999,-19999999999))

 
  max_pars_p_optim=optim(pars[1:2],function(x) (-1)*posterior_notau(pars=c(x[1],x[2]),y=y,tau=tau_prop,
                                                          PriorPars=PriorPars,parAdd=parAdd,discrete=discrete)$output,control=list(maxit=100000000))
 
 
  max_pars_p       = max_pars_p_optim$par

  eval_fun_p       =(-1)*max_pars_p_optim$value
  
  max_pars_i_optim=optim(pars[1:2],function(x) (-1)*posterior_notau(pars=c(x[1],x[2]),y=y,tau=te,
                                                          PriorPars=PriorPars,parAdd=parAdd,discrete=discrete)$output,control=list(maxit=100000000,trace=TRUE))
 

  max_pars_i       = max_pars_i_optim$par
  eval_fun_i       = (-1)*max_pars_i_optim$value

  

  return(list(max_prop=max_pars_p,max_it=max_pars_i,eval_fun_p=eval_fun_p,eval_fun_i=eval_fun_i))

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
##################################################


OptimizePars<- function(y=y,tau_prop=tau_prop,pars=c(pars,tau),PriorPars=PriorPars,parAdd=parAdd, cl){
	
	max_pars_p = max_pars_i=matrix(NA,8,2)
  eval_fun_p = eval_fun_i= rep(NA,8)
	te         = pars[length(pars)]
	pars       = pars[(1:(length(pars)-1))]
	N          = parAdd$N
	times      = parAdd$times


	
	clusterExport(cl, c("wrapp_fun",'Rhat',"make.SIR",'posterior','posterior_notau','times','loglik','logprior','rk4','N'))
	
	getlist=parLapplyLB(cl , 1:8,wrapp_fun,pars=pars,y=y,tau_prop=tau_prop,te=te,PriorPars=PriorPars,parAdd=parAdd)
	

	for (i in (1:8)){
		max_pars_p[i,] = getlist[[i]]$max_prop
	        eval_fun_p[i]  = getlist[[i]]$eval_fun_p
		max_pars_i[i,] = getlist[[i]]$max_it
	        eval_fun_i[i]  = getlist[[i]]$eval_fun_i
	}
	
	# to find the maximum value of I0, just find the I0 that corresponds to the maximum of the 
	# evaluated joint posterior posterior_notau
	#Do this for proposed and for current tau separatelly
	ind=which(eval_fun_p==max(eval_fun_p))
	# max_pars_prop keeps maximized (beta,alpha,S0,I0) at proposed tau
	max_prop=c(max_pars_p[ind,],N-ind,ind)
	
	
	ind1=which(eval_fun_i==max(eval_fun_i))
	# max_pars_it keeps maximized (beta,alpha,S0,I0) at current tau
	max_it=c(max_pars_i[ind1,],N-ind1,ind1)
	
	return(list(max_prop=max_prop, max_it=max_it))
	
	
}


make.SIR<-function(){
  SIR.fun <- function(times, y, p, more) {
    r = y
    r[, "S"] = -p["beta"] * y[, "S"]* y[, "I"]
    r[, "I"] =  p["beta"] * y[, "S"]* y[, "I"]-p["alpha"]*y[,"I"]
    r[, "R"] =                                 p["alpha"]*y[,"I"]        
    return(r)
  }
  SIR.fun.ode <- function(times, y, p) {
    r = y
    
    dimnames(r) = dimnames(y)
    
    r[ "S"] = -p["beta"] * y[ "S"]* y[ "I"]        
    r[ "I"] =  p["beta"] * y[ "S"]* y[ "I"]-p["alpha"]*y["I"]
    r[ "R"] =                               p["alpha"]*y["I"] 
    return(list(r))
  }
  SIR.dfdx <- function(times, y, p, more) {
    #print("SIR.dfdx")
    r = array(0, c(dim(y), dim(y)[2]))
    dimnames(r) = list(NULL, colnames(y), colnames(y))
    r[, "S", "S"] = -p["beta"] * y[, "I"]
    r[, "S", "I"] = -p["beta"] * y[, "S"]
    r[, "I","S"]  =  p["beta"] * y[, "I"]
    r[, "I","I"]  =  p["beta"] * y[, "S"]-p["alpha"]
    r[, "R","I"]  =                       p["alpha"]
    return(r)
  }
  SIR.dfdp <- function(times, y, p, more) {
    #print("SIR.dfdp")
    r = array(0, c(dim(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), names(p))
    
    r[, "S","beta"]  = - y[, "S"]* y[, "I"]
    r[, "I","beta"]  =   y[, "S"]* y[, "I"]
    r[, "I","alpha"] = - y[,"I"]
    r[, "R","alpha"] =   y[,"I"] 
    return(r)
  }
  SIR.d2fdx2 <- function(times, y, p, more) {
    r = array(0, c(dim(y), dim(y)[2], dim(y)[2]))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    r[, "S", "S","I"] = -p["beta"]
    r[, "S", "I","S"] = -p["beta"]
    r[, "I","S","I"]  =  p["beta"]
    r[, "I","I","S"]  =  p["beta"]
    return(r)
  }
  SIR.d2fdxdp <- function(times, y, p, more) {
    r = array(0, c(dim(y), dim(y)[2], length(p)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
    r[, "S", "S","beta"]  = -y[, "I"]
    r[, "S", "I","beta"]  = -y[, "S"]
    r[, "I", "S","beta"]  =  y[, "I"]
    r[, "I", "I","beta"]  =  y[, "S"]
    r[, "I", "I","alpha"] =  -1 
    r[, "R", "I","alpha"] =   1
    return(r)
  }
  return(list(fn = SIR.fun, fn.ode = SIR.fun.ode, dfdx = SIR.dfdx, 
              dfdp = SIR.dfdp, d2fdx2 = SIR.d2fdx2, d2fdxdp = SIR.d2fdxdp))
}
