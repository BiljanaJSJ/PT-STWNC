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
#install packages if they are not installed
###################################################################################
checkpackages=function(package){
  if (!package %in% installed.packages())
    install.packages(package)
}

checkpackages("gtools")
checkpackages("MCMCpack")
checkpackages("mvtnorm")
checkpackages("truncnorm")
checkpackages("optimx")
checkpackages("coda")
checkpackages("deSolve")
# This library lets us sample from a dirichlet
library(gtools)
library(MCMCpack)
library(mvtnorm)
library(truncnorm)
library(optimx)
library(coda)
library(deSolve)

###################################################################################
#Rhat - solution to the ODE 
#We modeled the data which corresponds to number of death people until time t
#we take the column of the ODE solution that corresponds to R (removed)
#input:     par - a vector of values (beta,alpha,S0,I0)
#output:    a matrix with 2 columns that contains solution for (R,I)
###################################################################################
Rhat = function(par){
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
###################################################################################
#SIRloglik  evaluate tempered likelihood, which is binomial
#also untempered likelihood is evaluated which can be used to estimate the 
#marginal likelihood
#input:          y            - data
#                R            - number of removed obtained from the numerical 
#                               solver of the ODE
#                I            - number of infected obtained from the numerical 
#                               solver of the ODE
#                tau          - inverse temperature
#                N            - population size           
#output:         a vector with evaluated tempered and untempered likelihood 
###################################################################################
SIRloglik = function(y,R,I,tau,N=261) {
  n       = length(I)
  llik    = tau*(sum(dbinom(y,N,R/N,log=T))+dbinom(0,N,I[n]/N,log=TRUE)+dbinom(1,N,I[n-1]/N,log=TRUE))
  mllik   = sum(dbinom(y,N,R/N,log=T))+dbinom(0,N,I[n]/N,log=TRUE)+dbinom(1,N,I[n-1]/N,log=TRUE)
  return(c(llik=llik,mllik=mllik))
}
###################################################################################


###################################################################################
#logprior  evaluate log prior
#input:        pars      - a vector of (beta,alpha, S0,I0)
#              I         - numerical solution for I
#              priorpars - prior parameters for alpha, beta and I0
#output:       evaluated log prior     
###################################################################################
logprior=function(pars,I,priorpars=c(1,1,N,5/N)){
  n      = length(I)
  term12 = sum(dgamma(pars[1:2],priorpars[1],priorpars[2],log=T))
  term3  = dbinom(pars[4],priorpars[3],priorpars[4],log=T)
  out    = term12+term3
  return(out)
}

###################################################################################
#SIRlogpost    evaluate joint log posterior  P(alpha,beta,I0 / Y, tau)
#input:        pars      - a vector of (beta,alpha,S0,I0)
#              y         - data
#              R         - numerical solution for R
#              I         - numerical solution for I
#              tau       - inverse temperature 
#              priorpars - prior parameters for alpha, beta and I0
#              N         - population size
#              discrete  - when not NULL used in the optimization routine
#                          to find the solution at each function evaluation
#output:       evaluated log posterior
###################################################################################
SIRlogpost = function(pars,y,R=NULL,I,tau,priorpars=c(1,1,N,5/N),N,discrete=NULL) {
 
  if (!(is.null(discrete))){
    pars=c(pars,N-discrete,discrete)
    R=Rhat(pars)
    if (length((which(is.nan(R[,1]))))>0) return(NA)
    I=R[,1]
    R=R[,2]
  }
  #if (is.nan(sum(dbinom(y,N,R/N,log=T)))){
 # print(pars)
  #print(R)   
  #}
 

  loglik    = SIRloglik(y=y,R=R,I=I,tau=tau)[1]
  logpriors = logprior(pars,I=I,priorpars=priorpars)
  out       = loglik+logpriors
 

  return(out)
}

###################################################################################
#SIRlogpost_tau   evaluate joint log posterior  P(alpha,beta,I0,tau / Y)
#input:        pars      - a vector of (alpha, beta,S0,I0)
#              y         - data
#              R         - numerical solution for R
#              I         - numerical solution for I
#              tau       - inverse temperature 
#              priorpars - prior parameters for alpha, beta and I0
#              N         - population size
#              max_pars  - maximum values of (beta,alpha,I0)
#output:       evaluated log posterior
###################################################################################
SIRlogpost_tau =function(pars,y,R,I,tau,priorpars=c(1,1,N,5/N),max_pars){
  
  loglik    = SIRloglik(y=y,R=R,I=I,tau=tau)[1]
  logpriors = logprior(pars,I=I,priorpars=priorpars)
  pr_tau    = prior_tau(max_pars,y,tau,priorpars=priorpars)
  out       = loglik+logpriors+pr_tau
  return(out)
}
###################################################################################
#prior_tau   evaluate the log prior for tau
#input:        max_pars  - maximum values of (beta,alpha,I0)
#              y         - data
#              tau       - inverse temperature 
#              priorpars - prior parameters for alpha, beta and I0
#output:       evaluated log prior of tau
###################################################################################
prior_tau=function(max_pars,y,tau,priorpars=c(1,1,N,5/N)){
  R         = Rhat(max_pars)
  loglik    = SIRloglik(y=y,R=R[,2],I=R[,1],tau=tau)[1]
  logpriors = logprior(max_pars,I=R[,1],priorpars=priorpars)
  return(-loglik-logpriors)
}


###################################################################################
#SamplePars   sample parmeters of interest while tau is held fixed to its 
#             current value
#input:        y               - data
#              theta           - a vector of current values of (beta,alpha,S0,I0)
#              tau             - inverse temperature
#              priorpars       - prior parameters for alpha, beta and I0
#              log_r_bot       - un-normalized log posterior evaluated at current values 
#                                of (beta,alpha,S0,I0) 
#              N               - population size
#              accepts         - a vector with counts of number of current accepted values
#                                for (beta,alpha,I0)
#              q1              - proposal st.dev for beta
#              q2              - proposal st.dev for alpha
#output:       theta           - a vector with updated parameters (beta,alpha,S0,I0)
#              log_r_bot       - updated un-normalized log posterior evaluated at last 
#                                accepted values of (beta,alpha,S0,I0)
#              accepts         - updated vector with counts of number of new accepted values
#                                for (beta,alpha,I0)
###################################################################################
SamplePars = function(y,theta=theta[,iter-1],tau=tau[iter-1],priorpars=priorpars,log_r_bot=NA,N,
                      accepts=accepts,q1,q2){


    #evaluate  log_r_bot at current theta
    R                 = Rhat(theta)
    log_r_bot         = SIRlogpost(theta,y,R=R[,2],I=R[,1],tau,priorpars=priorpars,N=N)

    #sample I0
    I0_prop        = rbinom(1,N,theta[4]/N)
    #range          = (theta[4]-1):(theta[4]+1)
    #(I0_prop       = sample(range,1))
    theta_prop     = theta
    theta_prop[4]  = I0_prop
    theta_prop[3]  = N-I0_prop
    R              = Rhat(theta_prop)
    if (length(R[R<0])>0| length(R[R>N]>0)) log_r_top=-Inf else log_r_top=SIRlogpost(theta_prop,y,R=R[,2],I=R[,1],tau,priorpars,N=N)
    ratio    = log_r_top-log_r_bot +dbinom(theta[4],N,I0_prop/N,log=T)-dbinom(I0_prop,N,theta[4]/N,log=T)
    prob     = min(1,exp(ratio))
    u        = runif(1)
    (u<=prob)
    if (u<=prob) {
      theta       = theta_prop
      log_r_bot   = log_r_top
      accepts[3]  = accepts[3]+1
    }
  
    #sample beta with log normal proposal
    delta           = rnorm(1,0,q1)
    logbeta_prop    = log(theta[1])+delta
    (beta_prop      = exp(logbeta_prop))
    theta_prop      = theta
    theta_prop[1]   = beta_prop
    R               = Rhat(theta_prop)
    if (length(R[R<0])>0 | length(R[R>N])>0) log_r_top = -Inf else log_r_top = SIRlogpost(theta_prop,y,R=R[,2],I=R[,1],tau,priorpars,N=N)
    ratio = log_r_top -log_r_bot+dlnorm(theta[1],logbeta_prop,q1,log=T)-dlnorm(beta_prop,log(theta[1]),q1,log=T)
    prob  = min(1,exp(ratio))
    u     = runif(1)
    (u<=prob)
    if (u<=prob) {
      theta       = theta_prop
      log_r_bot   = log_r_top
      accepts[1]  = accepts[1]+1
    }
  
    #sample alpha with log normal proposal
    delta           = rnorm(1,0,q2)
    logalpha_prop   = log(theta[2])+delta
    (alpha_prop     = exp(logalpha_prop))
    theta_prop      = theta
    theta_prop[2]   = alpha_prop

    R               = Rhat(theta_prop)
    if (length(R[R<0])>0 | length(R[R>N])>0) log_r_top = -Inf else log_r_top = SIRlogpost(theta_prop,y,R=R[,2],I=R[,1],tau,priorpars,N=N)
    ratio = log_r_top-log_r_bot+dlnorm(theta[2],logalpha_prop,q2,log=T)-dlnorm(alpha_prop,log(theta[2]),q2,log=T)
    prob  = min(1,exp(ratio))
    u     = runif(1)
 
    if (u<=prob) {
      theta           = theta_prop
      log_r_bot       = log_r_top
      accepts[2]      = accepts[2]+1
    }

  return(list(theta=theta,log_r_bot=log_r_bot,accepts=accepts))

}

####################################################################################################
# A wrapper function used for optimization called from the STstep_tau
#         optimizes alpha and beta for 8 different fixed values for I0
#         at proposed and current value of tau
#         this function is called 8 times in parallel
# input:  discrete     - when not NULL used in the optimization routine
#                        to find the solution at each function evaluation of SIRlogpost 
#         theta        - a vector of current values (beta,alpha, S0,I0)
#         y            - data
#         tau_prop     - proposed value of tau
#         te           - current value of tau
#         priorpars    - prior parameters for alpha, beta and I0
#         N            - population size
# output: max_pars_p   - maximized (beta, alpha) at proposed tau
#         eval_fun_p   - maximum value of joint posterior (beta,alpha, I0) at proposed tau
#         max_pars_i   - maximized (beta, alpha) at current tau
#         eval_fun_i   - maximum value of joint posterior (beta,alpha, I0) at current tau
####################################################################################################
wrapp_fun=function(discrete,theta,y,tau_prop,te,priorpars,N,max,iter){
 
 
 if ((is.na(SIRlogpost(pars=theta[1:2],y=y,tau=tau_prop,priorpars=priorpars,N=N,discrete=discrete))) || (is.na(SIRlogpost(pars=theta[1:2],y=y,tau=te,priorpars=priorpars,N=N,discrete=discrete))) )   return(c(-19999999999,-19999999999,-19999999999,-19999999999))

 
  max_pars_p_optim=optim(theta[1:2],function(x) (-1)*SIRlogpost(pars=c(x[1],x[2]),y=y,tau=tau_prop,
                                                          priorpars=priorpars,N=N,discrete=discrete),control=list(maxit=100000000))
 
 
  max_pars_p       = max_pars_p_optim$par

  eval_fun_p       =(-1)*max_pars_p_optim$value
  
  max_pars_i_optim=optim(theta[1:2],function(x) (-1)*SIRlogpost(pars=c(x[1],x[2]),y=y,tau=te,
                                                          priorpars=priorpars,N=N,discrete=discrete),control=list(maxit=100000000,trace=TRUE))
 

  max_pars_i       = max_pars_i_optim$par
  eval_fun_i       = (-1)*max_pars_i_optim$value

  

  return(c(max_pars_p,eval_fun_p,max_pars_i,eval_fun_i))

}

####################################################################################################
# STstep_tau    update the inverse temperature parameter while keeping the parameters of interest
#               fixed at their current values
#input:         te            - current value of tau
#               theta         - current value of (beta, alpha, S0,I0)
#               y             - data
#               priorpars     - prior pars (beta, alpha, I0)
#               tune          - if not NULL, tune the step if tau proposal
#               N             - population size
#               acc           - number of accepted tau
#               cl            - cluster
#output:        theta         - current value of (beta, alpha, S0,I0)
#               te            - updated value of tau
#               log_r_bot1    - un-normalized log posterior (thea,tau \Y) at current tau
#               accepts1      - updated number of accepted tau
#               tau_prop      - proposed tau
#               log_r_top1    - un-normalized log posterior (thea,tau \Y) at proposed tau
####################################################################################################
STstep_tau = function(te,theta,y,priorpars,tune=NULL,N,acc=accepts1[kp],cl,iter){
  
  #propose tau
  if (is.null(tune)){
    tau_prop       = runif(n=1,0,1) 
  }else{
    tau_prop       = rtruncnorm(1,a=0,b=1,te,tune)
  }


# optimize (beta, alpha) for 8 different values of I0
# run in parallel collect the results in the following matrices and vectors
max_pars_p =max_pars_i=matrix(NA,8,2)
eval_fun_p =eval_fun_i= rep(NA,8)


clusterExport(cl, c("wrapp_fun",'Rhat',"make.SIR",'SIRlogpost','times','SIRloglik','logprior','rk4'))

getlist=parLapplyLB(cl , 1:8,wrapp_fun,theta,y=y,tau_prop=tau_prop,te=te,priorpars=priorpars,N=N,iter=iter)

for (i in (1:8)){
max_pars_p[i,] =getlist[[i]][1:2]
eval_fun_p[i] =getlist[[i]][3]
max_pars_i[i,]=getlist[[i]][4:5]
eval_fun_i[i]=getlist[[i]][6]
}

  # to find the maximum value of I0, just find the I0 that corresponds to the maximum of the 
  # evaluated joint posterior SIRlogpost
  #Do this for proposed and for current tau separatelly
  ind=which(eval_fun_p==max(eval_fun_p))
  # max_pars_prop keeps maximized (beta,alpha,S0,I0) at proposed tau
  max_pars_prop=c(max_pars_p[ind,],N-ind,ind)

  
  ind1=which(eval_fun_i==max(eval_fun_i))
  # max_pars_it keeps maximized (beta,alpha,S0,I0) at current tau
  max_pars_it=c(max_pars_i[ind1,],N-ind1,ind1)

  
  # solve the ODE for current theta=(beta,alpha,S0,I0) needed to evaluate SIRlogpost_tau at proposed and current tau
  #(Note that solver is not run within SIRlogpost_tau, for computational efficiency)
  R             =   Rhat(theta)
  log_r_top     =   SIRlogpost_tau(theta,y=y,R=R[,2],I=R[,1],tau=tau_prop,priorpars=priorpars,max_pars=max_pars_prop)
  log_r_bot     =   SIRlogpost_tau(theta,y=y,R=R[,2],I=R[,1],tau=te,priorpars=priorpars,max_pars=max_pars_it)
    

  #calculate alpha
  (alpha = log_r_top - log_r_bot+log(dtruncnorm(te,a=0,b=1,tau_prop,tune))-log(dtruncnorm(tau_prop,a=0,b=1,te,tune)))
  # make a decision
  
  if (all(!is.na(alpha) , runif(1) < exp(alpha))){
    # accept the move
    acc=acc+1;
    te=tau_prop
    }
  
  return(list(te=te,
              log_r_bot1 =log_r_bot,
              accepts1=acc,
              tau_prop=tau_prop,
              log_r_top1=log_r_top)) 
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
