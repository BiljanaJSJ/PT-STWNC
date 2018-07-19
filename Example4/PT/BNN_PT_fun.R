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

###################################################################################
# BNN_ST_functions.R
# contains all the functions needed for the algorithm to run on the BNN example
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
#checkpackages("coda")
#checkpackages("deSolve")
# This library lets us sample from a dirichlet
library(gtools)
library(MCMCpack)
library(mvtnorm)
library(truncnorm)
library(optimx)
#library(coda)
#library(deSolve)
###################################################################################
#generateInput  Generates input for the BNN
#input:         y   data vector
#               p   number of input points  
#output:        data frame with dim length(y)*p 
###################################################################################

generateInput = function(y,p){

df=list()
for (i in (1:length(y))){
if (i==1) 
{df[[i]]=c(1,rep(NA,p))
}else if (i<=p) {df[[i]]=c(1,rev(y[(1:(i-1))]), rep(NA,p-i+1))
}else{
df[[i]]=c(1,y[((i-p+3):(i-p))])
}
}

Input=do.call(rbind,df)


return(Input)
}

###################################################################################
#psi    - evaluate activation function
#input  - z
#output - scalar
###################################################################################

psi=function(z){
  
  return((exp(z)-exp(-z))/(exp(z)+exp(-z)))
  
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
###################################################################################
#loglik  evaluate tempered likelihood
#also untempered likelihood is evaluated which can be used to estimate the 
#marginal likelihood
#input:        y data vector
#              tau  temperature
#              p number of input units 
#              M numebr of hidden units
#              lambda    shortcut connections dim(p)
#              beta_j    first layer connections dim(M)
#		 	   gamma_ij  second layer connections dim(M)
#			   sigmasq   sigma squared of the error dim(1)
#              fixed hyperparameters: sigma_lambda=10,sigma_beta=10,sigma_gamma=10,nu=0.05,delta=0.05          
#output:         a vector with evaluated tempered and untempered likelihood 
###################################################################################
loglik = function(pars, y,x,tau,p,M,hyperpars=c(rep(10,3),rep(0.05,2)),optim=NULL ){#sigma_lambda=10,sigma_beta=10,sigma_gamma=10,nu=0.05,delta=0.05) {
  
  lambda=pars[(1:(p+1))]
  beta_j=pars[((p+2):(p+1+M))]
  gamma_ij=matrix(pars[((p+M+2):(M*(p+1)+M+p+1))],ncol=M,byrow=T)
  sigmasq=pars[length(pars)]
  
  n=length(y)
  nr=nrow(x)
  
  #llik    = tau*(-(n/2+hyperpars[4]-1)*log(sigmasq)-(1/sigmasq)*(2*hyperpars[5]+sum( y-t(x)*lambda -sum(beta_j*sapply(1:M, function(l)  {psi(t(x)*gamma_ij[,l])}) ,na.rm=T) ,na.rm=T) ) )
  #llik    = tau*(-(n/2)*log(sigmasq)-(1/2*sigmasq)*(sum( y-t(x)*lambda -sum(beta_j*sapply(1:M, function(l)  {psi(t(x)*gamma_ij[,l])}) ,na.rm=T) ,na.rm=T)^2 ) )
  
  
  
  
  if (is.null(optim)){
    
    llik    = tau*sum(dnorm(y,mean=sapply(1:nr, function(l) {sum(t(x[l,])*lambda,na.rm=T)}) + apply(simplify2array(lapply(1:M, function(l) {sapply(1:nr, function(k) { beta_j[l]*psi( sum(t(x[k,])*gamma_ij[,l],na.rm=T) ) } )  })), 1, sum, na.rm=T)   ,sd=sqrt(sigmasq),log=T), na.rm=T)
  }
  if (!is.null(optim)){
    if (pars[length(pars)]<0){ 
      llik   = -9999999999
    }else{
      #llik    = tau*sum(dnorm(y,mean=sapply(1:nrow(x), function(l) {t(x[l,])*lambda}) + apply(simplify2array(lapply(1:M, function(l) {sapply(1:nrow(x), function(k) { beta_j[l]*psi(t(x[k,])*gamma_ij[,l] ) } )  })), c(1,2), sum, na.rm=T)   ,sd=sqrt(sigmasq),log=T), na.rm=T)
      llik    = tau*sum(dnorm(y,mean=sapply(1:nr, function(l) {sum(t(x[l,])*lambda,na.rm=T)}) + apply(simplify2array(lapply(1:M, function(l) {sapply(1:nr, function(k) { beta_j[l]*psi( sum(t(x[k,])*gamma_ij[,l],na.rm=T) ) } )  })), 1, sum, na.rm=T)   ,sd=sqrt(sigmasq),log=T), na.rm=T)
    }  
  }
  
  return(llik)
}
###################################################################################



###################################################################################
#mloglik  evaluate marginal log likelihood
#input:        y data vector
#              tau  temperature
#              p number of input units 
#              M numebr of hidden units
#              lambda    shortcut connections dim(p)
#              beta_j    first layer connections dim(M)
#		 	   gamma_ij  second layer connections dim(M)
#			   sigmasq   sigma squared of the error dim(1)
#              fixed hyperparameters: sigma_lambda=10,sigma_beta=10,sigma_gamma=10,nu=0.05,delta=0.05          
#output:         a vector with evaluated tempered and untempered likelihood 
###################################################################################
mloglik = function(y,x,p,M,pars, hyperpars=c(rep(10,3),rep(0.05,2)) ){#sigma_lambda=10,sigma_beta=10,sigma_gamma=10,nu=0.05,delta=0.05) {
  
  lambda=pars[(1:(p+1))]
  beta_j=pars[((p+2):(p+1+M))]
  gamma_ij=matrix(pars[((p+M+2):(M*(p+1)+M+p+1))],ncol=p+1,byrow=T)
  sigmasq=pars[length(pars)]
  
  n=length(y)
  nr=nrow(x)
  # x=generateInput(y,p)
  #mllik    = (-(n/2)*log(sigmasq)-(1/2*sigmasq)*(sum( y-t(x)*lambda -sum(beta_j*sapply(1:M, function(l)  {psi(t(x)*gamma_ij[,l])}) ,na.rm=T) ,na.rm=T) ) )
  mllik=  sum(dnorm(y,mean=sapply(1:nr, function(l) {sum(t(x[l,])*lambda,na.rm=T)}) + apply(simplify2array(lapply(1:M, function(l) {sapply(1:nr, function(k) { beta_j[l]*psi( sum(t(x[k,])*gamma_ij[,l],na.rm=T) ) } )  })), 1, sum, na.rm=T)   ,sd=sqrt(sigmasq),log=T), na.rm=T)
  
  
  return(mllik)
}
###################################################################################


###################################################################################
#logprior  evaluate log prior
#input:        lambda    shortcut connections dim(p)
#              beta_j    first layer connections dim(M)
#		 	   gamma_ij  second layer connections dim(M)
#			   sigmasq   sigma squared of the error dim(1)
#              fixed hyperparameters: sigma_lambda=10,sigma_beta=10,sigma_gamma=10,nu=0.05,delta=0.05
#output:       evaluated log prior     
###################################################################################
logprior=function(pars,p,M, hyperpars=c(rep(10,3),rep(0.05,2))){
  
  lambda=pars[(1:(p+1))]
  beta_j=pars[((p+2):(p+1+M))]
  gamma_ij=matrix(pars[((p+M+2):(M*(p+1)+M+p+1))],ncol=p+1,byrow=T)
  sigmasq=pars[length(pars)]
  
  
  # term1= sum(lambda^2/(2*hyperpars[1]^2))
  # term2= sum(beta_j^2/(2*hyperpars[2]^2))
  # term3= sum(gamma_ij^2/(2*hyperpars[3]^2))
  # term4= ((hyperpars[4]-1)*log(sigmasq)-(hyperpars[5]/sigmasq) )
  
  term1=sum(dnorm(lambda, mean=0,sd= hyperpars[1],log=T))
  term2=sum(dnorm(beta_j, mean=0,sd= hyperpars[2],log=T))
  term3=sum(dnorm(gamma_ij, mean=0,sd= hyperpars[3],log=T))
  term4=dgamma(1/sigmasq, shape=hyperpars[4],scale=1/hyperpars[5],log=T)
  
  out    = term1+term2+term3+term4
  return(out)
}

###################################################################################
#logpost     evaluate joint log posterior  
#input:        y data vector
#              tau  temperature
#              p number of input units 
#              M numebr of hidden units
#              lambda    shortcut connections dim(p)
#              beta_j    first layer connections dim(M)
#		 	   gamma_ij  second layer connections dim(M)
#			   sigmasq   sigma squared of the error dim(1)
#              fixed hyperparameters: sigma_lambda=10,sigma_beta=10,sigma_gamma=10,nu=0.05,delta=0.05          
#output:       evaluated log posterior
###################################################################################
logpost = function(pars, y,x,tau,p,M, hyperpars=c(rep(10,3),rep(0.05,2)),optim=NULL) {
  
  
  if (is.null(optim)){
    
    loglik_eval=loglik(pars=pars,y=y,x=x,tau=tau,p=p,M=M, hyperpars=hyperpars) 
    logpriors = logprior(pars,p=p,M=M, hyperpars=hyperpars)
    out       = loglik_eval+logpriors
    
  }else if (!is.null(optim)){
    if (pars[length(pars)]<0) 
    {out=-9999999999}else{
      loglik_eval=loglik(pars=pars,y=y,x=x,tau=tau,p=p,M=M,hyperpars=hyperpars) 
      logpriors = logprior(pars,p=p,M=M, hyperpars=hyperpars)
      out       = loglik_eval+logpriors
      
    }  
  }
  
  return(out)
}

###################################################################################
#SamplePars   sample parmeters of interest while tau is held fixed to its 
#             current value
#input:        y               - data
#              theta           - a vector of current values of (beta,alpha,S0,I0)
#              tau             - inverse temperature
#              hyperpars       - prior hyperparameters
#              log_r_bot       - un-normalized log posterior evaluated at current values                      
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
SamplePars = function(y,x,p,M,theta=theta[,iter-1],tau=tau[iter-1],hyperpars=hyperpars,log_r_bot=NA,
                      accepts=accepts,q1,q2,adaptvars=NULL){
  
  
  #sample sigma_sq from log normal
  nPars=length(theta)
  
  
  
  
  for (i in (1:(nPars-1))){
    #sample lambda, beta_j and gamma_ij from normal
    theta_prop     = theta
    theta_prop[i]  = rnorm(n=1,mean=theta[i],sd=q2)
    
    #the top part of alpha
    log_r_top= logpost( pars=theta_prop,y=y,x=x,tau=tau,p=p,M=M, hyperpars=hyperpars)
    
    #the last accepted parameter value
    log_r_bot=logpost(pars=theta,y=y,x=x,tau=tau,p=p,M=M,  hyperpars=hyperpars)
    
    
    #calculate alpha
    alpha = log_r_top - log_r_bot
    # make a decision
    
    if (all(!is.na(alpha) , runif(1) < exp(alpha))){
      # accept the move
      accepts[i] = accepts[i]+1;
      theta         = theta_prop;
      log_r_bot     = log_r_top;
      
    }
  }
  
  
  delta           = rnorm(1,0,q1)
  logtheta_prop   = log(theta[nPars])+delta
  theta_prop      = c(theta[1:(nPars-1)], exp(logtheta_prop))
  
  log_r_top       = logpost(pars=theta_prop,y=y,x=x,tau=tau,p=p,M=M,  hyperpars=hyperpars)
  
  log_r_bot       = logpost(pars=theta,y=y,x=x,tau=tau,p=p,M=M,  hyperpars=hyperpars)
  
  ratio = sum(log_r_top, -log_r_bot,dlnorm(theta[nPars],logtheta_prop[nPars],q1,log=T),-dlnorm(theta_prop[nPars],log(theta[nPars]),q1,log=T),na.rm=T)
  
  prob  = min(1,exp(ratio))
  u     = runif(1)
  
  if (u<=prob) {
    theta                 = theta_prop
    log_r_bot             = log_r_top
    accepts[36]           = accepts[36]+1
  }
  
  return(list(theta=theta,log_r_bot=log_r_bot,accepts=accepts))
  
}
