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
checkpackages("pryr")
#checkpackages("deSolve")
# This library lets us sample from a dirichlet
library(gtools)
library(MCMCpack)
library(mvtnorm)
library(truncnorm)
library(optimx)
library(pryr)
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


# logprior_sigmasq=function(pars, hyperpars=c(rep(10,3),rep(0.05,2))){
#   
#   
#   sigmasq=pars[length(pars)]
#   out    = dgamma(sigmasq,shape=hyperpars[4],scale=hyperpars[5])
#   return(out)
# }

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
#logpost_tau   evaluate joint log posterior  P(alpha,beta,I0,tau / Y)
#input:        
#output:       evaluated log posterior
###################################################################################

 
logpost_tau =function(y,x,tau,p,M,max_pars,pars,hyperpars=c(rep(10,3),rep(0.05,2))){
  
  loglik_eval = loglik(pars=pars,y=y, x=x,tau=tau,p=p,M=M, hyperpars=hyperpars) 
  logpriors   = logprior(pars=pars,p=p,M=M, hyperpars=hyperpars)
  pr_tau      = prior_tau(y=y,x=x,tau=tau,p=p,M=M,max_pars=max_pars, hyperpars=hyperpars)
  out         = loglik_eval+logpriors+pr_tau
  return(out)
}
###################################################################################
#prior_tau   evaluate the log prior for tau
#input:       
#output:       evaluated log prior of tau
###################################################################################
prior_tau=function(y,x,tau,p,M,max_pars, hyperpars=c(rep(10,3),rep(0.05,2))){
 
  loglik_eval=loglik(pars=max_pars,y=y,x=x,tau=tau,p=p,M=M, hyperpars=hyperpars) 
  logpriors = logprior(pars=max_pars,p=p,M=M, hyperpars=hyperpars)
  return(-loglik_eval-logpriors)
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

####################################################################################################
# A wrapper function used for optimization called from the STstep_tau
#         optimizes theta
#         at proposed and current value of tau
# input:  
#         theta        - a vector of current sampled values 
#         y            - data
#         tau_prop     - proposed value of tau
#         te           - current value of tau
#         hyperpars    - prior parameters for alpha, beta and I0
# output: max_pars_p   - maximized (beta, alpha) at proposed tau
#         eval_fun_p   - maximum value of joint posterior (beta,alpha, I0) at proposed tau
#         max_pars_i   - maximized (beta, alpha) at current tau
#         eval_fun_i   - maximum value of joint posterior (beta,alpha, I0) at current tau
####################################################################################################
wrapp_fun=function(theta,y,x,p,M,tau_prop,te,hyperpars,iter, initi){
  
 
     npars=ncol(initi)

 
 if ((is.na(logpost(pars=theta,y=y,x=x,tau=tau_prop,p=p,M=M,  hyperpars=hyperpars))) | (is.na(logpost(pars=theta,y=y,x=x,tau=te,p=p,M=M,  hyperpars=hyperpars))) )   return(c(-19999999999,-19999999999,-19999999999,-19999999999))

     # cluster = makeCluster(4, type = "SOCK")
     # registerDoParallel(cluster)
     # 
     # clusterCall(cluster,function(x) {source('BNN_ST.R'); source('BNN_ST_functions.R')})
     
     #t1=proc.time()[3]
     res = foreach(index=1:2) %dopar%{
       
       
       if(index == 1){
         #1. old optimizer
         max_pars_p_optim=optim(initi[1,(1:(npars-1))],function(l) (-1)*logpost(pars=l,y=y,x=x,tau=tau_prop,p=p,M=M,  hyperpars=hyperpars,optim=1),control=list(maxit=500))
         return(max_pars_p_optim)
       }else if(index==2){
         
         max_pars_i_optim=optim(initi[2,(1:(npars-1))],function(l) (-1)*logpost(pars=l, y=y,x=x,tau=te,p=p,M=M, hyperpars=hyperpars,optim=1),control=list(maxit=500))
         return(max_pars_i_optim)
       }
     }
     #stopCluster(cluster)
     
     # (t=(as.vector(proc.time())[3]-t1)/60)
     # print(paste('optim:',t,sep=''))
      
 
 
      max_pars_p       = res[[1]]$par
      eval_fun_p       =(-1)*res[[1]]$value
  
  
      max_pars_i       = res[[2]]$par
      eval_fun_i       = (-1)*res[[2]]$value

  

  return(c(max_pars_p,eval_fun_p,max_pars_i,eval_fun_i))

}

####################################################################################################
# STstep_tau    update the inverse temperature parameter while keeping the parameters of interest
#               fixed at their current values
#input:         te            - current value of tau
#               theta         - current value of (beta, alpha, S0,I0)
#               y             - data
#               hyperpars     - prior pars (beta, alpha, I0)
#               tune          - if not NULL, tune the step if tau proposal
#               acc           - number of accepted tau
#               cl            - cluster
#output:        theta         - current value of (beta, alpha, S0,I0)
#               te            - updated value of tau
#               log_r_bot1    - un-normalized log posterior (thea,tau \Y) at current tau
#               accepts1      - updated number of accepted tau
#               tau_prop      - proposed tau
#               log_r_top1    - un-normalized log posterior (thea,tau \Y) at proposed tau
####################################################################################################
STstep_tau = function(te,theta,y,x,p,M,hyperpars,tune=NULL,acc=accepts1[kp],iter, samples){
  
  #propose tau
 # if (is.null(tune)){
 #   tau_prop       = runif(n=1,0,1) 
 # }else{
    tau_prop       = rtruncnorm(1,a=0,b=1,te,tune)
 # }


# optimize (beta, alpha) for 8 different values of I0
# run in parallel collect the results in the following matrices and vectors
#max_pars_p =max_pars_i=matrix(NA,8,2)
#eval_fun_p =eval_fun_i= rep(NA,8)


   #clusterExport(cl, c("wrapp_fun",'logpost','loglik','logprior'))

   #getlist=parLapplyLB(cl , 1:8,wrapp_fun,theta,y=y,tau_prop=tau_prop,te=te,priorpars=priorpars,N=N,iter=iter, samples=samples)

   #samples=cbind(theta,tau)
   npars=ncol(samples)
   initProp=samples[which(abs(samples[,npars]-tau_prop)==min(abs(samples[,npars]-tau_prop),na.rm=T)),]
   initi=samples[which(abs(samples[,npars]-te)==min(abs(samples[,npars]-te),na.rm=T)),]
   initi=rbind(initProp,initi)



  getlist= wrapp_fun(theta=theta,y=y,x=x,p=p,M=M,tau_prop=tau_prop,te=te,hyperpars=hyperpars,iter=iter, initi=initi)

   
  max_pars_prop =getlist[1:(npars-1)]
	eval_fun_prop =getlist[npars]
	max_pars_it=getlist[((npars+1): (2*npars-1) )]
	eval_fun_it=getlist[2*npars]



 
  log_r_top= logpost_tau(y=y,x=x,tau=tau_prop,p=p,M=M, max_pars=max_pars_prop, pars=theta, hyperpars=hyperpars)
  log_r_bot=logpost_tau(y=y,x=x,tau=te,p=p,M=M, max_pars=max_pars_it, pars=theta, hyperpars=hyperpars)
  	
	

  #calculate alpha
 # (alpha = log_r_top - log_r_bot+log(dtruncnorm(te,a=0,b=1,tau_prop,tune))-log(dtruncnorm(tau_prop,a=0,b=1,te,tune)))
  # make a decision
   alpha = log_r_top - log_r_bot+log(dtruncnorm(te,a=0,b=1,tau_prop,tune))-log(dtruncnorm(tau_prop,a=0,b=1,te,tune))
  
  
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



