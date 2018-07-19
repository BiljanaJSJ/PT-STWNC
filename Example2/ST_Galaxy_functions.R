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
# loglike:  Tempered Likelihood is Gaussian
# here the tempered likelihood is evaluated, as well as
# the untempered log likelihood which later
# is used to calculate the log marginal likelihood
# input:  y           - Galaxy data from MASS package
#         Z           - matrix of indicator variable
#         p           - vector of mixture weights
#         means       - vector of means
#         sigma2      - vector of sigma2
#         tau         - value of tau
#         k           - number of components
# output:   a list with two elements: llik - tempered likelihood
#                                     mllik - untempered likelihood
#############################################################
loglike =function(y,Z,p,means,sigma2,tau,k){
  
 
  llik     = dnorm(y,apply(t(Z)*matrix(means,length(y),k,byrow=T),1,sum),
                     apply(t(Z)*matrix(sqrt(sigma2),length(y),k,byrow=T),1,sum),log=TRUE)
  
  tllik=tau*sum(llik)
 
 

 

   output = matrix(NA,length(y),length(p))
   for(kp in (1:length(p))){
   
        if (length(sigma2)==length(p)) 
           {sigma2_t=sigma2[kp]
            }else{
             sigma2_t=sigma2
            }
      
        output[,kp]=(p[kp])*dnorm(y,means[kp],sqrt(sigma2_t),log=F)
   
    }

    mllik=sum(log(apply(output,1,sum)))



  return(list(llik=tllik,mllik=mllik))
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
prior_means=function(means,Meanpriorpars){
  sum(dnorm(means,Meanpriorpars[1],Meanpriorpars[2],log=T))
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
prior_sigma2=function(sigma2,Sigmapriorpars){
  sum((dgamma(1/sigma2,shape=Sigmapriorpars[1],scale=1/Sigmapriorpars[2],log=T)))
}
#######################################################################
# posterior_mean: evaluate the conditional posterior distribution of sigma2
# input:  
#         means          - vector of maximized means
#         y              - Galaxy data from MASS package
#         Z              - matrix of indicator variable
#         p              - vector of maximized mixture weights
#         sigma2         - vector of maximized sigma2
#         tau            - value of tau
#         Meanpriorpars  - prior of means
#         k              - number of components
# output: evaluated conditional posterior of means
########################################################################

posterior_mean=function(means,y,Z,p,sigma2,tau,Meanpriorpars,k){
  llik     = loglike(y=y,Z=Z,p=p,means=means,sigma2=sigma2,tau=tau,k=k)$llik
  pmeans   = prior_means(means=means,Meanpriorpars=Meanpriorpars)
  #p_z     = apply(t(Z)*matrix(log(p),length(y),k,byrow=T),1,sum)
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
#         sigma2         - vector of maximized sigma2
#         y              - Galaxy data from MASS package
#         Z              - matrix of indicator variable
#         p              - vector of maximized mixture weights
#         means          - vector of maximized means
#         tau            - value of tau
#         Sigmapriorpars - prior of sigma2
#         k              - number of components
# output: evaluated conditional posterior of sigma2
########################################################################
posterior_sig=function(sigma2,y,Z,p,means,tau,Sigmapriorpars,k){
  llik     = loglike(y=y,Z=Z,p=p,means=means,sigma2=sigma2,tau=tau,k=k)$llik
  psig     = prior_sigma2(sigma2=sigma2,Sigmapriorpars=Sigmapriorpars)
  #p_z     = apply(t(Z)*matrix(log(p),length(y),k,byrow=T),1,sum)
  return(llik+psig)#+sum(p_z))
}
########################################################################


#######################################################################
# prior_tau: calculate the prior of tau
# input:  x              - Galaxy data from MASS package
#         means          - vector of maximized means
#         sigma2         - vector of maximized sigma2
#         p              - vector of maximized mixture weights
#         Z              - matrix of indicator variable
#         tau            - value of tau
#         Meanpriorpars  - prior of means
#         Sigmapriorpars - prior of sigma2
#         Ppriorpars     - prior of mixture weights
#         k              - number of components
# output: a value of calculated prior of tau
########################################################################
prior_tau=function(y,means,sigma2,p,Z,tau,Meanpriorpars,Sigmapriorpars,Ppriorpars,k){
  llik    = loglike(y=y,Z=Z,p=p,means=means,sigma2=sigma2,tau=tau,k=k)$llik
  pmeans  = prior_means(means=means,Meanpriorpars=Meanpriorpars)
  psig    = prior_sigma2(sigma2=sigma2,Sigmapriorpars=Sigmapriorpars)
  pp      = priorp(p=p,Ppriorpars=Ppriorpars)
  
  p_z     = apply(t(Z)*matrix(log(p),length(y),k,byrow=T),1,sum)
  #Nk         = apply(Z,1,sum)
  n=length(y)
  #output  = (n*tau/2)*log(sum(2*pi*t(Z)*matrix(sigma2,length(y),k,byrow=T)))-llik-pmeans-psig-pp-tau*sum(p_z)
  output  = -llik-pmeans-psig-pp-sum(p_z)
  return(output)
}
########################################################################
posterior_optim=function(y,Z,p,means,sigma2,tau,Meanpriorpars,Sigmapriorpars,Ppriorpars,k){
  llik    = loglike(y=y,Z=Z,p=p,means=means,sigma2=sigma2,tau=tau,k=k)$llik
  pmeans  = prior_means(means=means,Meanpriorpars=Meanpriorpars)
  psig    = prior_sigma2(sigma2=sigma2,Sigmapriorpars=Sigmapriorpars)
  pp      = priorp(p=p,Ppriorpars=Ppriorpars)
  p_z     = sum(apply(t(Z)*matrix(log(p),length(y),k,byrow=T),1,sum))
  output=llik+pmeans+psig+pp+p_z
  return(output)
}


#######################################################################
# posterior: evaluate the posterior distribution
# this function calls the prior of tau function
# so it needs current values of the means, sigma2, p and Z as well
# as maximized values of the parameters to supply the prior_t function
# input:  y              - Galaxy data from MASS package
#         Z              - matrix of indicator variable
#         p              - vector of maximized mixture weights
#         means          - vector of maximized means
#         sigma2         - vector of maximized sigma2
#         tau            - value of tau
#         max_means      - maximized value of means
#         max_sigma2     - maximized value of sigma2
#         max_p          - maximized value of p
#         max_Z          - Z sampled from the maximized values of means,sigma2 and p
#         Meanpriorpars  - prior of means
#         Sigmapriorpars - prior of sigma2
#         Ppriorpars     - prior of mixture weights
#         k              - number of components
# output: a value of the posterior evaluation
########################################################################
posterior=function(y,Z,p,means,sigma2,tau,max_means,max_sigma2,max_p,max_Z,Meanpriorpars,Sigmapriorpars,Ppriorpars,k){
  llik    = loglike(y=y,Z=Z,p=p,means=means,sigma2=sigma2,tau=tau,k=k)$llik
  ptau    = prior_tau(y=y,means=max_means,sigma2=max_sigma2,p=max_p,Z=max_Z,tau=tau,Meanpriorpars=Meanpriorpars,Sigmapriorpars=Sigmapriorpars,Ppriorpars=Ppriorpars,k=k)
  #p_z     = apply(t(Z)*matrix(log(p),length(y),k,byrow=T),1,sum)
  #print("c(llik,ptau,sum(p_z))")
  #print(c(llik,ptau))
  output=llik+ptau
  return(list(output=output,ptau=ptau))#+sum(p_z))
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
sampleMeans = function(y,Z,Nk,sigma2,Meanpriorpars,tau){
  sumY  = Z%*%y
  denom = sigma2+tau*Nk*Meanpriorpars[,2]
  mu    = (tau*sumY*Meanpriorpars[,2]+sigma2*Meanpriorpars[,1])/denom
  sd    = sqrt(sigma2*Meanpriorpars[,2]/denom)
  rnorm(n=length(Nk),mu,sd)
}


#############################################################
#Gibbs: Transition one of STWTDNC, sample parameters with fixed tau
# input:  y              - Galaxy data from MASS package
#         means          - vector of maximized means
#         sigma2         - vector of maximized sigma2
#         p              - vector of maximized mixture weights
#         Meanpriorpars  - prior of means
#         tau            - value of tau
#         Sigmapriorpars - prior of sigma2
#         Ppriorpars     - prior of mixture weights
#         accm           - vector of counts of accepted means
#         accs           - vector of counts of accepted sigma2
#         k              - number of components
#         adps           - tuning parameter for the transitional variance of sigma2
#         adpm           - tuning parameter for the transitional variance of means
#         equalVar       - indicator that switches to the model with equal variances 
#                          or to the model with unequal variance
#         post_sd        - transitional variance of the means
#output:  a list of sampled means, sigma2, p and Z
#############################################################
Gibbs=function(y,means,sigma2,p,Meanpriorpars,tau,Sigmapriorpars,Ppriorpars,accm,accs,k,adps,adpm,equalVar){
  

  fmpij      = mpij(y=y,means=means,sigma2=sigma2,p=p,tau=tau)
  pij        = exp(fmpij)/apply(exp(fmpij),1,sum)
  Z          = sampleZ(pij)

  Nk         = apply(Z,1,sum)

  p          = rdirichlet(1,Ppriorpars+Nk)[1,]
  
#   means_new   = sampleMeans(y=y,Z=Z,Nk=Nk,sigma2=sigma2,Meanpriorpars=Meanpriorpars,tau=tau)
#   
#   
#   # sample the variances
#   SSE        = SSE_f(y=y,Z=Z,means=means_new,k=k)
#   if (length(sigma2)==1){SSE=sum(SSE);Nk=sum(Nk);k=1}
#   sigma2_new    = 1/rgamma(k,shape=Nk*tau/2+Sigmapriorpars[1],scale=1/(SSE*tau/2+Sigmapriorpars[2]))    
  
 
  var_target=sigma2/(Nk*tau +sigma2) 
  sd=2.4*sqrt(var_target)

  
  for (l in (1:k)){
  
  means_prop    = means
  means_prop[l] = rtruncnorm(1,a=5,b=40,means[l],sd[l])

  alpha=posterior_mean(y=y,Z=Z,p=p,means=means_prop,sigma2=sigma2,tau=tau,Meanpriorpars=Meanpriorpars,k=k)-
        posterior_mean(y=y,Z=Z,p=p,means=means,sigma2=sigma2,tau=tau,Meanpriorpars=Meanpriorpars,k=k)
        +log(dtruncnorm(means,a=5,b=40,means_prop,sd[l]))-log(dtruncnorm(means_prop,a=5,b=40,means,sd[l]))
  
  if (all(!is.na(alpha) , runif(1) < exp(alpha))){
    means      = means_prop;
    accm[l]    = accm[l] +1
  }
  }
 

if (is.null(equalVar)){


for (l in (1:k)){

  sigma2_prop            = sigma2
  logsig_prop            = log(sigma2)
  logsig_prop[l]         = rnorm(1,log(sigma2[l]),adps[l])
  sigma2_prop[l]         = exp(logsig_prop[l])
  
  alpha=posterior_sig(y=y,Z=Z,p=p,means=means,sigma2=sigma2_prop,tau=tau,Sigmapriorpars=Sigmapriorpars,k=k)-
        posterior_sig(y=y,Z=Z,p=p,means=means,sigma2=sigma2,tau=tau,Sigmapriorpars=Sigmapriorpars,k=k)
        dlnorm(sigma2,logsig_prop,adps,log=T)-dlnorm(sigma2_prop,log(sigma2),adps,log=T)
  
  if (all(!is.na(alpha) , runif(1) < exp(alpha))){
    sigma2    = sigma2_prop;
    accs[l]   = accs[l]+1
  }
}
}else{
  
  sigma2_prop            = sigma2
  logsig_prop            = log(sigma2)
  logsig_prop            = rnorm(1,log(sigma2),adps)
  sigma2_prop            = exp(logsig_prop)
  
  alpha=posterior_sig(y=y,Z=Z,p=p,means=means,sigma2=sigma2_prop,tau=tau,Sigmapriorpars=Sigmapriorpars,k=k)-
        posterior_sig(y=y,Z=Z,p=p,means=means,sigma2=sigma2,tau=tau,Sigmapriorpars=Sigmapriorpars,k=k)
        +dlnorm(sigma2,logsig_prop,adps,log=T)-dlnorm(sigma2_prop,log(sigma2),adps,log=T)
  
  if (all(!is.na(alpha) , runif(1) < exp(alpha))){
    sigma2    = sigma2_prop;
    accs      = accs +1
  }
}

  
  return(list(means=means,sigma2=sigma2,p=p,Z=Z,Nk=Nk,accm=accm,accs=accs))
}



