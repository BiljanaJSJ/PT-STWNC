
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

####################
# Functions for the mixture model used in the Galaxy example
###################

# Checking the Availability of packages using a nice function 
checkpackages=function(package){
  if (!package %in% installed.packages())
    install.packages(package)
}

checkpackages("gtools")
checkpackages("MCMCpack")
checkpackages("mvtnorm")

library(gtools)
library(MCMCpack)
library(mvtnorm)


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


##############################################################################
# Log prior for mean
#input:      means            - value of the means at which the prior is evaluated
#            p                - mixture probabilities (last accepted values)
#            Meanpriorpars    - prior parameters mean and st.dev.
#output:     Evaluated prior for the means
##############################################################################
logpriormean = function(means,p,Meanpriorpars){
  # Gaussian Priors
  return( dnorm(means,Meanpriorpars[,1],sqrt(Meanpriorpars[,2]),log=TRUE))
}

##############################################################################
# Log Prior of the mixing probabilities 
#input:      p                - value of p at which the prior os evaluated
#            Ppriorpars       - prior parameters for p Dirichlet(1,,1)
#output:     evaluetd prior of p
##############################################################################
logpriorp = function(p,Ppriorpars,k=3){
  # Dirichlet
  log(ddirichlet(p,Ppriorpars))
}

##############################################################################
# Log prior for sigma2
#input:      Sigma            - value of the sigma2 at which the prior is evaluated
#            priorpars        - shape and scale of the prior of sigma2
#output:     Evaluated prior for the sigma2
##############################################################################
logpriorSigma = function(Sigma,priorpars){
  
      return(dgamma(1/Sigma,shape=priorpars[1],scale = 1/priorpars[2],log=TRUE))    
}


##############################################################################
# Log prior for sigma2
#input:      y                - Galaxy data
#            means            - value of the means at 
#                               which the posterior is evaluated
#            Sigma            - value of the sigma2 at 
#                               which the posterior is evaluated
#            p                - value of p 
#            Z                - the indicator variable that matches data points to the
#                               mixture components
#            beta             - temperature parameter
#            k                - number of components
#            Meanpriorpars    - prior parameters mean and st.dev.
#            priorpars        - shape and scale of the prior of sigma2
#            Ppriorpars       - prior parameters for p Dirichlet(1,,1)  
#output:     Evaluated log joint posterior  (means,sigma2,p,beta)
##############################################################################
posterior=function(y,means,Sigma,p,Z,beta,k=3,Meanpriorpars,priorpars,Ppriorpars){

llik=conditionalloglike(y,means,Sigma,p,Z,beta,k=k)
pmean=logpriormean(means,p,Meanpriorpars)
psigma=logpriorSigma(Sigma,priorpars)
pp=logpriorp(p,Ppriorpars,k=k)

return(sum(llik)+sum(pmean)+sum(psigma)+sum(pp))
}
##############################################################################
#SSEfun: Compute the SSE between the data and the group means
#input:    y         - Galaxy data
#          Z         - the indicator variable that matches data points to the
#                      mixture components
#          means     - value of the means at 
#                      which the posterior is evaluated
#          k         - number of components
##############################################################################
SSEfun = function(y,Z,means,k){  
  apply((t(Z)*y-t(Z)*(matrix(means,length(y),k,byrow=T)))^2,2,sum)
}



#############################################################
#fulldatajointloglike:  calculate full data log Likelihood 
# This is only used to sample Z, the indicator variable
# which is a matrix of (k X niter), indicating each of the
# data points to which mixture component belongs
# input:  y            - Galaxy data from MASS package
#         mean         - vector of means 
#         Sigma        - vector of sigma2  
#         p            - vector of mixture weights 
#output: matrix (niter X k) containiing probabilities
#        for each of the data point belong to each of the components          
#        which are used as probabilities for the multinomial distribution
#        when sampling the latent variable Z
#############################################################

fulldatajointloglike = function(y,mean,Sigma,p,beta){
  output = matrix(NA,length(y),length(p))
  
    for(kp in 1:length(p)){

        if (length(Sigma)==length(p)) 
         {sigma2_t=Sigma[kp]
         }else{
          sigma2_t=Sigma
         }


      output[,kp]=log(p[kp])+dnorm(y,mean[kp],sqrt(sigma2_t),log=TRUE)
    }
  
  return(output-apply(output,1,max))
}


#############################################################
# sampleZ - samples Z as multinomial using probabilities obtained from 
# fulldatajointloglike
# input:        pij - probabilities obrtained from 
#                     mpij
# output:       sampled indicator variable Z
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
#               beta                 - temperature
# output:       sampled means 
#############################################################
sampleMeans = function(y,Z,Nk,sigma2,Meanpriorpars,beta){
  sumY  = Z%*%y
  denom = sigma2+beta*Nk*Meanpriorpars[,2]
  mu    = (beta*sumY*Meanpriorpars[,2]+sigma2*Meanpriorpars[,1])/denom
  sd    = sqrt(sigma2*Meanpriorpars[,2]/denom)
  rnorm(n=length(Nk),mu,sd)
}


#############################################################
# conditionalloglike:  Tempered Likelihood is Gaussian
# The tempered likelihood is evaluated, as well as
# input:  y           - Galaxy data from MASS package
#         means       - vector of means
#         Sigma       - vector of sigma2
#         p           - vector of mixture weights
#         Z           - matrix of indicator variable
#         beta        - value of beta
#         k           - number of components
# output: evaluated tempered likelihood
#############################################################
conditionalloglike = function(y,means,Sigma,p,Z,beta,k=3){
    output = matrix(NA,length(y),length(p))
      
    output = dnorm(y,mean=apply(t(Z)*matrix(means,length(y),k,byrow=T),1,sum),
                            sd=apply(t(Z)*matrix(sqrt(Sigma),length(y),k,byrow=T),1,sum),log=TRUE)
    
     return(beta*sum(output))
}

#############################################################
# marginalllik:  untempered Likelihood 
# The untempered likelihood is evaluated, which is later
# used to obtain the marginal likelihood
# input:  y           - Galaxy data from MASS package
#         means       - vector of means
#         sigma2       - vector of sigma2
#         Z           - matrix of indicator variable
#         k           - number of components
# output: evaluated untempered likelihood
#############################################################
marginalllik=function(y,means,sigma2,Z,k,p){
  
#  mllik     = dnorm(y,apply(t(Z)*matrix(means,length(y),k,byrow=T),1,sum),
#                    apply(t(Z)*matrix(sqrt(sigma2),length(y),k,byrow=T),1,sum),log=TRUE)
#  return(sum(mllik))

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

  return(mllik)

}

#################################################################################
# sampleparams1chain: samples parameters of interest
# input:  y              - Galaxy data from MASS package
#         means          - vector of means
#         sigma2         - vector of sigma2
#         pmatrix        - vector of mixture weights
#         k              - number of components
#         Meanpriorpars  - mean and st.dev for the prior of means
#         Sigmapriorpars - shape and scale for the prior of sigma2
#         Ppriorpars     - prior of p
#         beta           - value of beta
# output: sampled parameters of interest at fixed temperature
#         mean           - sampled means
#         pmat           - sampled mixture probabilities
#         sigma2         - sampled sigma2
#         SSE            - SSE used to obtain samples for sigma2 
#         Z=Z            - indicator variable
##################################################################################

sampleparams1chain = function(y,means,sigma2,pmatrix,k=3,Meanpriorpars,Sigmapriorpars,Ppriorpars,beta){
  
  
  # get the probabilities for the Zij values
  temp = fulldatajointloglike(y,means,sigma2 ,pmatrix,beta) 
  pij  = exp(temp)/matrix(apply(exp(temp),1,sum),length(y),k)  
  # Now sample a Z from a multinomial distribution.  And find the number of galaxies per group
  Z    = sampleZ(pij)
  
  # sample the p
  Nk         = apply(Z,1,sum)
  pmatrixnew = rdirichlet(1,Ppriorpars+Nk)
  
  
  # sample the means
  meansnew   = sampleMeans(y,Z,Nk,sigma2,Meanpriorpars,beta)
  
  
  # sample the variances
  SSE        = SSEfun(y,Z,meansnew,k)
  if (length(sigma2)==1){SSE=sum(SSE);Nk=sum(Nk);k=1}
  sigma2new  = 1/rgamma(k,shape=Nk*beta/2+Sigmapriorpars[1],scale=1/(SSE*beta/2+Sigmapriorpars[2]))    
  
  return(list(mean=meansnew,pmat=pmatrixnew,sigma2=sigma2new,SSE=SSE,Z=Z))
}

##################################################################################
#END!!!
##################################################################################