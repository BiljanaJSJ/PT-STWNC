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




##########################################################################
#run the PAWN algorithm
#The code for the algorithm is by Dr. Luke Bornn
##########################################################################

############################################
#set up the prior for mu
############################################
dprior_mu=function(x,mean_mu=0,k=1,log=T){
  return(sum(dnorm(x,mean=mean_mu,sd=k,log=log)))
}
############################################
#set up the log likelihood function
############################################

loglik=function(x,mean_mu,sigma=1,tau,log=T){
  
  if (log==T)
  {
    out=(1/tau)*sum(dnorm(x,mean=abs(mean_mu),sd=sqrt(sigma),log=log))
  }else{
    out=(prod(dnorm(x,mean=abs(mean_mu),sd=(sigma),log=log))^(1/tau))
  }
  mllik=sum(dnorm(x,mean=abs(mean_mu),sd=sqrt(sigma),log=log)) 
  return(list(out=out,mllik=mllik))
}

############################################
#set up the posterior for mu
############################################

posterior_mu=function(n=25,y,tau,mu,k=1,sigma=1,log=T){
  
  llik  = loglik(y,mean_mu=mu,sigma=sigma,tau=tau,log=log)$out
  pr_mu = dprior_mu(x=mu,k=k,log=log)
  
  
  if (log==T){ 
    out_mu=llik+pr_mu
  }else{
    out_mu=llik*pr_mu
  }
  
  return(out_mu)
}
############################################
#get the data
############################################
y=get(load('ST35000_power.RData'))$y


#make plot of posterior of mu
mugrid=seq(-3,3,len=1000)
post_mu=rep(0,1000)
for (i in (1:1000)){
  post_mu[i]=posterior_mu(n=25,y=y,tau=1,mu=mugrid[i],k=1,sigma=1,log=T)
}
plot(mugrid,post_mu, type="l")

#library(parallel)
#library(foreach)
#registerDoMC()

############################################################################
#The PAWN algorithm main function
############################################################################
big_loop <- function(N, T_max, L, Stoch_Approx, wang_landau, SA_seq, c, adapt, M=10){
  
  
  
  WL_seq <- SA_seq
  k <- 1
  last_step_change <- 0
  
  log_theta <- matrix(NA, nrow=N, ncol=T_max)
  log_theta[1,] <- rep(1,T_max)
  target_theta_dist <- rep(1/T_max,T_max)
  
  # Sampler
  X_sample <- matrix(NA,nrow=N, ncol=M)
  T_sample <- matrix(NA,nrow=N, ncol=M)
  X_sample[1,] <- rep(-10,M)
  if(M==1){
    T_sample[1,] <- 1
  }else{
    T_sample[1,] <- c(1,sample(1:T_max,1))
  }
  
  acceptX <- 0
  acceptT <- 0
  
  step_size = rep(NA,N)
  step_size[1] = 10
  
  for(i in 2:N){
    # Update X
    X_proposed <- X_sample[i-1,] + rnorm(M,0,sd=step_size[i-1])
    
    for(j in 1:M){
      
      top=posterior_mu(n=25,y=y,tau=T_sample[i-1,j],mu=X_proposed[j],k=1,sigma=1,log=T)
      bottom=posterior_mu(n=25,y=y,tau=T_sample[i-1,j],mu=X_sample[i-1,j],k=1,sigma=1,log=T)
      
      if(log(runif(1))< top-bottom){
        #if(log(runif(1))< target(X_proposed[j],T_sample[i-1]) - target(X_sample[i-1,j],T_sample[i-1])){
        X_sample[i,j] <- X_proposed[j]
        acceptX <- acceptX+1
      }else{
        X_sample[i,j] <- X_sample[i-1,j]
      }}
    
    if(adapt){
      step_size[i] =  max(10^(-10), step_size[i-1] + SA_seq[i] * (2 * ((acceptX)/(i*M) > 0.234) - 1))
      #cat(paste((acceptX)/(i*M),"   "))
    }else{
      step_size[i] = step_size[i-1]
    }
    
    # Update T
    for(j in 1:M){
      if(T_sample[i-1,j]==1){
        T_proposed <- 2
        proposal_ratio <- 2
      }else{
        if((T_sample[i-1,j]==T_max)){
          T_proposed <- T_max - 1
          proposal_ratio <- 2
        }else{
          T_proposed <- T_sample[i-1,j] + sample(c(-1,1),1)
          proposal_ratio <- 1
        }
      }
      
      top_tau=posterior_mu(n=25,y=y,tau=T_proposed,mu=X_sample[i,j],k=1,sigma=1,log=T)
      bottom_tau=posterior_mu(n=25,y=y,tau=T_sample[i-1,j],mu=X_sample[i,j],k=1,sigma=1,log=T)
      
      ratio <- exp(top_tau)/exp(bottom_tau)*proposal_ratio
      #ratio <- target(X_sample[i,j],T_proposed) - target(X_sample[i,j],T_sample[i-1]) - log(proposal_ratio)
      if(Stoch_Approx){ratio <- ratio + log_theta[i-1,T_sample[i-1,j]] - log_theta[i-1,T_proposed]}
      
      if (runif(1)< ratio){
        T_sample[i,j] <- T_proposed
        acceptT <- acceptT+1
      }else{
        T_sample[i,j] <- T_sample[i-1,j]
      }
    }
    #this is updating bias and checking the flat histogram rule  
    if(Stoch_Approx){
      indicator <- tabulate(T_sample[i,], nbins=10)/M
      if(wang_landau == TRUE){
        log_theta[i,] <- log_theta[i-1,] + WL_seq[k] * (indicator - target_theta_dist)
        if(i%%100==0){
          #cat(paste(max(abs(table(T_sample[last_step_change:i])/(i-last_step_change) - target_theta_dist)))," ")
          #on every 100 iteration check the flat histogram rule for the last 100 iterations
          if(max(abs(tabulate(T_sample[last_step_change:i], nbins=T_max)/(i-last_step_change) - target_theta_dist))< (c/T_max)){
            k <- k+1
            last_step_change <- i
          }
        }
      }else{
        log_theta[i,] <- log_theta[i-1,] + SA_seq[i] * (indicator - target_theta_dist)
      }
    }
    
    
  }
  return(list(acceptT=acceptT,acceptX=acceptX,X_sample=X_sample,T_sample=T_sample))
} 


sequence_builder <- function(N, t0)    t0/pmax(rep(t0,N),1:N)



do_it_step <- function(i){
  results <- rep(NA,7)
  cat(paste("Iteration", i, "of", length(N_vec), "\n"))
  results_WL <- big_loop(N=N_vec[i], T_max=T_max, L=L, Stoch_Approx=TRUE, wang_landau=TRUE, sequence_builder(N_vec[i],1),.01, FALSE)

  results[5] <-   sqrt(sum(results_WL$mean_vec^2)/L)

  return(results)
}


T_max <- 10
(N_vec = 273.45*2^(0:7))

#run the PAWN algorithm
do_it_adapt <- function(i){
  results <- rep(NA,6)
  cat(paste("Iteration", i, "of", length(N_vec), "\n"))
  results_WL_par <- big_loop(N=N_vec[i], T_max=T_max, L=L, Stoch_Approx=TRUE, wang_landau=TRUE, sequence_builder(N_vec[i],1),.1, TRUE, M=10)
 
  return(results_WL_par)
}

#collect the results and write in the file
result_adapt_500 <- do_it_adapt(8)
save(result_adapt_500,file='result_adapt_500.RData')

#generate plots
M=10
targetPost=list()
X_sample_smpled=result_adapt_500$X_sample[1000:N,]
T_samples_sampled=result_adapt_500$T_sample[1000:N,]
for (i in (1:M)){
  targetPost[[i]]=X_sample_smpled[which(T_samples_sampled[,1]==1),i]
  
}
tPost=do.call(rbind,targetPost)
png('margmutau.png')
par(mfrow=c(1,2))
hist(tPost,nclass=100,main='Marginal of mu',xlab=expression(mu),prob=TRUE)
hist(1/T_samples_sampled[,1],nclass=35,main='Marginal of tau',xlab=expression(tau),prob=TRUE)
dev.off()

png('tau_2_10_.png')
par(mfrow=c(3,3))
for (i in (2:10)){
  hist(1/T_samples_sampled[,i],nclass=35,main='Marginal of tau',xlab=expression(tau),prob=TRUE) 
}
dev.off()
