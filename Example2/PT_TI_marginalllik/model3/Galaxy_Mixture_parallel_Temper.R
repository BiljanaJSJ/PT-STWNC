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


####################################################################################################
#Calculate the thermodinamic integral via Parallel Tempering(PT)
#run PT with optimal temperature schedule as in Calderhead and Girolami (2009) (ti=(i/N)^5) with N=30 chains
#This code runs a model with 3 components with unequal variances
####################################################################################################



#parallel chains temperature parameters



runTI_PT=function(niter=50000,data=galaxies,k=k1,Meanpriorpars= cbind(rep(20,k1),rep(100,k1)),
                  Sigmapriorpars=c(3,1/20),
                  Ppriorpars= rep(1,k1),beta=beta,equalVar=NULL,eval_name){


npar    = length(beta)
M      = npar
y     = data/1000 ## normalize the data 
n     = length(y)


# Set up the matrices
# from all tempered chains:
means   = list()
sigma2  = list()
pmatrix = list()
Z       = list()
mllik   = matrix(NA,niter,npar)
for(kp in 1:npar){
  means[[kp]]       = matrix(NA,k,niter)
  means[[kp]][,1]   = rnorm(k, mean=mean(y), sd=2*sd(y))
  
  if (is.null(equalVar)){

  sigma2[[kp]]      = matrix(NA,k,niter)
  sigma2[[kp]][,1]  = rep(.5,k)

  }else{
  sigma2[[kp]]      = rep(NA,niter)
  sigma2[[kp]][1]   = .5

 
  }

  pmatrix[[kp]]     = matrix(NA,k,niter)
  pmatrix[[kp]][,1] = rep(1/k, k)

  if (is.null(equalVar)){ 
   sigma2_1=sigma2[[kp]][,1]
  }else{
  sigma2_1=sigma2[[kp]][1]
  }

 # temp = fulldatajointloglike(y,means[[kp]][,1],sigma2_1 ,pmatrix[[kp]][,1],beta=beta) 
#  pij  = exp(temp)/matrix(apply(exp(temp),1,sum),length(y),k)  
}

# fill a matrix to keep track of swap proposal / acceptance
swappers = matrix(0,niter,3)


for(iter in 2:niter){

  # potentially select two chains and propose a swap between parameters
  # we don't really want to swap every time, but almost every time would be nice.  
  # I control this using a lazy algorithm:
  maybeswap = sample(0:(npar+1),2,replace=TRUE)
  
if(all(length(unique(maybeswap))==2,min(maybeswap)>0,max(maybeswap)<=npar,iter>2)){
    #propose a parameter swap


    if (is.null(equalVar)){
    sigma2_st      = sigma2[[maybeswap[1]]][,iter-1]
    }else{
    sigma2_st      = sigma2[[maybeswap[1]]][iter-1]  
    }

    
    C1like1 = conditionalloglike(y,means[[maybeswap[1]]][,iter-1],sigma2_st,pmatrix[[maybeswap[1]]][,iter-1],Z[[maybeswap[1]]],beta[maybeswap[1]],k)
    C2like2 = conditionalloglike(y,means[[maybeswap[2]]][,iter-1],sigma2_st,pmatrix[[maybeswap[2]]][,iter-1],Z[[maybeswap[2]]],beta[maybeswap[2]],k)
    
    C2like1 = conditionalloglike(y,means[[maybeswap[1]]][,iter-1],sigma2_st,pmatrix[[maybeswap[1]]][,iter-1],Z[[maybeswap[1]]],beta[maybeswap[2]],k)
    C1like2 = conditionalloglike(y,means[[maybeswap[2]]][,iter-1],sigma2_st,pmatrix[[maybeswap[2]]][,iter-1],Z[[maybeswap[2]]],beta[maybeswap[1]],k)
    
    if(runif(1)< exp(C1like2+C2like1  -C1like1 - C2like2)){
      #accept the swap
      #exchange parameters information between the swapping chains
      means[[maybeswap[1]]][,iter]   = means[[maybeswap[2]]][,iter-1]
      means[[maybeswap[2]]][,iter]    = means[[maybeswap[1]]][,iter-1]

      pmatrix[[maybeswap[1]]][,iter] = pmatrix[[maybeswap[2]]][,iter-1]
    
      pmatrix[[maybeswap[2]]][,iter]  = pmatrix[[maybeswap[1]]][,iter-1]
      
      if (is.null(equalVar)){

      sigma2[[maybeswap[1]]][,iter]  = sigma2[[maybeswap[2]]][,iter-1]
      sigma2[[maybeswap[2]]][,iter]   = sigma2[[maybeswap[1]]][,iter-1]
       }else{

      sigma2[[maybeswap[1]]][iter]  = sigma2[[maybeswap[2]]][iter-1]
      sigma2[[maybeswap[2]]][iter]   = sigma2[[maybeswap[1]]][iter-1]
 
       }

      
      swappers[iter,]=c(sort(maybeswap),1)
	  
      tmp               = Z[[maybeswap[1]]]
      Z[[maybeswap[1]]] = Z[[maybeswap[2]]]
      Z[[maybeswap[2]]] = tmp
	 
      
    }else{
      # no swap accepted
      # copy the current values of the parameters of interest
      # from the last accepted value
      means[[maybeswap[1]]][,iter]   = means[[maybeswap[1]]][,iter-1]
      means[[maybeswap[2]]][,iter]   = means[[maybeswap[2]]][,iter-1]

      pmatrix[[maybeswap[1]]][,iter] = pmatrix[[maybeswap[1]]][,iter-1]
      pmatrix[[maybeswap[2]]][,iter] = pmatrix[[maybeswap[2]]][,iter-1]
      
    
      if (is.null(equalVar)){

      sigma2[[maybeswap[1]]][,iter]  = sigma2[[maybeswap[1]]][,iter-1]
      sigma2[[maybeswap[2]]][,iter]  = sigma2[[maybeswap[2]]][,iter-1]
      
      }else{

      sigma2[[maybeswap[1]]][iter]  = sigma2[[maybeswap[1]]][iter-1]
      sigma2[[maybeswap[2]]][iter]  = sigma2[[maybeswap[2]]][iter-1]

      }



      swappers[iter,]=c(sort(maybeswap),0)
    }
    
    #Mutation step
    #The other chains still need to move around
    index=1:npar
    index = index[-maybeswap]
    for(chain in index){


      if (is.null(equalVar)){
       sigma2_t1  = sigma2[[chain]][,iter-1]
      }else{
       sigma2_t1  = sigma2[[chain]][iter-1]
      }



      output = sampleparams1chain(y,means=means[[chain]][,iter-1],sigma2= sigma2_t1 ,pmatrix=pmatrix[[chain]][,iter-1],k=k,Meanpriorpars,Sigmapriorpars,Ppriorpars,beta=beta[chain])
   

        means[[chain]][,iter]   = output$mean
   
         if (is.null(equalVar)){
         sigma2[[chain]][,iter]  = output$sigma2
  
         }else{
         sigma2[[chain]][iter]  = output$sigma2
         }
  


      pmatrix[[chain]][,iter] = output$pmat
      Z[[chain]]              = output$Z
    }
    
    
  }else{
    #skip the parameter swap and update as usual
 
    for(chain in 1:npar){


        if (is.null(equalVar)){
           sigma2_s  =sigma2[[chain]][,iter-1]
        }else{
           sigma2_s  = sigma2[[chain]][iter-1]
        }

      output = sampleparams1chain(y,means=means[[chain]][,iter-1],sigma2= sigma2_s ,pmatrix[[chain]][,iter-1],k=k,Meanpriorpars,Sigmapriorpars,Ppriorpars,beta[chain])
    
      means[[chain]][,iter]   = output$mean


         if (is.null(equalVar)){
         sigma2[[chain]][,iter]  = output$sigma2  
         }else{
         sigma2[[chain]][iter]  = output$sigma2
         }
  

      
      pmatrix[[chain]][,iter] = output$pmat
      Z[[chain]]              = output$Z
    }
  }
 
   for(chain in 1:npar){

      if (is.null(equalVar)){
         sigma2_mllik  = sigma2[[chain]][,iter]
       }else{
         sigma2_mllik  = sigma2[[chain]][iter]
       }



      mllik[iter,chain]=marginalllik(y,means[[chain]][,iter],sigma2_mllik,Z[[chain]],k=k,p=pmatrix[[chain]][,iter])
   }

if ((iter %% 1000) == 0) {
out_ls=list(means=means,sigma2=sigma2,pmatrix=pmatrix,Z=Z,mllik=mllik)
save(out_ls,file=paste('PT_ThermodynInt_',k,'_',iter,'.RData',sep=''))
}
}

 save(out_ls,file=paste('mllik',eval_name,'.RData',sep=''))

#evaluate the posterior needed for calculations of the Kullberg Libler divergence (equation (41) in the paper)

sequence=1:M
beta    = (sequence/M)^5
npar    = length(beta)



pn_1=pn=matrix(NA,niter,npar)
pn_1s=pns=matrix(NA,niter,npar)
for (i in (2:npar)){
for (l in (1:niter)){  

        if (is.null(equalVar)){
         sigma2_kl1  =out_ls$sigma2[[i]][,l]
         }else{
         sigma2_kl1  =out_ls$sigma2[[i]][l]
         }
         pn_1s[l,i] = posterior(y,means=out_ls$means[[i]][,l],Sigma=sigma2_kl1,p=out_ls$pmatrix[[i]][,l],Z=out_ls$Z[[i]],beta=beta[i],k=k,Meanpriorpars=Meanpriorpars,priorpars=Sigmapriorpars,Ppriorpars=Ppriorpars)
}
}
#save the evaluations of the posterior in EvallogPost.RData.
#Evaluated posterior saved in EvallogPost1.RData. corresponds to one run, for a given model
eval=list(pn_1s=pn_1s)
save(eval,file=paste('EvallogPost',eval_name,'.RData',sep=''))
return(out_ls)
}