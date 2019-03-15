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

#########################################################
#List of functions in this file

#1.runST            - runs the PT-STWNC algorithm
#2.initialization   - initializes the PT-STWNC algorithm
#3.STstep_tau       - samples the auxiliary parameter tau  
#########################################################


#########################################################
#User specified functions needed for the runST function

#1.loglik          - tempered log likelihood
#2.mloglik         - unitempered log likelihood
#3.STstep_pars     - sample parameters of interest
#4.OptimizePars    - optimize parameters of interest at 
#                    current tau  and proposed tau 
#5.posterior_notau - tempered posterior distribution
#                    without prior for tau
#########################################################




######################################################################################
#STstep_tau: Transition step two of the STWTDNC:
#            update tau with fixed parameters
#input:       

#          pars      - a vector of sampled parameters
#          tau       - a scalar of sampled tau
#          y         - data vector
#          PriorPars - vector of prior parameters given 
#                      in a sequence of appearence of the sampled parameters in the vector pars
#          acc       - a scalar with counts of the accepted tau 
#          tune      - a scalar of tuning for the transition variance of tau
#          parAdd    - a list of additional parameters, can be either sampled 
#                      parameteres or additional parameters
   
#output:   a list of   
#         theta      - sampled tau in the vector c(pars, tau)
#         log_r_bot  - denominator of the acceptance ratio
#         accepts1   - counts of accepted tau
#         log_r_top1 - numerator of the acceptance ratio
#         tau_prop        - proposed tau
######################################################################################
STstep_tau = function(pars,tau,y,PriorPars,acc=accepts1,tune=NULL,parAdd=NULL, cl=NULL){
	
	if (is.na(tune)){
		tau_prop       = runif(n=1,0,1) 
	}else{
		tau_prop       = rtruncnorm(1,a=0,b=1,tau,tune)
	}
	
	
	Optim=OptimizePars(y=y,tau_prop=tau_prop,pars=c(pars,tau),PriorPars=PriorPars,parAdd = parAdd, cl=cl)
	
	max_prop=Optim$max_prop
	max_it=Optim$max_it
	
	if ((!is.null(parAdd)) & (!is.null(Optim$parAdd_max_prop))) {
	
	#calculate the MH acceptance prob in case 
	#when there are additional parameters in parAdd
	#that also need to be sampled	
	parAdd_max_prop = Optim$parAdd_max_prop
	parAdd_max_it   = Optim$parAdd_max_it
	log_r_top       = posterior(y=y,pars=pars,tau=tau_prop,PriorPars=PriorPars,par_max=max_prop,parAdd=parAdd,parAdd_max=parAdd_max_prop)
	log_r_bot       = posterior(y=y,pars=pars,tau=tau,PriorPars=PriorPars,par_max=max_it,parAdd=parAdd,parAdd_max=parAdd_max_it)
	}else if (((!is.null(parAdd)) & (is.null(Optim$parAdd_max_prop)))){
		
			#calculate the MH acceptance prob in case 
		#when there are additional parameters in parAdd
		#but do not need to be sampled	
		log_r_top       = posterior(y=y,pars=pars,tau=tau_prop,PriorPars=PriorPars,par_max=max_prop,parAdd=parAdd)
		log_r_bot       = posterior(y=y,pars=pars,tau=tau,PriorPars=PriorPars,par_max=max_it,parAdd=parAdd)
		
	}else{
	
		#calculate the MH acceptance prob in case 
		#when there are no additional parameters in parAdd
		log_r_top       = posterior(y=y,pars=pars,tau=tau_prop,PriorPars=PriorPars,par_max=max_prop)
		log_r_bot       = posterior(y=y,pars=pars,tau=tau,PriorPars=PriorPars,par_max=max_it)
		
	}
	
	pars                          = c(pars,tau)

	#calculate alpha
        if (is.na(tune)){
	alpha                         = log_r_top - log_r_bot
        }else{
	alpha                         = log_r_top - log_r_bot+log(dtruncnorm(tau,a=0,b=1,tau_prop,tune))-log(dtruncnorm(tau_prop,a=0,b=1,tau,tune))
        }

        # make a decision
	
	if (all(!is.na(alpha) , runif(1) < exp(alpha))){
		# accept the move
		acc                   = acc+1;
		pars[length(pars)]    = tau_prop
		#pars                  = c(pars,tau_prop)
		
	}
	return(list(theta=pars,	log_r_bot = log_r_bot,accepts1=acc,
			log_r_top1=log_r_top,tau_prop=tau_prop)) 
}
######################################################################################
#initialization: Initialize the algorithm 
#Input:
#                niter   - number of iterations for algorithm to run
#                nstops  - needed for calculation of the acceptance rate when tuning
#                kp1     - kp1 starts the counts of accepted values divided into nstop
#                          number of stops to tau
#                q1      - a vector of length two, first element to initialize the 'tune' parameter 
#                          that tunes the transition variance of the first parameter and second element is TRUE/FALSE indicating whether to tune the variance 
#                q2      - a vector of length two, first element to initialize the 'tune' parameter 
#                          that tunes the transition variance of the first parameter and second element is TRUE/FALSE indicating whether to tune the variance 
#                IniPar  - vector of initial values for the parameters 
#                          with tau added as a last element of the vector
#                nchains - number of chains
#                ttau    - a vector of length two, first element initalizes tune variance 
#                          for tau and second element is TURE/FALSE whihc indicates whether to tune variance of tau 
#                parAdd  - a list of additional parameters, can be either sampled 
#                          parameteres or additional parameters
#Output:      a list
#             log_r_bot        - a chain of posterior evaluations for updating the parameters
#             accepts          - counts of accepted updates of the parameters
#             acceptsTau       - counts of accepted updates of tau
#             acc_rate         - acceptance rate of the parameters
#             acc_rate1        - acceptance rate of the tau
#             tau_prop         - chain of proposals for tau
#             mllik            - maximum likelihood estimates for tempered and target chains
#             PT_chain         - chains of updated parameters and tau for tempered and target chains
#             swappers         - a matrix of swapped chains and accepted swaps in the third column
#             npar             - number of parameters of interest
#             kp               - starts the counts of accepted values divided into nstop
#                                number of stops for parameters of interest
#             tune_q1          - tunning chain for the variance of the transition kernerl of the first parameter
#             tune_q2          - tunning chain for the variance of the transition kernerl of the second parameter
#             parAdd_chain     - initialized chain for parAdd
#             tunetau          - tunning chain for the variance of the transition kernerl of tau
######################################################################################

initialization=function(niter,nstops,kp1,tune_pars_init,IniPar,nchains,ttau,parAdd){
	
    npar              = length(IniPar)-1
		#intialize the algorithm for the first run
    log_r_bot         = list()
 	  log_r_bot[[1]]    = rep(0,niter)
 	  log_r_bot[[2]]    = rep(0,niter)
	
		# log_r_bot1       = rep(0,niter)
		# log_r_top1       = rep(0,niter)
		
		accepts          = list()
		accepts[[1]]     = matrix(1,nstops,npar)
		accepts[[2]]     = matrix(1,nstops,npar)
		
		
		acc_rate          = list()
		acc_rate[[1]]     = matrix(1,nstops,npar)
		acc_rate[[2]]     = matrix(1,nstops,npar)
		
         kp                  = list()
         kp[[1]]             = rep(kp1,npar)
         kp[[2]]             = rep(kp1,npar)
    
		#kp                = rep(kp1,npar)
    

		acceptsTau        = rep(1,nstops)
		acc_rate1         = rep(1,nstops)
		
	#	tau_prop          = rep(0,niter) 
		
		# tune_q2            = list()
		# tune_q2[[1]]       = rep(NA,niter)
		# tune_q2[[2]]       = rep(NA,niter)
		# tune_q2[[1]][1]    = q2[1]
		# tune_q2[[2]][1]    = q2[1]
		# 
		# tune_q1            = list()
		# tune_q1[[1]]       = rep(NA,niter)
		# tune_q1[[2]]       = rep(NA,niter)
		# tune_q1[[1]][1]    = q1[1]
		# tune_q1[[2]][1]    = q1[1]
		
		# tune_pars            = list()
		# tune_pars[[1]]       = matrix(NA,niter,ncol=npar)
		# tune_pars[[2]]       = matrix(NA,niter,ncol=npar)
		# 
		# tune_pars[[1]][1,]   = tune_pars_init
		# tune_pars[[2]][1,]   = tune_pars_init
		
		
		tune_pars          = list()
		tune_pars[[1]]     = matrix(NA,nstops,npar)
		tune_pars[[2]]     = matrix(NA,nstops,npar)
		
		tune_pars[[1]][1,]   = tune_pars_init
		tune_pars[[2]][1,]   = tune_pars_init
		
		
		
		
		PT_chain           = list()
		PT_chain[[1]]      = matrix(NA,niter,ncol=length(IniPar))
		PT_chain[[2]]      = matrix(NA,niter,ncol=length(IniPar))
		PT_chain[[1]][1,]  = IniPar
		PT_chain[[2]][1,]  = IniPar
		PT_chain[[2]][,npar+1]  = 1
		

		swappers           = matrix(0,niter,nchains+1)
		mllik              = list()
		mllik[[1]]         = rep(NA,niter)
		mllik[[2]]         = rep(NA,niter)
		
		parAdd_chain=list()
		parAdd_chain[[1]]  = parAdd_chain[[2]]=list()
		parAdd_chain[[1]]  = parAdd
		parAdd_chain[[2]]  = parAdd
		
		
		# tunetau            = rep(NA,niter)
		# tunetau[1]         = ttau[1]
		
		tunetau            = rep(NA,nstops)
		tunetau[1]         = ttau[1]

	return(list(log_r_bot        = log_r_bot,
							# log_r_bot1       = log_r_bot1,
							# log_r_top1       = log_r_top1,
							accepts          = accepts,
							acceptsTau       = acceptsTau,
							acc_rate         = acc_rate,
							acc_rate1        = acc_rate1,
					#		tau_prop         = tau_prop,
							mllik            = mllik,
							PT_chain         = PT_chain,
							swappers         = swappers,
							npar             = npar,
							kp               = kp,
							#tune_q1          = tune_q1,
							#tune_q2          = tune_q2,
							tune_pars        = tune_pars,
							parAdd_chain     = parAdd_chain,
							tunetau          = tunetau))
}
######################################################################################



# runST runs ST algorithm.
# input:    niter       - number of iterations for algorith to run
#           y           - data
#           PriorPars   - priors for all of the parameters  
#           IniPar      - vector of initial values for the parameters 
#                         with tau added as a last element of the vector
#           q1          - a vector of length two, first element to initialize the 'tune' parameter 
#                         that tunes the transition variance of the first parameter and second element is TRUE/FALSE indicating whether to tune the variance 
#           q2          - a vector of length two, first element to initialize the 'tune' parameter 
#                         that tunes the transition variance of the first parameter and second element is TRUE/FALSE indicating whether to tune the variance 
#           nstops      - needed for calculation of the acceptance rate when tuning
#           kp1         - kp1 starts the counts of accepted values divided into nstop
#                         number of stops to tau
#           nchains     - number of chains
#           ttau        - a vector of length two, first element initalizes tune variance 
#                         for tau and second element is TURE/FALSE whihc indicates whether to tune variance of tau 
#           parAdd      - a list of additional parameters, can be either sampled 
#                         parameteres or additional parameters

# output:      a list
#             PT_chain         - chains of updated parameters and tau for tempered and target chains
#             rate             - acceptance rate of the parameters
#             rate1            - acceptance rate of tau
#             accepts          - counts of accepted updates of the parameters
#             acceptsTau       - counts of accepted updates of tau
#             kp               - starts the counts of accepted values divided into nstop
#                                number of stops for parameters of interest
#             y                - data
#             tau_prop         - chain of proposals for tau
#             log_r_top        - a chain of posterior evaluations for updating tau
#             mllik            - maximum likelihood estimates for tempered and target chains
#             swappers         - a matrix of swapped chains and accepted swaps in the third column
#             tunetau          - tunning chain for the variance of the transition kernerl of tau
#             tune_q1          - tunning chain for the variance of the transition kernerl of the first parameter
#             tune_q2          - tunning chain for the variance of the transition kernerl of the second parameter
######################################################################################

runST=function(niter = 500,y,PriorPars,IniPar,tune_pars_init,ttune_pars,nstops=20,kp1=1,kptau=1, nchains=2,ttau=c(NA, FALSE),parAdd=NULL, cl=NULL, whichTune=NULL){


	#start the PT with two chains
	
	#initializee the algorithm
	init=initialization(niter=niter,nstops=nstops,kp1=kp1,tune_pars_init=tune_pars_init,IniPar=IniPar,parAdd=parAdd,nchains=nchains,ttau=ttau)
	
	     
	
	log_r_bot          = init$log_r_bot
	# log_r_bot1         = init$log_r_bot1
	# log_r_top1         = init$log_r_top1

	accepts            = init$accepts
	acceptsTau         = init$acceptsTau
	acc_rate           = init$acc_rate
	acc_rate1          = init$acc_rate1
#	tau_prop           = init$tau_prop  
	mllik              = init$mllik
	PT_chain           = init$PT_chain
	swappers           = init$swappers
	npar               = init$npar
	kp                 = init$kp
#	tune_q1            = init$tune_q1
#	tune_q2            = init$tune_q2
	tune_pars          = init$tune_pars
	parAdd_chain       = init$parAdd
	tunetau            = init$tunetau 

	
	#set up the starting index 
	index            = which(is.na(PT_chain[[1]][,2]))[1]
	

	#initialize log_r_bot for chain 1
	log_r_bot[[1]][1]  =  	posterior_notau(y=y,pars=PT_chain[[1]][1,(1:npar)],tau=PT_chain[[1]][1,npar+1],log=T,PriorPars=PriorPars,parAdd=parAdd)$output

	#initialize log_r_bot for chain 2
	log_r_bot[[2]]     =  log_r_bot[[1]]
	
	out = list()

	
	for(iter in index:niter){
		print(iter)
		#propose chains for swap
		maybeswap = sample(0:(nchains+1),2,replace=TRUE)
		
		if(all(length(unique(maybeswap))==2,min(maybeswap)>0,max(maybeswap)<=nchains,iter>2)){
			#calcluate tempered likelihood
			l11=loglik(x=y,pars=PT_chain[[maybeswap[1]]][iter-1,(1:npar)],tau=PT_chain[[maybeswap[1]]][iter-1,npar+1],log=T,parAdd=parAdd_chain[[maybeswap[1]]])$out
			l22=loglik(x=y,pars=PT_chain[[maybeswap[2]]][iter-1,(1:npar)],tau=PT_chain[[maybeswap[2]]][iter-1,npar+1],log=T,parAdd=parAdd_chain[[maybeswap[2]]])$out
			
			l21=loglik(x=y,pars=PT_chain[[maybeswap[1]]][iter-1,(1:npar)],tau=PT_chain[[maybeswap[2]]][iter-1,npar+1],log=T,parAdd=parAdd_chain[[maybeswap[1]]])$out
			l12=loglik(x=y,pars=PT_chain[[maybeswap[2]]][iter-1,(1:npar)],tau=PT_chain[[maybeswap[1]]][iter-1,npar+1],log=T,parAdd=parAdd_chain[[maybeswap[2]]])$out
			
			
			#accept or reject the proposed swap
			if(runif(1)< exp(l12 + l21 -l11 - l22)){
				#accept the swap
				#exchange the parameter values of the two chains 
				PT_chain[[maybeswap[1]]][iter,(1:npar)] =PT_chain[[maybeswap[2]]][iter-1,(1:npar)]
				PT_chain[[maybeswap[2]]][iter,(1:npar)] =PT_chain[[maybeswap[1]]][iter-1,(1:npar)]
				
				
				if (!is.null(parAdd)){
					tmp                = parAdd_chain[[1]]
					parAdd_chain[[1]]  = parAdd_chain[[2]]
					parAdd_chain[[2]]  = tmp
				}
				
				
				swappers[iter,]=c(maybeswap,1)
				
			}else{
				# no swap accepted
				#copy the parameter values from the last accepted values
				PT_chain[[maybeswap[1]]][iter,] =PT_chain[[maybeswap[1]]][iter-1,]
				PT_chain[[maybeswap[2]]][iter,] =PT_chain[[maybeswap[2]]][iter-1,]
		
				swappers[iter,]=c(maybeswap,0)
			}
			
			#setup the current tau1 to copy the last accepted value of tau from the 'tempered' chain
			PT_chain[[1]][iter,npar+1] = PT_chain[[1]][iter-1,npar+1]
			#setup the current tau2 =1
			#PT_chain[[2]][iter,npar+1] = 1
			
		}else{
			
			#mutation step:
			#update the 'tempered' chain via STWTDNC     
			chain=1
			
			
			#update the parameters of interest at fixed value of tau1 (current tau1)
			out          =  STstep_pars(pars=PT_chain[[chain]][iter-1,(1:npar)],
										            	tau=PT_chain[[chain]][iter-1,npar+1],
											            PriorPars=PriorPars,
											            y=y,log_r_bot=log_r_bot[[chain]][iter-1],
											            acc= sapply(1:npar,function(x){accepts[[chain]][kp[[chain]][x],x]}),
											            tune_pars=sapply(1:npar,function(x){tune_pars[[chain]][kp[[chain]][x],x]}),
											            ttune_pars = ttune_pars,
											            parAdd=parAdd_chain[[chain]])
			
			
			
			(PT_chain[[chain]][iter,]                    = out$theta)
			for (i in (1:npar)){
			accepts[[chain]][kp[[chain]][i], i]          = out$accepts[i]
			}
			log_r_bot[[chain]][iter]                     = out$log_r_bot
			
			if (!is.null(parAdd)){
				parAdd_chain[[chain]]                    = out$parAdd
			}
			#update tau1 while parameters of interest are kept at their last accepted values
			
	
	out1                = tryCatch({
			               STstep_tau(pars=PT_chain[[chain]][iter,(1:npar)],
	 								              tau=PT_chain[[chain]][iter,npar+1],
									              y=y,
									              PriorPars=PriorPars,
									              acc=acceptsTau[kptau],
									              parAdd=parAdd_chain[[chain]],
									              tune=tunetau[kptau],
									              cl=cl)
				
			},
			error=function(e){
				message(paste("er:",e,sep=''))
				return(NA)
			},                                
			warning=function(wr) {
				message(paste("wr:",wr,sep=''))
				return(NULL)
			}
			)   
			if (is.list(out1)) {
				
			(PT_chain[[chain]][iter,]  = out1$theta)
			
			acceptsTau[kptau]          = out1$accepts1
	#		log_r_bot1[iter]           = out1$log_r_bot
	#		tau_prop[iter]             = out1$tau_prop
	#		log_r_top1[iter]           = out1$log_r_top1
			
			
			}else{
				

				(PT_chain[[chain]][iter,]  = PT_chain[[chain]][iter-1,])
				
			#	log_r_bot1[iter]           = log_r_bot1[iter-1]
			#	tau_prop[iter]             = tau_prop[iter-1]
			#	log_r_top1[iter]           = log_r_top1[iter-1]
				
			}
			
			
			#mutation step for the 'target' chain
			
			chain=2
			#update the parameters of ineterst at tau2=1
			out2 =  STstep_pars(pars=PT_chain[[chain]][iter-1,(1:npar)],
					    						tau=1,
							    				PriorPars=PriorPars,
									   		  y=y,log_r_bot=log_r_bot[[chain]][iter-1],
											    acc=sapply(1:npar,function(x){accepts[[chain]][kp[[chain]][x],x]}),
											    tune_pars=sapply(1:npar,function(x){tune_pars[[chain]][kp[[chain]][x],x]}),
											    ttune_pars = ttune_pars,
											    parAdd=parAdd_chain[[chain]])
			
			
			(PT_chain[[chain]][iter,]                    = out2$theta)
			#mu_last                                      = PT_chain[[chain]][iter-1,1]
		#	accepts[[chain]][kp[[chain]][[chain]],]               = out$accepts
			for (i in (1:npar)){
				accepts[[chain]][kp[[chain]][i], i]                = out2$accepts[i]
			}
			#accepts[[chain]][kp[[chain]][[chain]][2],2]           = out$accepts[2]
			log_r_bot[[chain]][iter]                     = out2$log_r_bot
					
	  	    if (!is.null(parAdd)){
			 	   parAdd_chain[[chain]]                 = out2$parAdd
				
	   	    }
		}	

		#calculate the marginal likelihood for the two chains
		for (chain in (1:nchains)){
		mllik[[chain]][iter]    = loglik(x=y,pars=PT_chain[[chain]][iter,(1:npar)],tau=PT_chain[[chain]][iter,npar+1],log=T,parAdd=parAdd_chain[[chain]], mllik_eval=T)$mllik
		}
		# chain=2
		# mllik[[chain]][iter]    = loglik(x=y,pars=PT_chain[[chain]][iter,(1:npar)],tau=PT_chain[[chain]][iter,npar+1],log=T,parAdd=parAdd_chain[[chain]])$mllik
		# 
		
		
		for (chain in (1:nchains)){
			
		
			for (i in (1:npar)){
				
				
			if ((iter== kp[[chain]][i]*niter/nstops) & (iter!=niter)) {  
			 	  
				
			   acc_rate[[chain]][kp[[chain]][i],i]    = (accepts[[chain]][kp[[chain]][i],i]+1)/(2+(niter/nstops));				
			   kp[[chain]][i]                         =  kp[[chain]][i]+1;
			   tune_pars[[chain]][kp[[chain]][i],i]   = tune_pars[[chain]][kp[[chain]][i]-1,i]
			   
			   if ((ttune_pars[i]) & (iter< niter/2)) {
			   	
			   #if(iter< niter/2){
					
					if(acc_rate[[chain]][kp[[chain]][i]-1,i]>.29 || acc_rate[[chain]][kp[[chain]][i]-1,i]<.19 ){
						# adjust the transition density if we are in the first half of the iterations. 
                        				
						#tune transition step  tune_pars
						#tune_pars[[chain]][iter,i] = tune_pars[[chain]][iter,i]*acc_rate[[chain]][kp[[chain]][[chain]][i]-1,i]/.24
						tune_pars[[chain]][kp[[chain]][i],i] = tune_pars[[chain]][kp[[chain]][i],i]*acc_rate[[chain]][kp[[chain]][i]-1,i]/.24
						}
						
					}
					
			  # }
			}
			}
			
			
		}
				
				

		
		
		if ((iter == kptau*niter/nstops) & (iter!=niter)) {  
			
			
			acc_rate1[kptau]                       = (acceptsTau[kptau]+1)/(2+(niter/nstops));
			kptau                                  =  kptau+1;
			tunetau[kptau]                         = tunetau[kptau-1]
			# adjust the transition density if we are in the first half of the iterations. 
			if (ttau[2]){
			if(iter< niter/2){
				if(acc_rate1[kptau-1]>.28 || acc_rate1[kptau-1]<.18 ){
					tunetau[kptau]          = tunetau[kptau]*acc_rate1[kptau-1]/.23
				}
			}
		}
		}
		#     write the output at every 1000-th iterations
		#     save all the sampled parameters and make some plots
		if ((iter %% 1000) == 0) {
			out_ls=list(PT_chain=PT_chain,rate=acc_rate,rate1=acc_rate1,accepts=accepts,acceptsTau=acceptsTau,log_r_bot=log_r_bot,
									kp=kp,y=y,mllik=mllik,swappers=swappers,tunetau=tunetau,tune_pars=tune_pars,parAdd_chain=parAdd_chain)
			
			save(out_ls, file=paste("ST",iter,".RData",sep=""))
		 }
		 
	}
	
	return(list(PT_chain=PT_chain,rate=acc_rate,rate1=acc_rate1,accepts=accepts,acceptsTau=acceptsTau,log_r_bot=log_r_bot,
							kp=kp,y=y,mllik=mllik,swappers=swappers,tunetau=tunetau,tune_pars=tune_pars,parAdd_chain=parAdd_chain))
}
######################################################################################
#END!!