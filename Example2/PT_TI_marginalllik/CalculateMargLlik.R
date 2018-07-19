
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

#######################################################################################################
#Calculate the thermodynamic integral
#from the samples obtained from 20 independent PT runs on 30 chains
#######################################################################################################
#checkpackages(entropy)
library(entropy)



listfolders=list.files('/Users/bjonoska/Documents/BiasGalaxyNewPriors_7_4_2016/PT_mllik', full.names = TRUE)
n=listfolders[14:18]

#setup the temperature schedule, it is the same scedule used by PT with 30 chains
sequence = 1:30
beta     = (sequence/30)^5
npar     = length(beta)
#discard the first half of the iterations
#low      = round(35000-35000/30)
low      = 17500
up       = 35000

#calculate the approximation to the thermodynamic integral (equation (41) from the paper)
approx_int=matrix(NA,10,npar)
appThInt=matrix(NA,10,length(n))

for (o in (1:length(n))){
foldername=n[o]
approx_int=matrix(NA,10,npar)
for(m in (1:10)){

   #subdirName=paste0(foldername,'/run',m,sep='')
  


   out_ls=get(load(paste(foldername,'/mllik',m,'.RData',sep='')))
     
    

  for (i in (2:npar)){
    Ei_1=mean(out_ls$mllik[low:up,i-1][-c(1,which(out_ls$mllik[low:up,i-1]==-Inf))])
    E_i=mean(out_ls$mllik[low:up,i][-c(1,which(out_ls$mllik[low:up,i]==-Inf))])
    delta_t=beta[i]-beta[i-1]
    #approximation to the integral
    approx_int[m,i]=delta_t*(Ei_1+E_i)
    #KL-divergence
  }
  
  (appThInt[m,o]=(1/2)*sum(approx_int[m,-1]))
  
  
}
}
#save the calculated approximation for each of the 20 runs for each of the four models
#setwd('/stat/stat2/bjonoska/ST/newST/Gallaxy/PTThermodynamicInt')
save(appThInt,file='appThInt.RData')
#appThInt=get(load('appThInt.RData'))

#calculate the Kullbarg-Leibler divergence 
#Note that for each of the 10 runs, for each of the four models
#we already have calculate the needed evaluations of posterior distribution
#which are stored in the file EvallogPost.RData in each for each of the 4 models
#for each of the 10 runs. EvallogPost.RData is written by the Galaxy_Mixture_parallel_Temper.R script

KL1=KL2=matrix(NA,10,npar)
KL_TI=matrix(NA,10,length(n))

for (o in (1: length(n))){
  foldername=n[o]
  for (l in (1:10)){
       out_ls=get(load(paste(foldername,'/EvallogPost_upd',l,'.RData',sep='')))
      
       for (i in (2:npar)){
          KL1[l,i]=KL.plugin(out_ls$pn_1s[,i-1], out_ls$pn_1s[,i])
          KL2[l,i]=KL.plugin(out_ls$pn_1s[,i], out_ls$pn_1s[,i-1])
       }
      (KL_TI[l,o]=(1/2)*sum(KL1[l,-c(1,2)]-KL2[l,-c(1,2)]))
   }
  
}   
#save the file with calculated KL divergence
#The saved object is KL_TI which is a matrix with dimension (20 runs X 4 models)
save(KL_TI,file='KL_list.RData')

#appThInt=get(load('appThInt.RData'))
#KL_TI=get(load('KL_list.RData'))

#calculate the marginal likelihood
margLik=margLik_sd=matrix(NA,2,length(n))
for (i in (1:length(n))){
#TI according to Caldearhead and Girolami (2009)
(margLik[1,i]=mean(appThInt[,i]+KL_TI[,i]))
(margLik_sd[1,i]= sqrt(var(appThInt[,i]+KL_TI[,i]))    )
#TI according to Friel and Pettitt (2008)
(margLik[2,i]=mean(appThInt[,i]))
(margLik_sd[2,i]= sqrt(var(appThInt[,i]))    )

}


out_ls=list(margLik=margLik,margLik_sd=margLik_sd)
#setwd('/stat/stat2/bjonoska/ST/newST/Gallaxy/PTThermodynamicInt')
save(out_ls,file='marglikelihood.RData')




tbl=list(out_ls$margLik,out_ls$margLik_sd)
colnames(tbl[[1]])=c('model 2 comp. equal var','model 3 comp.','model 3 comp. equal var','model 4 comp.','model 5 comp.')
rownames(tbl[[1]])=c('marg. llik. CG','marg. llik. FP')
t(tbl[[1]])

colnames(tbl[[2]])=c('model 2 comp. equal var','model 3 comp.','model 3 comp. equal var','model 4 comp.','model 5 comp.')
rownames(tbl[[2]])=c('SE CG','SE FP')
t(tbl[[2]])



> library(xtable)
> xtable(t(tbl[[1]]))
% latex table generated in R 3.2.2 by xtable 1.8-2 package
% Mon Dec 12 10:08:57 2016
\begin{table}[ht]
\centering
\begin{tabular}{rrr}
  \hline
 & marg. llik. CG & marg. llik. FP \\
  \hline
model 2 comp. equal var & -238.02 & -238.03 \\
  model 3 comp. & -224.26 & -224.28 \\
  model 3 comp. equal var & -224.20 & -224.23 \\
  model 4 comp. & -223.88 & -223.89 \\
  model 5 comp. & -223.85 & -223.85 \\
   \hline
\end{tabular}
\end{table}


 t(tbl[[2]])
                             SE CG      SE FP
model 2 comp. equal var 0.02533089 0.02323573
model 3 comp.           0.03003759 0.03306666
model 3 comp. equal var 0.05704219 0.04962223
model 4 comp.           0.02497470 0.02185081
model 5 comp.           0.01609339 0.01604633

