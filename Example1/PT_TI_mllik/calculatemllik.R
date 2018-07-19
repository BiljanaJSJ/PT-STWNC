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


y=get(load('ST35000_power.RData'))$y
n=length(y)

#obtain analitical  solution
a=n+1
b=sum(y)
b1=-b
c=sum(y^2)
I=log((1/(2*sqrt(a)*sqrt(2*pi)^n))*( exp(-c/2 + b^2/(2*a)  ) +  exp(-c/2 + b1^2/(2*a)  ) ))

#I=log(((2*pi)^(-n/2))*sqrt((n+1)^(-1))*exp( (1/2)*((sum(y)^2)/(n+1)) - (1/2)*(sum(y^2))    ))

save(I,file='analiticalSol.RData')


#run termodinamic integration via parallel tempering



source("Bimodal_parallel_Temper.R")
library(MASS)
library(entropy)

M=30
sequence = 1:M
#beta     = 0.05*(sequence/M)+0.95*(sequence/M)^3
beta     = (sequence/M)^5
npar     = length(beta)



PT_mllik=rep(NA,20)
mean_se=rep(NA,20)
n=25
low=1000
up=35000
npar=M


approx_int=rep(NA,npar)
appThInt=rep(NA,20)
KL1=KL2=matrix(NA,20,npar)
KL_TI=rep(NA,20)

Em_1=E_m=delta_t=rep(NA,npar)

for (i in (1:20)){
 out_ls=runPT(niter = 35000,M=M)
 save(out_ls,file=paste('TI_mllik',i,'.RData',sep=''))
 #out_ls=get(load(paste('TI_mllik',i,'.RData',sep='')))

   for (m in (2:npar)){
     Em_1[m]=mean(out_ls$mllik[low:up,m-1])
     E_m[m]=mean(out_ls$mllik[low:up,m])
     delta_t[m]=beta[m]-beta[m-1]
     #approximation to the integral
     approx_int[m]=delta_t[m]*(Em_1[m]+E_m[m])
     #KL-divergence
   }
  
   (appThInt[i]=(1/2)*sum(approx_int[-1]))

   for (m in (2:npar)){

     KL1[i,m]=KL.plugin(out_ls$pn_1s[low:up,m], out_ls$pn_1s[low:up,m-1])
     KL2[i,m]=KL.plugin(out_ls$pn_1s[low:up,m-1], out_ls$pn_1s[low:up,m])

   }
   (KL_TI[i]=(1/2)*sum(KL2[i,-1]-KL1[i,-1]))
     
#save the file with calculated KL divergence
#The saved object is KL_TI which is a matrix with dimension (20 runs X 4 models)



}
save(KL_TI,file='KL_list.RData')

save(appThInt,file='appThInt.RData')


margLik=margLik_sd=rep(NA,2)
#TI according to Calderhead and Girolami (2009)
(margLik[1]=mean(appThInt+KL_TI))
(margLik_sd[1]= sqrt(var(appThInt+KL_TI))  )   


#TI according to Friell and Pettit (2008)
(margLik[2]=mean(appThInt))
(margLik_sd[2]= sqrt(var(appThInt))  )   



mllik_out_ls=list(margLik=margLik,margLik_sd=margLik_sd)

save(mllik_out_ls,file='PT_mllik_out_ls.RData')
