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

source('../Functions/ST_functions.R')
source("ST_pants.r")
library(truncnorm)
library("MASS")
library(coda)
library(MCMCpack)
#library(LaplacesDemon)


#generate data
n=25
p_mu=1.5
sigma=1
set.seed(2555)

y=data(n=n,mu=abs(p_mu),sigma=sigma)


out_ls=runST(niter=50000,y=y,PriorPars=c(0,1,1,1),IniPar=c(-0.8,0.2,0.5),tune_pars_init=c(0.15,0.1),ttune_pars=c(rep(TRUE,2)))


