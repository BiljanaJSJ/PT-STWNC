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


#####################################################################################
#perform diagnostics on the chains
#####################################################################################

ch30=ch60=ch100=list()
library(coda)

for (i in (1:20)){
ch30[[i]]=get(load(paste('chains30/TI_mllik', i,'.RData',sep='')))
ch60[[i]]=get(load(paste('chains60/TI_mllik', i,'.RData',sep='')))
ch100[[i]]=get(load(paste('chains100/TI_mllik', i,'.RData',sep='')))
}



p_est=raft=gelman=list()

# 30 chains
p_est[[1]]=gelman[[1]]=raft[[1]]=matrix(NA,ncol=2,nrow=20)

for (i in (1:20)){

multiple=do.call(cbind,ch30[[i]]$mu)
chains=lapply(1:ncol(multiple),function(x) {multiple[,x]=as.mcmc(multiple[,x])})

chains=mcmc.list(chains)

grdg=gelman.diag(chains, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
            multivariate=TRUE)

gelman[[1]][i,]=grdg$psrf 

#the R is not bigger than 1, and the CI is very small
#plot(chains)


raftdiag=raftery.diag(chains)
raft[[1]][i,]=raftdiag[[2]]$resmatrix[,c(1,4)]

#no dependence factor >5 hence no autocorrelation

p_est_pos=length(ch30[[i]]$mu[[30]][which(ch30[[i]]$mu[[30]]>0)])
p_est_neg=length(ch30[[i]]$mu[[30]][which(ch30[[i]]$mu[[30]]<0)])
p_est[[1]][i,]=c(p_est_neg,p_est_pos)/(p_est_neg+p_est_pos)
}



# 60 chains
p_est[[2]]=gelman[[2]]=raft[[2]]=matrix(NA,ncol=2,nrow=20)

for (i in (1:20)){

multiple=do.call(cbind,ch60[[i]]$mu)
chains=lapply(1:ncol(multiple),function(x) {multiple[,x]=as.mcmc(multiple[,x])})

chains=mcmc.list(chains)

grdg=gelman.diag(chains, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
            multivariate=TRUE)

gelman[[2]][i,]=grdg$psrf 

#the R is not bigger than 1, and the CI is very small
#plot(chains)


raftdiag=raftery.diag(chains)
raft[[2]][i,]=raftdiag[[2]]$resmatrix[,c(1,4)]

#no dependence factor >5 hence no autocorrelation

p_est_pos=length(ch60[[i]]$mu[[60]][which(ch60[[i]]$mu[[60]]>0)])
p_est_neg=length(ch60[[i]]$mu[[60]][which(ch60[[i]]$mu[[60]]<0)])
p_est[[2]][i,]=c(p_est_neg,p_est_pos)/(p_est_neg+p_est_pos)
}


# 100 chains
p_est[[3]]=gelman[[3]]=raft[[3]]=matrix(NA,ncol=2,nrow=20)

for (i in (1:20)){

multiple=do.call(cbind,ch100[[i]]$mu)
chains=lapply(1:ncol(multiple),function(x) {multiple[,x]=as.mcmc(multiple[,x])})

chains=mcmc.list(chains)

grdg=gelman.diag(chains, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
            multivariate=TRUE)

gelman[[3]][i,]=grdg$psrf 

#the R is not bigger than 1, and the CI is very small
#plot(chains)


raftdiag=raftery.diag(chains)
raft[[3]][i,]=raftdiag[[2]]$resmatrix[,c(1,4)]

#no dependence factor >5 hence no autocorrelation

p_est_pos=length(ch100[[i]]$mu[[100]][which(ch100[[i]]$mu[[100]]>0)])
p_est_neg=length(ch100[[i]]$mu[[100]][which(ch100[[i]]$mu[[100]]<0)])
p_est[[3]][i,]=c(p_est_neg,p_est_pos)/(p_est_neg+p_est_pos)
}


chain_diagnostics=list(p_est=p_est,gelman=gelman,raft=raft,chains=c('30','60','100'))
save(chain_diagnostics,file='chain_diagnostics.RData')


#100 chains
range(abs(chain_diagnostics[[1]][[3]][,1]-chain_diagnostics[[1]][[3]][,2])*100)
#  0.5257143 19.6800000
#60 chains
range(abs(chain_diagnostics[[1]][[2]][,1]-chain_diagnostics[[1]][[2]][,2])*100)
#0.04571429 13.05142857
#30 chains
range(abs(chain_diagnostics[[1]][[1]][,1]-chain_diagnostics[[1]][[1]][,2])*100)
# 0.5257143 9.0400000

