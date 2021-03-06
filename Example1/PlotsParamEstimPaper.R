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


#this script creates the plots from the paper

library(coda)
source('../Functions/ST_functions.R')
source('filledContourFunction1.R')
source('ST_pants.r')
#burn the first half of the iterations

out_ls= get(load('ST50000.RData'))
y     = out_ls$y
length(which(out_ls$swappers[,3]==1))

########################################
#traceplots

#library(ks)

########################################
par(mfrow=c(3,2))
traceplot(as.mcmc(out_ls$PT_chain[[1]][,1]))
traceplot(as.mcmc(out_ls$PT_chain[[2]][,1]))
traceplot(as.mcmc(out_ls$PT_chain[[1]][,2]))
traceplot(as.mcmc(out_ls$PT_chain[[2]][,2]))
traceplot(as.mcmc(out_ls$PT_chain[[1]][,3]))
traceplot(as.mcmc(out_ls$PT_chain[[2]][,3]))
############################################################################################

low=1000
up=50000

############################################################################################
#calculate the theoretical mean of mu
n=length(y)
#y=out_ls$y
y_sum=sum(y)
mean_th=y_sum*(1/(n+var(y)))
sd_th=sqrt(var(y)/(n+var(y)))
mugrid=seq(-3,3,length=up-low)
thmean=rep(NA,up-low)

for (i in (1:(up-low))){
  thmean[i]=(10^(17.177))*(1/(sqrt(2*pi)^(n+1)))*exp(-(1/2)*(mugrid[i]^2))*gamma(n/2+1)*(( (1/2)*sum((y-abs(mugrid[i]))^2)+1)^(-n/2-1))
}


plot(mugrid,thmean)
############################################################################################
# plot(chain2mu,col=rgb(1,0,0),ylim=c(0,1.2),xlim=c(-2,2),xlab='',ylab='',main=expression(mu),
# 		 cex.main=3.5,cex.axis=3,cex.lab=3,lwd=6.5,lty=1)
# 
# lines(mugrid,thmean,col='darkblue',lwd=3,lty=18)
# 


############################################################################################
#calculate the theoretical mean of sigma2
a=n/2+1
SSE = SSEfun(y,1.4) 
b=1/((SSE/2)+1)
sigmagrid=seq(0.0001,10,length=1000)
thsigma =rep(NA,1000)
for (i in (1:1000)){
  thsigma[i]=(10^(17.56))*(1/(sqrt(2*pi)^(n+1)))*(sigmagrid[i]^(-(n+4)/2))*exp(-1/sigmagrid[i])*exp(-(1/2)*( (sum(y^2)/sigmagrid[i]) - ( (sum(y)/sigmagrid[i])^2)/((n/sigmagrid[i])+1))     ) *sqrt(2*pi)*sqrt(1/((n/sigmagrid[i])+1))
}

plot(sigmagrid,thsigma)
############################################################################################

# plot(densSg2_1$x,densSg2_1$y,col=rgb(1,0,0),xlim=c(0,10),ylim=c(0,1.7),xlab='',
# 		 ylab='',main=expression(sigma^2),cex.main=3.5,cex.axis=3,cex.lab=3,xaxt='n',lwd=5,lty=1)
# lines(sigmagrid,thsigma,col='darkblue',lwd=3,lty=18)




tau_samples=out_ls$PT_chain[[1]][low:up,3]
mu_samples=out_ls$PT_chain[[1]][low:up,1]
sigma_samples=out_ls$PT_chain[[1]][low:up,2]

#fold the samples 
tau_samples=c(tau_samples, (-1)*tau_samples,2-tau_samples)
mu_samples=rep(mu_samples,3)
sigma_samples=rep(sigma_samples,3)

#kernel density estimate of the folded samples
#z=kde2d(mu_samples,tau_samples)

#get results from the ST-PAWN 
result_adapt_500=get(load('WL_pants/result_adapt_500.RData'))

M=10
(N_vec = 273.45*2^(0:7))
N=N_vec[8]
targetPost=list()
X_sample_smpled=result_adapt_500$X_sample[1000:N,]
T_samples_sampled=result_adapt_500$T_sample[1000:N,]
# for (i in (1:M)){
#   targetPost[[i]]=X_sample_smpled[which(1/T_samples_sampled[,1]==1),i]
#   
# }
tPost=do.call(rbind,targetPost)
tPost=as.vector(X_sample_smpled[1/T_samples_sampled==1])



#marginal posteriors of all parameters from tempered and target chain
nbin=1000000

tau_samples=c(tau_samples, (-1)*tau_samples,2-tau_samples)
denstau=density(tau_samples,from=0,to=1,n=nbin)


#geometric temperature schedule
#take kernel density estimate of the folded values of geometric temperature 
x=1:30
GeomSamples=(x/30)^5
densGeom=density(c(GeomSamples, (-1)*GeomSamples,2-GeomSamples),from=0,to=1,n=nbin)



#density estimates of mu and sig2
densSg2=density(out_ls$PT_chain[[1]][low:up,2],adj=1.5,n=nbin)
densSg2_1=density(out_ls$PT_chain[[2]][low:up,2],adj=1.2,n=nbin)
chain2mu=density(out_ls$PT_chain[[2]][low:up,1],adj=0.5,n=nbin)
chain1mu=density(out_ls$PT_chain[[1]][low:up,1],adj=1.15,n=nbin)


#png('pairsplot1.png',width=1200,height=800)

setEPS()
postscript("FIG3.eps",horizontal=FALSE, paper="special",height=14,width=19, colormodel = "cmyk", 
           family = "Helvetica")



layout(matrix(c(1,2,3,4),2,2))
par(oma=c(1.5,1.5,1.5,1.5),mar=c(5,5,5,5))

#plot mu
plot(chain1mu,col='darkgray',ylim=c(0,1.2),xlim=c(-2,2),xlab=expression(mu),ylab='Density',main=expression(mu),
     cex.main=3.5,cex.axis=3,cex.lab=3,lwd=3,lty=1)
par(new=T)

plot(chain2mu,col=rgb(1,0,0),ylim=c(0,1.2),xlim=c(-2,2),xlab='',ylab='',main=expression(mu),
     cex.main=3.5,cex.axis=3,cex.lab=3,lwd=6.5,lty=1)

lines(mugrid,thmean,col='darkblue',lwd=3,lty=18)
lines(density(tPost), col='green',lwd=5,lty=33)


#plot sig2
plot(densSg2$x,densSg2$y,col='darkgray',xlim=c(0,10),ylim=c(0,1.7),xlab=expression(sigma^2),
     ylab='Density',main=expression(sigma^2),cex.lab=3,cex.main=3.5,cex.axis=3,lwd=3,lty=1)

par(new=T)
plot(densSg2_1$x,densSg2_1$y,col=rgb(1,0,0),xlim=c(0,10),ylim=c(0,1.7),xlab='',
     ylab='',main=expression(sigma^2),cex.main=3.5,cex.axis=3,cex.lab=3,xaxt='n',lwd=5,lty=1)
lines(sigmagrid,thsigma,col='darkblue',lwd=3,lty=18)




#add legend
plot(1, type="n", axes=F, xlab="", ylab="")
par(mar=rep(0.003,4))
legend('center', 
       c("'tempered' chain","'target' chain","theoretical distribution","PAWL within ST"), 
       lty=c(1,1,18,33), cex=3.5,
       lwd=c(3,5,3,4),col=c("darkgray",rgb(1,0,0),"darkblue",'green')) 




#tau
par(mar=rep(5,4))
#h=hist(1/T_samples_sampled[,1],plot=F,nclass=35)
h=hist(as.vector(1/T_samples_sampled),plot=F,nclass=35)
h$counts<-(h$counts/sum(h$counts))


#plot(denstau$x[which(denstau$x>0 & denstau$x<1)],denstau$y[which(denstau$x>0 & denstau$x<1)],lwd=3,xlim=c(0,1),col='darkgray',xlab=expression(tau),ylab='Density',main=expression(tau),lty=1,cex.main=3.5,cex.axis=2.5,cex.lab=2.5,type='l')
plot(denstau$x,denstau$y,lwd=3,xlim=c(0,1),col='darkgray',xlab=expression(tau),ylab='Density',main=expression(tau),lty=1,cex.main=3.5,cex.axis=3,cex.lab=3,type='l')
legend('topright', 
       c('Geometric schedule'), 
       lty=34, cex=3,
       lwd=2,col='purple') 


par(new=T)
#add discrete  tau from PAWN
plot(h,col='green',xlab='',xlim=c(0,1),main='',xaxt='n',yaxt='n',lty=33,lwd=4,ylab='')
#add discrete geometric schedule 
#lines(densGeom$x[which(densGeom$x>0 & densGeom$x<1)],10^-0.5*densGeom$y[which(densGeom$x>0 & densGeom$x<1)],xlim=c(0,1),col='purple',lty=34,lwd=2,xlab='',ylab='',main='')
lines(densGeom$x,10^-0.5*densGeom$y,xlim=c(0,1),col='purple',lty=34,lwd=2,xlab='',ylab='',main='')


dev.off()



############################################################################################

#obtain parameter estimates of mu and sigma2 from the 'target' chain

############################################################################################
chain=2
#posterior mean of mu for the mode on the positive side
(pos_mean=mean(out_ls$PT_chain[[chain]][low:up,1][which(out_ls$PT_chain[[chain]][low:up,1]>0)]))

acfmean=acf(out_ls$PT_chain[[chain]][low:up,1][which(out_ls$PT_chain[[chain]][low:up,1]>0)])$acf
acfmean[which(acfmean<0.05)]=0
acfmean=2*sum(acfmean)-1
#standard error adjusted for autocorrelation of the posterior mean estimate of mu (the positive side)
(pos_mean_se=sqrt(var(out_ls$PT_chain[[chain]][low:up,1][which(out_ls$PT_chain[[chain]][low:up,1]>0)])/(up-low))*sqrt(acfmean))

#posterior mean of mu for the mode on the negative side
(neg_mean=mean(out_ls$PT_chain[[chain]][low:up,1][which(out_ls$PT_chain[[chain]][low:up,1]<0)]))

acfmean=acf(out_ls$PT_chain[[chain]][low:up,1][which(out_ls$PT_chain[[chain]][low:up,1]<0)])$acf
acfmean[which(acfmean<0.05)]=0
acfmean=2*sum(acfmean)-1
#standard error adjusted for autocorrelation of the posterior mean estimate of mu (the negative side)
(neg_mean_se=sqrt(var(out_ls$PT_chain[[chain]][low:up,1][which(out_ls$PT_chain[[chain]][low:up,1]<0)])/(up-low))*sqrt(acfmean))


chain=2
#posterior mean of sigma2
mean(out_ls$PT_chain[[chain]][low:up,2])

acfsigma=acf(out_ls$PT_chain[[chain]][low:up,2])$acf
acfsigma[which(acfsigma<0.05)]=0
acfsigma=2*sum(acfsigma)-1
#standard error adjusted for autocorrelation of the posterior mean estimate of sigma2
(sigma2_se=sqrt(var(out_ls$PT_chain[[chain]][low:up,1])/(up-low))*sqrt(acfsigma))


#marginal likelihood and acceptance rates
chain=1
mean(out_ls$mllik[[chain]][low:up])

mean(c(out_ls$mllik[[1]][low:up],out_ls$mllik[[2]][low:up]))

acc_mu_mu_chain1=sum(out_ls$accepts[[1]][,1]/up)
acc_mu_sigma_chain1=sum(out_ls$accepts[[1]][,2]/up)
acc_tau=sum(out_ls$acceptsTau/up)
acc_mu_mu_chain2=sum(out_ls$accepts[[2]][,1]/up)
acc_mu_sigma_chain2=sum(out_ls$accepts[[2]][,2]/up)

#calculate posterior mean for mu
n=length(y)
y_bar           = mean(y)
n*y_bar*((1*1^2)/(n*1*1^2+1))

#calculate posterior mean for sigma2
SSE                = SSEfun(y,1.409704) 
d0                 = n/2+1
v0                 = SSE/2+1
(postmean          = v0/(d0-1))

