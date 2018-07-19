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


source('ST_functions_new.r')
source('filledContourFunction1.r')
#burn the first half of the iterations

y=get(load('ST50000_power.RData'))$y



out_ls=get(load('ST50000_power.RData'))
########################################
#traceplots
library(coda)
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

low=15000
up=50000

############################################################################################
#calculate the theoretical mean of mu
n=25
y=out_ls$y
y_sum=sum(y)
mean_th=y_sum*(1/(n+var(y)))
sd_th=sqrt(var(y)/(n+var(y)))
mugrid=seq(-3,3,length=up-low)
thmean=rep(NA,up-low)


for (i in (1:(up-low))){
  thmean[i]=(10^(17.2))*(1/(sqrt(2*pi)^(n+1)))*exp(-(1/2)*(mugrid[i]^2))*gamma(n/2+1)*(( (1/2)*sum((y-abs(mugrid[i]))^2)+1)^(-n/2-1))
}


plot(mugrid,thmean)
############################################################################################

############################################################################################
#calculate the theoretical mean of sigma2
a=n/2+1
SSE = SSEfun(y,1.4) 
b=1/((SSE/2)+1)
sigmagrid=seq(0.0001,10,length=1000)
thsigma =rep(NA,1000)
for (i in (1:1000)){
  thsigma[i]=(10^(17.47))*(1/(sqrt(2*pi)^(n+1)))*(sigmagrid[i]^(-(n+4)/2))*exp(-1/sigmagrid[i])*exp(-(1/2)*( (sum(y^2)/sigmagrid[i]) - ( (sum(y)/sigmagrid[i])^2)/((n/sigmagrid[i])+1))     ) *sqrt(2*pi)*sqrt(1/((n/sigmagrid[i])+1))
}

plot(sigmagrid,thsigma)
############################################################################################





tau_samples=out_ls$PT_chain[[1]][low:up,2]
mu_samples=out_ls$PT_chain[[1]][low:up,1]
sigma_samples=out_ls$PT_chain[[1]][low:up,3]

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
chain1mu=density(out_ls$PT_chain[[1]][low:up,1],n=nbin)

tau_samples=c(tau_samples, (-1)*tau_samples,2-tau_samples)
denstau=density(tau_samples,from=0,to=1,n=nbin)


#geometric temperature schedule
#take kernel density estimate of the folded values of geometric temperature 
x=1:30
GeomSamples=(x/30)^5
densGeom=density(c(GeomSamples, (-1)*GeomSamples,2-GeomSamples),from=0,to=1,n=nbin)



#density estimates of mu and sig2
densSg2=density(out_ls$PT_chain[[1]][low:up,3],n=nbin)
densSg2_1=density(out_ls$PT_chain[[2]][low:up,3],n=nbin)
chain2mu=density(out_ls$PT_chain[[2]][low:up,1],adj=0.2,n=nbin)


#png('pairsplot1.png',width=1200,height=800)

setEPS()
postscript("FIG3.eps",horizontal=FALSE, paper="special",height=14,width=19, colormodel = "cmyk", 
           family = "Helvetica")



layout(matrix(c(1,2,3,4),2,2))
par(oma=c(1.5,1.5,1.5,1.5),mar=c(5,5,5,5))

#plot mu
plot(chain1mu,col='darkgray',ylim=c(0,1),xlim=c(-2,2),xlab=expression(mu),ylab='Density',main=expression(mu),
     cex.main=3.5,cex.axis=3,cex.lab=3,lwd=3,lty=1)
par(new=T)

plot(chain2mu,col=rgb(1,0,0),ylim=c(0,1),xlim=c(-2,2),xlab='',ylab='',main=expression(mu),
     cex.main=3.5,cex.axis=3,cex.lab=3,lwd=5,lty=1)

lines(mugrid,thmean,col='darkblue',lwd=3,lty=18)
lines(density(tPost), col='green',lwd=5,lty=33)


#plot sig2
plot(densSg2$x,densSg2$y,col='darkgray',xlim=c(0,10),ylim=c(0,1.5),xlab=expression(sigma^2),
     ylab='Density',main=expression(sigma^2),cex.lab=3,cex.main=3.5,cex.axis=3,lwd=3,lty=1)

par(new=T)
plot(densSg2_1$x,densSg2_1$y,col=rgb(1,0,0),xlim=c(0,10),ylim=c(0,1.5),xlab='',
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

#perspective and contour plots

############################################################################################




library(MASS)
#library(rgl)


taugrid = seq(0.000001,1,length=100)
mugrid  = seq(-4,4,length=100)
#obtain kernel density estimate of the joint posterior distribution of (mu,tau)
#by using reflected boundary condition for smoothing to avoid the edge effects of smoothing


#stack together data sets [theta, tau], [theta, - tau] and [theta, 2- tau] 
#use the reflected samples to obtain kernel density estimate
#thus avoiding the edge effects

tau_samples=out_ls$PT_chain[[1]][low:up,2]
theta_samples=out_ls$PT_chain[[1]][low:up,1]

tau_samples=c(tau_samples, (-1)*tau_samples,2-tau_samples)
theta_samples=rep(theta_samples,3)

f1=get(load('f1.RData'))
# f1=kde2d(theta_samples,tau_samples,n=c(1000,1000))
# save(f1,file='f1.RData')


 nbcol = 50
 color = terrain.colors(nbcol)
 mycut = function(x, breaks) as.numeric(cut(x=x, breaks=breaks)) 
 zcol2 = as.numeric(apply(f1$z[,which(f1$y<1 & f1$y>0)],2, mycut, breaks=nbcol)) 


#using rgl package create a perspective 3D plot of the joint posterior distribution
#of (mu,tau)
#par(mfrow=c(1,1))
persp3d(f1$x,f1$y[which(f1$y<1 & f1$y>0)],f1$z[,which(f1$y<1 & f1$y>0)],
       col = color[zcol2],xlab=expression(mu),
       ylab=expression(tau),zlab="",main=bquote(paste("a.Perspective plot of the joint posterior of" ~ mu ,'and' ~ tau,sep='')))




#codes for Greek letters in R are borrowed from http://stackoverflow.com/questions/27690729/greek-letters-symbols-and-line-breaks-inside-a-ggplot-legend-label

#install.packages("extrafont")


#jet.colors <- colorRampPalette( c("#cc0000","#ffcccc") ) 

# Generate the desired number of colors from this palette
#nbcol <- 100
#color <- jet.colors(nbcol)

nrz <- length(f1$x)
ncz <- length(f1$y[which(f1$y<1 & f1$y>0)])
zfacet <- f1$z[,which(f1$y<1 & f1$y>0)][-1, -1] + f1$z[,which(f1$y<1 & f1$y>0)][-1, -ncz] + f1$z[,which(f1$y<1 & f1$y>0)][-nrz, -1] + f1$z[,which(f1$y<1 & f1$y>0)][-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)


library(extrafont)
font_import()

#create a contour plot of the respective perspective plot
#pdf('Persp1.pdf',height = 9,width=9,family = 'serif')
# png('Persp1.png',width=550,height = 484)
# persp(f1$x,f1$y[which(f1$y<1 & f1$y>0)],f1$z[,which(f1$y<1 & f1$y>0)],expand=0.8,col= color[zcol2],
#       xlab="\u03BC",ltheta = 130, shade = 0.4,border=NA,mgp=c(3,3,3),
#       ylab="\n \u03c4",family = 'serif',ticktype='detailed',theta=-150,phi=0,zlab="",
#       main=bquote(paste("A. Perspective plot of the joint posterior of " ~ mu ,' and ' ~ tau,sep='')),cex.axis=2,cex.lab=2,cex.main=2)
# dev.off()

#for the presentation
#png('Persp1.png',width=550,height = 484)
setEPS()
postscript("FIG2.eps",horizontal=FALSE, paper="special",height=10,width=20, colormodel = "cmyk",
           family='Helvetica')

par(mfrow=c(1,2))
#windows(family='serif')
persp(f1$x,f1$y[which(f1$y<1 & f1$y>0)],f1$z[,which(f1$y<1 & f1$y>0)],expand=0.8,col= color[facetcol],
      xlab="\u03BC",ltheta = 8, shade = 0.6,border=NA,mgp=c(3,3,3),
      ylab="\u03c4",ticktype='detailed',theta=-150,phi=0,zlab="",
      main=bquote(paste("A. Perspective plot of the joint posterior of " ~ mu ,' and ' ~ tau,sep='')),cex.axis=2,cex.lab=2,cex.main=2.5)

# persp(f1$x,f1$y[which(f1$y<1 & f1$y>0)],f1$z[,which(f1$y<1 & f1$y>0)],expand=0.8,col= color[facetcol],
#       xlab="\u03BC",ltheta = 8, shade = 0.6,border=NA,mgp=c(3,3,3),
#       ylab="\n \u03c4",ticktype='detailed',theta=-150,phi=0,zlab="",
#       main=bquote(paste("A. Perspective plot of the joint posterior of " ~ mu ,' and ' ~ tau,sep='')),cex.axis=2,cex.lab=2,cex.main=2.5)

#dev.off()


#windows(family="Helvetica")
#png('Contourplot.png',width=550,height = 484)
#pdf('Contourplot.pdf',height = 9,width=9,points=12,family = 'serif')
#filled.contour4(f1$x,f1$y[which(f1$y<1 & f1$y>0)],f1$z[,which(f1$y<1 & f1$y>0)],col=terrain.colors(10),nlevel=7,axes=F,xlim=c(-2.5,2.5),xlab=expression(mu),ylab=expression(tau),cex.axis=2,cex.lab=2,cex.main=2.5,family = 'serif',main=paste0('B. Contours of the joint posterior of \u03BC and \u03c4',sep=''))#, expression(mu),' and ',expression(tau),sep=''))
filled.contour4(f1$x,f1$y[which(f1$y<1 & f1$y>0)],f1$z[,which(f1$y<1 & f1$y>0)],col=terrain.colors(10),nlevel=7,axes=F,xlim=c(-2.5,2.5),xlab=expression(mu),ylab=expression(tau),cex.axis=2,cex.lab=2,cex.main=2.5,main=bquote(paste("B. Contours of the joint posterior of " ~ mu ,' and ' ~ tau,sep='')))#, expression(mu),' and ',expression(tau),sep=''))
#filled.contour4(f1$x,f1$y[which(f1$y<1 & f1$y>0)],f1$z[,which(f1$y<1 & f1$y>0)],col=terrain.colors(10),nlevel=7,axes=F,xlim=c(-2.5,2.5),xlab=expression(mu),ylab=expression(tau),cex.axis=2,cex.lab=2,cex.main=2.5,main=paste0('B. Contours of the joint posterior of \u03BC and \u03c4',sep=''))#, expression(mu),' and ',expression(tau),sep=''))
contour(f1$x,f1$y[which(f1$y<1 & f1$y>0)],f1$z[,which(f1$y<1 & f1$y>0)],add=TRUE,nlevel=6,labcex=1)
axis(2,seq(0,1,length=6),cex.axis=2,at=seq(0,1,length=6))
axis(1,seq(-2,2,length=5),cex.axis=2,at=seq(-2,2,length=5))
dev.off()


break()
############################################################################################
#make some animations of the contour plot above
#save the animated plots in the subfolder Animation_new
############################################################################################
setwd('D:/PhdProjects/PT_STWNC/SubmitCode/SubmitCode/Example1')
animations=1:400 
for (i in animations){
  png(paste('Animation_new/Animation',i,'.png',sep=''),width=550,height = 484)
  
  filled.contour4(f1$x,f1$y[which(f1$y<1 & f1$y>0)],f1$z[,which(f1$y<1 & f1$y>0)],col=terrain.colors(10),nlevel=7,axes=F,xlim=c(-2.5,2.5),xlab=expression(mu),ylab=expression(tau),cex.axis=2,cex.lab=2,cex.main=2,family = 'serif',main=paste0('Contours of the joint posterior of \u03BC and \u03c4',sep=''))#, expression(mu),' and ',expression(tau),sep=''))
  axis(2,seq(0,1,length=6),cex.axis=2,at=seq(0,1,length=6))
  axis(1,seq(-2,2,length=5),cex.axis=2,at=seq(-2,2,length=5))
  points(out_ls$PT_chain[[1]][1,1],out_ls$PT_chain[[1]][1,2])
  
  for (j in (1:i)){
  points(out_ls$PT_chain[[1]][j,1],out_ls$PT_chain[[1]][j,2])
  lines(out_ls$PT_chain[[1]][1:j,1],out_ls$PT_chain[[1]][1:j,2])
  }
  
  dev.off()
}
library(animation)
#to create a .gif animation type
#system("convert  *.png PT-STWNC.gif")
#in the cmd console directly


############################################################################################
#make new plots os sample paths
#sample path plot of the two PT-STWTDNC chains
#'tempered' chain is in 'black'
#'target' chain is in 'red'
############################################################################################
#library(coda)
pdf('smplpath.pdf',width=10,height=6)

par(mfrow=c(1,1))
chain=1
plot(out_ls$PT_chain[[chain]][low:up,1],xlab='iterations',ylab='',main='Sampled paths of the two chains')
chain=2
points(out_ls$PT_chain[[chain]][low:up,1],col='red')
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
pos_mean_se=sqrt(var(out_ls$PT_chain[[chain]][low:up,1][which(out_ls$PT_chain[[chain]][low:up,1]>0)])/(up-low))*sqrt(acfmean)

#posterior mean of mu for the mode on the negative side
(neg_mean=mean(out_ls$PT_chain[[chain]][low:up,1][which(out_ls$PT_chain[[chain]][low:up,1]<0)]))

acfmean=acf(out_ls$PT_chain[[chain]][low:up,1][which(out_ls$PT_chain[[chain]][low:up,1]<0)])$acf
acfmean[which(acfmean<0.05)]=0
acfmean=2*sum(acfmean)-1
#standard error adjusted for autocorrelation of the posterior mean estimate of mu (the negative side)
neg_mean_se=sqrt(var(out_ls$PT_chain[[chain]][low:up,1][which(out_ls$PT_chain[[chain]][low:up,1]<0)])/(up-low))*sqrt(acfmean)


chain=2
#posterior mean of sigma2
mean(out_ls$PT_chain[[chain]][low:up,3])

acfsigma=acf(out_ls$PT_chain[[chain]][low:up,3])$acf
acfsigma[which(acfsigma<0.05)]=0
acfsigma=2*sum(acfsigma)-1
#standard error adjusted for autocorrelation of the posterior mean estimate of sigma2
sigma2_se=sqrt(var(out_ls$PT_chain[[chain]][low:up,1])/(up-low))*sqrt(acfsigma)


#marginal likelihood and acceptance rates
chain=1
mean(out_ls$mllik[[chain]][low:up])
#-40.43459 from chain=1

acc_mu_mu_chain1=sum(out_ls$accepts[[1]][,1]/up)
#0.43464
acc_mu_sigma_chain1=sum(out_ls$accepts[[1]][,2]/up)
#0.51374
acc_tau=(out_ls$accepts1[1]/up)
#0.69096

acc_mu_mu_chain2=sum(out_ls$accepts[[2]][,1]/up)
#0.3882
acc_mu_sigma_chain2=sum(out_ls$accepts[[2]][,2]/up)
#0.49304



#plot of the geometric schedule
#ti=(i/M)^5
# ti=sapply(1:30, function(x) {(x/30)^5})
# pdf('Geom_schedule.pdf')
# hist(ti,nclass=50,col=rgb(1,0,0,0.5),prob=T,main='Geometric schedule',xlab=expression(tau))
# dev.off()