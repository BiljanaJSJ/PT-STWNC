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
source('../Functions/ST_functions.r')
source('filledContourFunction1.r')
source('st_pants.R')
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

low=15000
up=50000

############################################################################################



############################################################################################

#perspective and contour plots

############################################################################################




library(MASS)
library(rgl)


taugrid = seq(0.000001,1,length=100)
mugrid  = seq(-4,4,length=100)
#obtain kernel density estimate of the joint posterior distribution of (mu,tau)
#by using reflected boundary condition for smoothing to avoid the edge effects of smoothing


#stack together data sets [theta, tau], [theta, - tau] and [theta, 2- tau] 
#use the reflected samples to obtain kernel density estimate
#thus avoiding the edge effects

tau_samples=out_ls$PT_chain[[1]][low:up,3]
theta_samples=out_ls$PT_chain[[1]][low:up,1]

tau_samples=c(tau_samples, (-1)*tau_samples,2-tau_samples)
theta_samples=rep(theta_samples,3)

#f1=get(load('f1.RData'))
f1=kde2d(theta_samples,tau_samples,n=c(1000,1000))
save(f1,file='f1.RData')


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



#for the presentation
#png('Persp1.png',width=550,height = 484)
setEPS()
postscript("FIG2.eps",horizontal=FALSE, paper="special",height=10,width=20, colormodel = "cmyk",
           family='Helvetica')

par(mfrow=c(1,2))
#windows(family='serif')
persp(f1$x,f1$y[which(f1$y<1 & f1$y>0)],f1$z[,which(f1$y<1 & f1$y>0)],expand=0.8,col= color[facetcol],
      xlab="\u03BC",ltheta = 8, shade = 0.6,border=NA,mgp=c(3,3,3),
      ylab="\u03c4",ticktype='detailed',theta=140,phi=0,zlab="",
      main=bquote(paste("A. Perspective plot of the joint posterior of " ~ mu ,' and ' ~ tau,sep='')),cex.axis=2,cex.lab=2,cex.main=2.5)


#filled.contour4(f1$x,f1$y[which(f1$y<1 & f1$y>0)],f1$z[,which(f1$y<1 & f1$y>0)],col=terrain.colors(10),nlevel=7,axes=F,xlim=c(-2.5,2.5),xlab=expression(mu),ylab=expression(tau),cex.axis=2,cex.lab=2,cex.main=2.5,family = 'serif',main=paste0('B. Contours of the joint posterior of \u03BC and \u03c4',sep=''))#, expression(mu),' and ',expression(tau),sep=''))
filled.contour4(f1$x,f1$y[which(f1$y<1 & f1$y>0)],f1$z[,which(f1$y<1 & f1$y>0)],col=terrain.colors(10),nlevel=7,axes=F,xlim=c(-2.5,2.5),xlab=expression(mu),ylab=expression(tau),cex.axis=2,cex.lab=2,cex.main=2.5,main=bquote(paste("B. Contours of the joint posterior of " ~ mu ,' and ' ~ tau,sep='')))#, expression(mu),' and ',expression(tau),sep=''))
#filled.contour4(f1$x,f1$y[which(f1$y<1 & f1$y>0)],f1$z[,which(f1$y<1 & f1$y>0)],col=terrain.colors(10),nlevel=7,axes=F,xlim=c(-2.5,2.5),xlab=expression(mu),ylab=expression(tau),cex.axis=2,cex.lab=2,cex.main=2.5,main=paste0('B. Contours of the joint posterior of \u03BC and \u03c4',sep=''))#, expression(mu),' and ',expression(tau),sep=''))
contour(f1$x,f1$y[which(f1$y<1 & f1$y>0)],f1$z[,which(f1$y<1 & f1$y>0)],add=TRUE,nlevel=6,labcex=1)
axis(2,seq(0,1,length=6),cex.axis=2,at=seq(0,1,length=6))
axis(1,seq(-2,2,length=5),cex.axis=2,at=seq(-2,2,length=5))
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


