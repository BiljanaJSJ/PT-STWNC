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

#setwd('D:/PhdProjects/PT_STWNC/SubmitCode/SubmitCode/Example2')



out_ls=get(load('ST35000.RData'))


#make a pairs plot from the paper
low=1000
up=35000

length(which(out_ls$swappers[,3]==1))/up

library(MASS)

library(coda)
########################################
pdf('traceplots.pdf')
par(mfrow=c(3,2))
traceplot(as.mcmc(out_ls$PT_chain[[1]][,1]))
traceplot(as.mcmc(out_ls$PT_chain[[2]][,1]))
traceplot(as.mcmc(out_ls$PT_chain[[1]][,2]))
traceplot(as.mcmc(out_ls$PT_chain[[2]][,2]))
traceplot(as.mcmc(out_ls$PT_chain[[1]][,3]))
traceplot(as.mcmc(out_ls$PT_chain[[2]][,3]))

traceplot(as.mcmc(out_ls$PT_chain[[1]][,4]))
traceplot(as.mcmc(out_ls$PT_chain[[2]][,4]))
traceplot(as.mcmc(out_ls$PT_chain[[1]][,11]))
traceplot(as.mcmc(out_ls$PT_chain[[2]][,11]))
dev.off()


#kernel density estimates of mu1 and mu2 from tempered and target chains
mu_chain1=density(out_ls$PT_chain[[1]][low:up,1],n=1000000)
mu_chain2=density(out_ls$PT_chain[[2]][low:up,1],n=1000000)
mu2_chain1=density(out_ls$PT_chain[[1]][low:up,2],n=1000000)
mu2_chain2=density(out_ls$PT_chain[[2]][low:up,2],n=1000000)
mu3_chain1=density(out_ls$PT_chain[[1]][low:up,3],n=1000000)
mu3_chain2=density(out_ls$PT_chain[[2]][low:up,3],n=1000000)

#kernel density estimates of folded tau samples 
tau_samples=out_ls$PT_chain[[1]][low:up,11]
tau_samples=c(tau_samples, (-1)*tau_samples,2-tau_samples)
denstau=density(tau_samples,n=1000000)

#kernel density estimates of sigma1 from tempered and target chains
sigma2_chain1=density(out_ls$PT_chain[[1]][low:up,4],n=100000)
sigma2_chain2=density(out_ls$PT_chain[[2]][low:up,4],n=100000)

#make a plot
#png('GalaxyDataPlots.png',width=1200,height=800)

setEPS()
postscript("FIG4.eps",horizontal=FALSE, paper="special",height=13,width=19, colormodel = "cmyk", 
					 family = "Helvetica")


par(mgp=c(5,1,0))

layout(matrix(c(1,2,3,4,5,6),2,3))
par(oma=c(4,4,rep(4,2))+0.05,mar=c(7,7,7,7))

plot(mu_chain2,col=rgb(1,0,0),ylim=c(0,0.3),main=expression(mu [1]),
		 cex.main=3.5,cex.axis=2,lwd=5,lty=1,ylab='Density',xlab=expression(mu [1]),cex.lab=3)
lines(mu_chain1,col='darkgrey',lty=1,lwd=3)

bivn.kde <- kde2d(out_ls$PT_chain[[1]][low:up,1], out_ls$PT_chain[[1]][low:up,3])

plot(x = 0, y = 0, type = "n",cex.axis=2, xlim = range(bivn.kde$x), ylim = range(bivn.kde$y), xlab=expression(mu [3]),ylab=expression(mu [3]),main=expression(paste(mu [1],' and ',mu [3],sep='')),cex.main=3.5,cex.lab=3)
contour(bivn.kde,col='darkgrey',add=T,axes=F,lwd=3,lty=1,labels='')

bivn.kde_target <- kde2d(out_ls$PT_chain[[2]][low:up,1], out_ls$PT_chain[[2]][low:up,3],n=20)
contour(bivn.kde_target,col=rgb(1,0,0),add=TRUE,axes=F,lwd=5,labels='',nlevels = 80,levels = pretty(range(bivn.kde_target$z), 70))


plot(mu3_chain2,col=rgb(1,0,0),ylim=c(0,0.3),main=expression(mu [3]),
		 cex.main=3.5,cex.axis=2,lwd=5,lty=1,ylab='Density',xlab=expression(mu [3]),cex.lab=3)
lines(mu3_chain1,col='darkgrey',lty=1,lwd=3)


plot(sigma2_chain2,col=rgb(1,0,0),main=expression(sigma [1]^2),
		 cex.main=3.5,cex.axis=2,lwd=5,lty=1,ylab='Density',xlab=expression(sigma [1]^2),cex.lab=3)

lines(sigma2_chain2,col='darkgrey',lty=1,lwd=3)


par(mar=rep(4,4))
plot(1, type="n", axes=F, xlab="", ylab="")

legend("topright", 
			 c("'tempered' chain","'target' chain"), 
			 lty=c(1,1), cex=3,
			 lwd=c(3,5),col=c("darkgrey",rgb(1,0,0))) 
par(mar=rep(7,4))

plot(denstau$x[denstau$x>0 & denstau$x<1],denstau$y[denstau$x>0 & denstau$x<1],ylim=c(0,1.7),col='darkgrey',cex.axis=2,main=expression(tau),ylab='Density',xlab=expression(tau),cex.lab=3,cex.main=3.5,lwd=3,lty=1)

dev.off()

