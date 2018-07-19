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


####################################################################
#generate plots
####################################################################


setwd('D:/Publications/PT-STWNC/Code/Example4')
source('SIR_ST_functions.R')
source('pairs2_function.R')
out_ls=get(load('ST35000.RData'))
low=1000
up=35000

## put histograms on the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  #par(usr = c(0,1, 0, 1))
  if (unique(x) %in% c(1,2,3,4,5,6,7)){
    h <- hist(x, plot = FALSE,breaks=seq(range(x)[1],range(x)[2],by=0.2))  
  }else{
  h <- hist(x, plot = FALSE,breaks=150)
  }
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}


## put histograms on the diagonal
#pdf('pairsplot_SIR1.pdf',width=15,height=12,pointsize=15)
#png('pairsplot_SIR1.png',width = 1200,height = 700)
setEPS()
postscript("FIG6.eps",horizontal=FALSE, paper="special",height=12,width=19, colormodel = "cmyk", 
           family = "Helvetica")

par(mfrow=c(2,2))
hist(out_ls$PT_chain[[1]][low:up,2],main=expression(alpha),cex.axis=2.5,cex.main=3.5,nclass=150,col='cyan', probability = T, xlab   ='', ylab   ='')
hist(out_ls$PT_chain[[1]][low:up,1],main=expression(beta),cex.axis=2.5,cex.main=3.5,nclass=150,col='cyan', probability = T, xlab   ='', ylab   ='')
hist(out_ls$PT_chain[[1]][low:up,4],main=expression(I(0)),cex.axis=2.5,cex.main=3.5,nclass=150,col='cyan', probability = T, xlab   ='', ylab   ='')
hist(out_ls$tau[low:up],main=expression(tau),cex.axis=2.5,cex.main=3.5,nclass=150,col='cyan', probability = T, xlab   ='', ylab   ='')
dev.off()

#pdf('pairsplot_SIR.pdf',width=15,height=12,pointsize=15)
#png('pairsplot_SIR.png',width = 1200,height = 700)
setEPS()
postscript("FIG5.eps",horizontal=FALSE, paper="special",height=15,width=19, colormodel = "cmyk", 
           family = "Helvetica")
pairs2(out_ls$PT_chain[[2]][low:up,c(2,1,4)],labels=c(expression(alpha),expression(beta),expression(I(0)),expression(tau)), bg="cyan",
      diag.panel=panel.hist,cex.labels = 4, font.labels=4,cex.axis=4,upper.panel=NULL, oma=c(2,3,2,2),gap=8,label.pos =0.9,mgp = c(5, 3, 1))
dev.off()






library(coda)
traceplot(as.mcmc(out_ls$PT_chain[[chain]][,1]))
traceplot(as.mcmc(out_ls$PT_chain[[chain]][,2]))
traceplot(as.mcmc(out_ls$PT_chain[[chain]][,4]))


par(mfrow=c(2,1))
traceplot(as.mcmc(out_ls$tau))
traceplot(as.mcmc(out_ls$PT_chain[[1]][,4]))
str(out_ls)


  chain=2
   #plot of joint prior of alpha and beta compared to the joint posterior 
   #distribution of alpha and beta

  # pdf('PriorPosteriorab.pdf',width = 9, height = 8)
   #png('PriorPosteriorab.png',width=650,height=550)
   setEPS()
   postscript("FIG7.eps",horizontal=FALSE, paper="special",height=10,width=19, colormodel = "cmyk", 
              family = "Helvetica")
   par(mfrow=c(1,1), mar=c(4.5,6.4, 5.1, 5.1))
   #oldpar=par(mfrow=c(1,2), mar=c(6.5,5.4, 4.1, 4.1))
   f1= kde2d(rgamma(500,shape=1,scale=1),rgamma(500,shape=1,scale=1))
   #par(mfrow=c(1,1), oma=c(1,1,1,1))
   contour(f1$x,f1$y,f1$z,col='dimgrey',cex.axis=3)

   legend('topright', 
          c("prior","posterior"), inset=c(-0.0002,0),
          lty=c(1,1), cex=3,
          lwd=c(2,2),col=c("lightgrey","red")) 
   
   par(new=TRUE)
   plot(out_ls$PT_chain[[chain]][low:up,2],out_ls$PT_chain[[chain]][low:up,1],col='red',xlab=expression(alpha),ylab=expression(beta),main=expression(paste('A. Joint prior and posterior of ', alpha,' and ', beta,sep='')),xaxt='n',yaxt='n',xlim=range(f1$x),ylim=range(f1$y),cex.main=4,cex=3,cex.lab=3)
   par(new=TRUE, oma=c(6,3,1,2))
   ## create a layout to plot the subplot in the right bottom corner
   layout(matrix(1:4,2))
   plot(out_ls$PT_chain[[chain]][low:up,2],out_ls$PT_chain[[chain]][low:up,1],col='red',xlab=expression(alpha),ylab=expression(beta),main='Posterior (zoomed in)',cex.main=4,cex=3,cex.lab=2.5,cex.axis=2.5)
   dev.off()
   #layout(matrix(1:2,ncol=2))
   #pdf('PriorI0.pdf',width = 9, height = 8)
   #png('PriorI0.png',width=650,height=550)
   setEPS()
   postscript("FIG8.eps",horizontal=FALSE, paper="special",height=10,width=19, colormodel = "cmyk", 
              family = "Helvetica")
    
    
   par(mfrow=c(1,1), mar=c(6,6.4, 5.1, 5.1), xpd=TRUE)
   N=261
   pl=rbinom(100000,N,5/N)
   h=hist(pl,nclass=50,plot=FALSE)
   h$counts=h$counts/sum(h$counts)
   
   h1=hist(out_ls$PT_chain[[chain]][low:up,4],breaks=seq(range(pl)[1],range(pl)[2],by=0.5), plot=FALSE)
   h1$counts=h1$counts/sum(h1$counts)
   
   plot(h,main='B. Prior and posterior distribution of I(0)',col='lightgrey',xlab='',ylab='',cex.axis=3,cex.main=3.5,xlim= range(pl),ylim=range(h1$counts))
   par(new=T)
   plot(h1,xlab='I(0)',xaxt='n',col='red',yaxt='n',main='',cex.lab=3,xlim= range(pl),ylim=range(h1$counts))
   box(lwd=2, lty=1,col = 'grey')
   
   
   legend('topright', 
          c("prior","posterior"), inset=c(-0.0002,0),
          lty=c(1,1), cex=3,border='lightgrey',
          lwd=c(2,2),col=c("lightgrey","red")) 
   dev.off()
   
   
   
   
   
   
 
   
  