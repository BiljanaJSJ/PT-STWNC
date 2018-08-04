
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

library(coda)

pred=get(load('Predictions1_455.RData'))
pred_PT=get(load('/zfs/users/b.stojkova/b.stojkova/BNN/BNN_predictions/PT/Predictions1.RData'))
up1=150000

nPred=296;up=455000
mode_y_hat=c()
sd_y_hat=matrix(NA, nrow=nPred,ncol=2)
sd_y_hat_PT=matrix(NA, nrow=nPred,ncol=2)


data         = read.table("gas-furnace.csv",sep = ",",header =  T)
y=data$CO2


for (i in (207:(nPred) )){
#sd_y_hat[i,]=quantile(pred$y_hats[,i][((up/2):up)],c(0.05, 0.975), na.rm=T) 
sd_y_hat_PT[i,]= HPDinterval(mcmc(pred_PT$y_hats[((up1/2):up1),i]),0.95, na.rm=T) 

}


for (i in (207:(nPred) )){
#sd_y_hat[i,]=quantile(pred$y_hats[,i][((up/2):up)],c(0.05, 0.975), na.rm=T) 
sd_y_hat[i,]= HPDinterval(mcmc(pred$y_hats[((up/2):up),i]),0.95, na.rm=T) 

}



trainpred=get(load('PredST455000.RData'))
trainpred_PT=get(load('/zfs/users/b.stojkova/b.stojkova/BNN/BNN_predictions/PT/Pred_PT_150000.RData'))

y_hat_train=rep(NA,206)
y_hat_train[1:206]=apply(simplify2array(trainpred)[,(up/2):up], 1,mean, na.rm=T)

sd_y_hat_train=sd_y_hat=matrix(NA, nrow=dim(simplify2array(trainpred))[1],ncol=2)

mat_trainpred=t(simplify2array(trainpred))
for (i in (1:206 )){
##sd_y_hat_train[i,]=quantile(mat_trainpred[((up/2):up),i],c(0.05, 0.975), na.rm=T) 
sd_y_hat_train[i,]=HPDinterval(mcmc(mat_trainpred[((up/2):up),i]),0.95, na.rm=T) 
}


y_hat_train_PT=rep(NA,206)
y_hat_train_PT[1:206]=apply(simplify2array(trainpred_PT)[,(up1/2):up1], 1,mean, na.rm=T)

sd_y_hat_train_PT=sd_y_hat=matrix(NA, nrow=dim(simplify2array(trainpred_PT))[1],ncol=2)

mat_trainpred_PT=t(simplify2array(trainpred_PT))
for (i in (1:206 )){
##sd_y_hat_train[i,]=quantile(mat_trainpred_PT[((up1/2):up1),i],c(0.05, 0.975), na.rm=T) 
sd_y_hat_train_PT[i,]=HPDinterval(mcmc(mat_trainpred_PT[((up1/2):up1),i]),0.95, na.rm=T) 
}


#posterior predictions and highest posterior density credible interval 

x=207:209;x1=1:206
setEPS()
postscript("FIG10.eps",horizontal=FALSE, paper="special",height=22,width=22, colormodel = "cmyk", 
           family = "Helvetica")
par(mfrow=c(2,2),mar=c(7,12,12,7))

plot(x1,y_hat_train[x1], ylim=c(45,62),cex=0.5 ,xaxt = "n",col='blue', cex.axis=3.5, cex.main=4, cex.lab=3.5, xlab="time", ylab='Posterior prediction',main='A. PT-STWNC: posterior predictions \nand 95% HPD credible intervals, \n training data set', mgp= c(4.5, 2, 1))
axis(1, at=round(seq(x1[1],x1[206],length=5),0), labels=round(seq(x1[1],x1[206],length=5),0),cex.axis=3.5)
polygon(c(rev(x1), x1), c(rev(sd_y_hat_train[,1]), sd_y_hat_train[,2]),border='grey48', col = 'grey80', pch=19, lty=2)
points(x1,y[x1],col='red',cex=3.5, pch=19)
lines(x1,y_hat_train[x1],col='blue',cex=2.5)
points(x1,y_hat_train[x1],col='blue',cex=3.5, pch=19)

plot(x,pred$y_hat[x], ylim=c(56,62),cex=0.5 ,xaxt = "n",col='blue', cex.axis=3.5, cex.main=4, cex.lab=3.5, xlab="time", ylab='Posterior prediction',main='B. PT-STWNC: posterior predictions \nand 95% HPD credible intervals, \none, two and three steps ahead', mgp= c(4.5, 2, 1))
axis(1, at=x, labels=x,cex.axis=3.5)
polygon(c(rev(x), x), c(rev(c(sd_y_hat[x,1])), c(sd_y_hat[x,2])),border='grey48', col = 'grey80', pch=19, lty=2)
points(x,y[x],col='red',cex=3.5, pch=19)
lines(x,pred$y_hat[x],col='blue',cex=2.5)
points(x,pred$y_hat[x],col='blue',cex=3.5, pch=19)

plot(x1,y_hat_train_PT[x1], ylim=c(45,62),cex=0.5 ,xaxt = "n",col='blue', cex.axis=3.5, cex.main=4, cex.lab=3.5, xlab="time", ylab='Posterior prediction',main='C. PT: posterior predictions \nand 95% HPD credible intervals, \n trainig data set', mgp= c(4.5, 2, 1))
axis(1, at=round(seq(x1[1],x1[206],length=5),0), labels=round(seq(x1[1],x1[206],length=5),0),cex.axis=3.5)
polygon(c(rev(x1), x1), c(rev(sd_y_hat_train_PT[,1]), sd_y_hat_train_PT[,2]),border='grey48', col = 'grey80', pch=19, lty=2)
points(x1,y[x1],col='red',cex=3.5, pch=19)
lines(x1,y_hat_train_PT[x1],col='blue',cex=2.5)
points(x1,y_hat_train_PT[x1],col='blue',cex=3.5, pch=19)

plot(x,pred_PT$y_hat[x], ylim=c(56,62),cex=0.5 ,xaxt = "n",col='blue', cex.axis=3.5, cex.main=4, cex.lab=3.5, xlab="time", ylab='Posterior prediction',main='D. PT: posterior predictions \nand 95% HPD credible intervals, \none, two and three steps ahead ', mgp= c(4.5, 2, 1))
axis(1, at=x, labels=x,cex.axis=3.5)
polygon(c(rev(x), x), c(rev(c(sd_y_hat_PT[x,1])), c(sd_y_hat[x,2])),border='grey48', col = 'grey80', pch=19, lty=2)
points(x,y[x],col='red',cex=3.5, pch=19)
lines(x,pred_PT$y_hat[x],col='blue',cex=2.5)
points(x,pred_PT$y_hat[x],col='blue',cex=3.5, pch=19)

dev.off()








out_ls=get(load('/zfs/users/b.stojkova/b.stojkova/BNN/ParallelInputRestartNew/ST455000.RData'))
source('pairs2_function.R')
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


#l=c(2,10,35)
#setEPS()
#pairs(samples$PT_chain[[2]][((up/2):up),l])

#setEPS()
#postscript("FIG9.eps",horizontal=FALSE, paper="special",height=12,width=19, colormodel = "cmyk", 
#            family = "Helvetica")

#par(mfrow=c(2,2))
#hist(out_ls$PT_chain[[1]][((up/2):up),l[1]],main=expression(lambda),cex.axis=2.5,cex.main=3.5,nclass=150,col='cyan', probability = T, xlab   ='', ylab   ='')
#hist(out_ls$PT_chain[[1]][((up/2):up),l[2]],main=expression(beta),cex.axis=2.5,cex.main=3.5,nclass=150,col='cyan', probability = T, xlab   ='', ylab   ='')
#hist(out_ls$PT_chain[[1]][((up/2):up),l[3]],main=expression(gamma),cex.axis=2.5,cex.main=3.5,nclass=150,col='cyan', probability = T, xlab   ='', ylab   ='')
#hist(out_ls$tau[((up/2):up)],main=expression(tau),cex.axis=2.5,cex.main=3.5,nclass=150,col='cyan', probability = T, xlab   ='', ylab   ='')
#dev.off()

l=c(4,10,32)
setEPS()
postscript("FIG9.eps",horizontal=FALSE, paper="special",height=15,width=19, colormodel = "cmyk", 
           family = "Helvetica")
pairs2(out_ls$PT_chain[[2]][((up/2):up),l],labels=c(expression(lambda),expression(beta),expression(gamma)), bg="cyan",
      diag.panel=panel.hist,cex.labels = 4, font.labels=4,cex.axis=4,upper.panel=NULL, oma=c(2,3,2,2),gap=8,label.pos =0.9,mgp = c(5, 3, 1))
dev.off()






#png("predicProbPlots.png", width = 10,  height = 10,  units = "in",  res = 1200, pointsize = 4)
#color=c(rep('black', 207), rep('cyan',3))
#plot(c(y_hat_train[1:206],pred$y_hat[207:210]), ylim=c(45,60),cex=1.5, col=color)
#polygon(c(rev(1:210), 1:210), c(rev(c(sd_y_hat_train[,1],sd_y_hat[207:210,1])), c(sd_y_hat_train[,2],sd_y_hat[207:210,2])), col = 'grey80', border = NA)
#points(y[1:210],col='red',cex=1.5, pch=19)
#points(c(y_hat_train[1:206],pred$y_hat[207:210]),col=color,cex=2, pch=19)
#lines(c(sd_y_hat_train[,1],sd_y_hat[207:210,1]),col=color)
#lines(c(sd_y_hat_train[,2],sd_y_hat[207:210,2]),col=color)
#dev.off()


#x=207:209;x1=1:206
#setEPS()
#postscript("FIG10.eps",horizontal=FALSE, paper="special",height=19,width=16, colormodel = "cmyk", 
#           family = "Helvetica")
#par(mfrow=c(2,1),mar=c(7,7,7,7))

#plot(c(x1,x),c(y_hat_train[x1],pred$y_hat[x]), ylim=c(45,62),cex=0.5 ,xaxt = "n",col='blue', cex.axis=3.5, cex.main=4, cex.lab=3.5, xlab="time", ylab='Posterior prediction',main='A. PT-STWNC: posterior predictions and \n95% HPD credible intervals', mgp= c(4.5, 2, 1))
#axis(1, at=c(seq(1,206, by=7),seq(206,209, by=2)), labels=c(seq(1,206, by=7),seq(206,209, by=2)),cex.axis=3.5)
#polygon(c(rev(c(x1,x)), c(x1,x)), c(rev(c(sd_y_hat_train[,1],sd_y_hat[x,1])), c(sd_y_hat_train[,2],sd_y_hat[x,2])),border='grey48', col = 'grey80', pch=19, lty=2)
#points(c(x1,x),y[c(x1,x)],col='red',cex=3.5, pch=19)
#lines(x,pred$y_hat[x],col='blue',cex=2.5)
#points(x,pred$y_hat[x],col='blue',cex=3.5, pch=19)


#plot(c(x1,x),c(y_hat_train_PT[x1],pred_PT$y_hat[x]), ylim=c(45,62),cex=0.5 ,xaxt = "n",col='blue', cex.axis=3.5, cex.main=4, cex.lab=3.5, xlab="time", ylab='Posterior prediction',main='B. PT: posterior predictions and \n95% HPD credible intervals', mgp= c(4.5, 2, 1))
#axis(1, at=c(seq(1,206, by=7),seq(206,209, by=2)), labels=c(seq(1,206, by=7),seq(206,209, by=2)),cex.axis=3.5)
#polygon(c(rev(c(x1,x)), c(x1,x)), c(rev(c(sd_y_hat_train_PT[,1],sd_y_hat_PT[x,1])), c(sd_y_hat_train_PT[,2],sd_y_hat_PT[x,2])),border='grey48', col = 'grey80', pch=19, lty=2)
#points(c(x1,x),y[c(x1,x)],col='red',cex=3.5, pch=19)
#lines(x,pred_PT$y_hat[x],col='blue',cex=2.5)
#points(x,pred_PT$y_hat[x],col='blue',cex=3.5, pch=19)

#dev.off()






