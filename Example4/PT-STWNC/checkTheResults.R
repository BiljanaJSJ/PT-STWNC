source('BNN_ST.R')
source('BNN_ST_functions.R')
data         = read.table("gas-furnace.csv",sep = ",",header =  T)
y=data$CO2[1:206]
p=4
M=5
up=526000
low=1

library(coda)
out_ls=get(load('/zfs/users/b.stojkova/b.stojkova/BNN/ParallelInputRestartNew/ST526000.RData'))

library(doSNOW)
cl = makeCluster(getOption("cl.cores", 4),type="SOCK")
gamma_ij_hat=list()

gamma_ij_hat=lapply(1:dim(out_ls$PT_chain[[2]][low:up,])[1], function(l) {matrix(out_ls$PT_chain[[2]][l,((p+M+2):(M*(p+1)+M+p+1))],ncol=M,byrow=T)})
x=generateInput(data$CO2[1:206],p=4)
                                                         
beta_j = out_ls$PT_chain[[2]][low:up,6:10]
lambda  = out_ls$PT_chain[[2]][low:up,1:5]

clusterExport(cl, list=c("x",'beta_j','M','gamma_ij_hat','psi','lambda','generateInput'), envir = .GlobalEnv)

#this takes tooo loooong! Try to Optimize the code, parallel
predictions <- parLapply(cl,1:dim(beta_j)[1], function(e)  {sapply(1:nrow(x), function(l) {sum(t(x[l,])*lambda[e,],na.rm=T)}) + apply(simplify2array(lapply(1:M, function(l) {sapply(1:nrow(x), function(k) { beta_j[e,l]*psi( sum(t(x[k,])*gamma_ij_hat[[ e]][,l],na.rm=T) ) } )  })), 1, sum, na.rm=T)})

save(predictions,file=paste('PredST',up,'.RData',sep=''))



y_hat=rep(NA,296)
y_hat[1:206]=apply(simplify2array(predictions)[,(up/2):up], 1,mean, na.rm=T)

#plot(data$CO2[1:206])
#lines(y_hat[1:206], col='red')

(trainingError=mean((y_hat[1:206]-data$CO2[1:206])^2,na.rm=T))


nPred=296
MSPE=rep(NA,3)
pred_oneAhead=list()
y_hat=c()


#cl = makeCluster(getOption("cl.cores", 4),type="SOCK")

beta_j = out_ls$PT_chain[[2]][low:up,6:10]
lambda  = out_ls$PT_chain[[2]][low:up,1:5]
gamma_ij_hat=lapply(1:dim(out_ls$PT_chain[[2]][low:up,])[1], function(l) {matrix(out_ls$PT_chain[[2]][l,((p+M+2):(M*(p+1)+M+p+1))],ncol=M,byrow=T)})
clusterExport(cl, list=c("beta_j","M","gamma_ij_hat","psi","lambda","generateInput"), envir = .GlobalEnv)  
#one-step ahead from multiple steps ahead
l=1;MSPE=c()
nPred=296
y_hat=c()
pred_oneAhead=list()
y_hats=matrix(NA,nrow=up,ncol=296)

for (i in (207:(nPred) )){
  print(i)
   #pred_oneAhead[[i]]=c()
   # y_hat[[l]][i]=NULL
  if (i==207){
    x=generateInput(data$CO2[1:i],p=4)[i,]
    clusterExport(cl, list=c("x","beta_j","M","gamma_ij_hat","psi","lambda","generateInput"), envir = .GlobalEnv)  
    pred_oneAhead[[i]]=parLapply(cl,1:dim(beta_j)[1], function(e)  {sum(t(x)*lambda[e,],na.rm=T) + sum(unlist(lapply(1:M, function(l) {beta_j[e,l]*psi( sum(t(x)*gamma_ij_hat[[ e]][,l],na.rm=T) ) } )), na.rm=T)})
 
  }else{ #208+1,209+1,..,296-1+1
    #here we will need a matrix object for y[samples,datapoints] 
    clusterExport(cl, list=c("y_hats","i","y_hat","data","beta_j","M","gamma_ij_hat","psi","lambda","generateInput"), envir = .GlobalEnv) 
    xmat=parLapply(cl,1:up, function(x) { generateInput(c(data$CO2[1:206],y_hats[x,207:i],1),p=4)[i,] } )

    
    clusterExport(cl, list=c("xmat","y_hats","i","y_hat","data","beta_j","M","gamma_ij_hat","psi","lambda","generateInput"), envir = .GlobalEnv) 
    pred_oneAhead[[i]]=unlist( parLapply(cl, 1:dim(beta_j)[1], function(e)  {sum(t(xmat[[e]])*lambda[e,],na.rm=T) + sum(unlist(lapply(1:M, function(l) {beta_j[e,l]*psi( sum(t(xmat[[e]])*gamma_ij_hat[[ e]][,l],na.rm=T) ) } )), na.rm=T)}))
  }
  y_hats[,i]=unlist(pred_oneAhead[[i]])
  y_hat[i]=mean(unlist(pred_oneAhead[[i]])[((up/2):up)],na.rm=T)  
}

for (l in (1:3)){
MSPE[l]=mean((data$CO2[(207+(l-1)):(nPred)]-y_hat[(207+(l-1)):(nPred)])^2,na.rm=T)
}

savelst=list(trainingError=trainingError,y_hats=y_hats,y_hat=y_hat,MSPE=MSPE, pred_oneAhead= pred_oneAhead)
save(savelst,file='Predictions1.RData')


