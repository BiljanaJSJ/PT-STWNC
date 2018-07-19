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


mlik=get(load('marglikelihood.RData'))

mllikB=mlik[[1]][1,]


mllikNB=mlik[[1]][2,]


logBij=matrix(NA,5,5)
for (i in (1:5)){
  for (j in (1:5)){
    logBij[i,j]=mllikB[i]-mllikB[j]
  }
}
rownames(logBij)=colnames(logBij)=1:5
exp(logBij)


logBijNB=matrix(NA,5,5)
for (i in (1:5)){
  for (j in (1:5)){
    logBijNB[i,j]=mllikNB[i]-mllikNB[j]
  }
}
rownames(logBijNB)=colnames(logBijNB)=1:5
exp(logBijNB)




# > logBij
# 1            2            3            4            5
# 1  0.00000 -13.75466104 -13.82066705 -14.13993976 -14.16593029
# 2 13.75466   0.00000000  -0.06600601  -0.38527872  -0.41126926
# 3 13.82067   0.06600601   0.00000000  -0.31927271  -0.34526324
# 4 14.13994   0.38527872   0.31927271   0.00000000  -0.02599053
# 5 14.16593   0.41126926   0.34526324   0.02599053   0.00000000

# > logBijNB
# 1            2            3            4            5
# 1  0.00000 -13.74392240 -13.79634901 -14.14119816 -14.17512827
# 2 13.74392   0.00000000  -0.05242661  -0.39727576  -0.43120587
# 3 13.79635   0.05242661   0.00000000  -0.34484916  -0.37877926
# 4 14.14120   0.39727576   0.34484916   0.00000000  -0.03393011
# 5 14.17513   0.43120587   0.37877926   0.03393011   0.00000000


newPT_STWNC=c(-240.36,-226.96,-232.68,-224.51,-222.91)

logBijPT_STWNC=matrix(NA,5,5)
for (i in (1:5)){
  for (j in (1:5)){
    logBijPT_STWNC[i,j]=newPT_STWNC[i]-newPT_STWNC[j]
  }
}
rownames(logBijPT_STWNC)=colnames(logBijPT_STWNC)=1:5
exp(logBijPT_STWNC)
