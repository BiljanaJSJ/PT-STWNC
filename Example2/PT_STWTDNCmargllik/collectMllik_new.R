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


#collect the marginal likelihood

#setwd('/stat/stat2/bjonoska/ST/newST/Gallaxy/NewGalaxyCode9/multipleruns')

low = 17500
up  = 35000

listfolders=list.files('/zfs/users/b.stojkova/b.stojkova/PT-STWNC/Example12_March12/PT_STWTDNCmargllik', full.names = TRUE)
n=listfolders[4:8]
mllik=matrix(NA,up-low+1,20)
mllik_se=mllik_se1=rep(NA,20)


table=matrix(NA,length(n),2)
for (i in (1:length(n))){
out_ls=get(load(paste(n[i],'/savelist_allChains.RData',sep='' )))
table[i,1]=mean(out_ls[[1]])
table[i,2]=sd(out_ls[[1]])
}
colnames(table)=c('marignal likelihood','SD')
rownames(table)=c('model2e','model3','model3e','model4','model5')



#Bayes Factors
BF=rep(NA,9)
for (i  in (2:5)){
BF[i]=table[i,1]-table[1,1]
}

for (i  in (3:5)){
BF[5+i-2]=table[i,1]-table[2,1]
}
for (i  in (4:5)){
BF[8+i-3]=table[i,1]-table[3,1]
}
i=5
BF[10+i-4]=table[i,1]-table[4,1]

#	$\log{BF21}$	$\log{BF31}$ $\log{BF41}$ $\log{BF51}$ $\log{BF32}$ $\log{BF42}$ & 2.45 & 0.39 & 0.38 \\ 
#		$\log{BF52}$ $\log{BF43}$ $\log{BF53}$ 	$\log{BF54}$

retls=list(table=table,BF=BF)
save(retls,file='Table_paper.RData')


library(xtable)
xtable(table)
print(BF)

