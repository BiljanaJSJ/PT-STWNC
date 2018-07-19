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

listfolders=list.files('/Users/bjonoska/Documents/BiasGalaxyNewPriors_7_4_2016/STWNC_mllik', full.names = TRUE)
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


save(table,file='Table_paper.RData')



library(xtable)
xtable(table)


% latex table generated in R 3.2.2 by xtable 1.8-2 package
% Mon Dec 12 08:55:00 2016
\begin{table}[ht]
\centering
\begin{tabular}{rrr}
  \hline
 & marignal likelihood & SD \\
  \hline
model2e & -247.51 & 0.76 \\
  model3 & -234.93 & 0.65 \\
  model3e & -240.28 & 0.67 \\
  model4 & -231.04 & 0.89 \\
  model5 & -227.46 & 0.76 \\
   \hline
\end{tabular}
\end{table}
