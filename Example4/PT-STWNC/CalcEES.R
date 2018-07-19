library(coda)
library(mcmcse)
out_ls=get(load('/zfs/users/b.stojkova/b.stojkova/BNN/ParallelInputRestartNew/ST455000.RData'))
out_ls_PT=get(load('/zfs/users/b.stojkova/b.stojkova/BNN/BNN_PT_fixedInput/PT_ThermodynInt_150000.RData'))

up=455000
low=up/2


ess_multi_coda=effectiveSize(as.mcmc.list(lapply(as.data.frame(out_ls$PT_chain[[2]][low:up,]), mcmc)))
109476.9



up=150000
low=up/2


ess_multi_coda_PT=effectiveSize(as.mcmc.list(lapply(as.data.frame(out_ls_PT$mu[[15]][low:up,]), mcmc)))
2619.904

