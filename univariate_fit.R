## This function contain fitting of GEV distribution to maxima series
library(evd)
library(kSamples)
library(snpar)
library(ggplot2)
marginal_fit=function(data){
  fit=fgev(as.matrix(data), std.err = TRUE)
  par=fit$param
  F_x=pgev(data, loc=as.numeric(par[1]), scale=as.numeric(par[2]), shape=as.numeric(par[3]))
  rt=1/(1-F_x)
  fun.cdef=ecdf(data)
  my_ecdf=fun.cdef(data)
  data_ecdf=data.frame(data, my_ecdf)
  data_test=cbind(F_x,data_ecdf$my_ecdf)
  
  
  # paired t test
  t=t.test(as.matrix(data_test[,1]),as.matrix(data_test[,2]), paired=TRUE, alterantive= "two-sided")
  # wilcoxin-sign test for paired sample
  w=wilcox.test(as.matrix(data_test[,1]),as.matrix(data_test[,2]), paired=TRUE, alterantive= "two-sided")
  # AD test
  ad=ad.test(as.matrix(data_test[,1]),as.matrix(data_test[,2]), method="exact", dist=FALSE, Nsim=1000)
  # KS Test
  ks=ks.test(as.matrix(data_test[,1]),as.matrix(data_test[,2]))
  # quantile test
  t_0.75=quant.test(as.matrix(data_test[,1]),as.matrix(data_test[,2]), paired=TRUE, p=0.75)
  
  t_0.5=quant.test(as.matrix(data_test[,1]),as.matrix(data_test[,2]),  paired=TRUE, p=0.5)
  
  t_0.25=quant.test(as.matrix(data_test[,1]),as.matrix(data_test[,2]), paired=TRUE, p=0.25)
  # summary of tests
  options(scipen=999)
  test=as.data.frame(c("T-test","Wilcoxin-Sign Rank Test","AD test","KS Test","Quantile Test 0.75","Quantile Test 0.5","Quantile Test 0.25"))
  test["t"]=as.data.frame(c("T-test","Wilcoxin-Sign Rank Test","AD test","KS test","Qunatile Test 0.75","Qunatile Test 0.5","Quantile Test 0.25"))
  test["p-value"]=c(t$p.value,w$p.value,ad[["ad"]][[5]],ks[["p.value"]],t_0.75$p.value,t_0.5$p.value,t_0.25$p.value)
  test=subset(test, select = c("t","p-value"))
  
  d=list(par,F_x,rt, test)
  names(d)=c("parameters","CDF", "RP", "hypothesis")
  
  
  return(d)
}
