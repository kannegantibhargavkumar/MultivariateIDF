### hypothesis testing function for paired samples
library(kSamples)
library(snpar)
library(ggplot2)

ht=function(data){
  # paired t test
  t=t.test(as.matrix(data[,1]),as.matrix(data[,2]), paired=TRUE, alterantive= "two-sided")
  # wilcoxin-sign test for paired sample
  w=wilcox.test(as.matrix(data[,1]),as.matrix(data[,2]), paired=TRUE, alterantive= "two-sided")
  # AD test
  ad=ad.test(as.matrix(data[,1]),as.matrix(data[,2]), method="exact", dist=FALSE, Nsim=1000)
  # KS Test
  ks=ks.test(as.matrix(data[,1]),as.matrix(data[,2]))
  # quantile test
  t_0.75=quant.test(as.matrix(data[,1]),as.matrix(data[,2]), paired=TRUE, p=0.75)
  
  t_0.5=quant.test(as.matrix(data[,1]),as.matrix(data[,2]),  paired=TRUE, p=0.5)
  
  t_0.25=quant.test(as.matrix(data[,1]),as.matrix(data[,2]), paired=TRUE, p=0.25)
  # summary of tests
  options(scipen=999)
  test=as.data.frame(c("T-test","Wilcoxin-Sign Rank Test","AD test","KS Test","Quantile Test 0.75","Quantile Test 0.5","Quantile Test 0.25"))
  test["t"]=as.data.frame(c("T-test","Wilcoxin-Sign Rank Test","AD test","KS test","Qunatile Test 0.75","Qunatile Test 0.5","Quantile Test 0.25"))
  test["p-value"]=c(t$p.value,w$p.value,ad[["ad"]][[5]],ks,t_0.75$p.value,t_0.5$p.value,t_0.25$p.value)
  test=subset(test, select = c("t","p-value"))
  
  return(test)
}

