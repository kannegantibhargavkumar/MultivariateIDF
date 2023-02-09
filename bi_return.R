bi_return=function(ams){
  # input: Annual maxima series of two variables
  ## required packages
  library(evd)
  library(copula)
  source("ecdf_kendall.r")
  source("hypothesis_test.r")
  # univariate distribution
  ## soil moisture
  x1=ams[,1]
  x1_dist=fgev(x1)
  F_x1=pgev(x1, loc=as.numeric(x1_dist$param[1]), scale = as.numeric(x1_dist$param[2]), shape = as.numeric(x1_dist$param[3]))
  rt_x1=1/(1-F_x1)
  ## Streamflow
  x2=ams[,2]
  x2_dist=fgev(x2)
  F_x2=pgev(x2, loc=as.numeric(x2_dist$param[1]), scale = as.numeric(x2_dist$param[2]), shape = as.numeric(x2_dist$param[3]))
  rt_x2=1/(1-F_x2)
  
  
  modl=BiCopSelect(F_x1,F_x1, familyset = NA, selectioncrit = "AIC", indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE, se = TRUE, presel = TRUE, method = "mle" )
  
  
  
  # bivariate fitting
  #U=cbind(F_x1,F_x2)
  # kendall correlation
  #tau=cor(x1,x2, method="kendall")
  # copula parameter
  #alpha=iTau(gumbelCopula(), tau=tau)
  # dimension 2
  #d=2
  # gumbel model
  #gc=gumbelCopula(alpha, dim=d)
  # bivariate distribution function
  #F_x1x2=pCopula(U, gc)
  F_x1x2=BiCopCDF(F_x1,F_x2,obj=modl)
  # 
  ## kendall distribution/measure, taken from VineCopula:::obs.stat
  
  kendall.Gumbel <- function(copula, t){
    #par = copula@parameters
    par=copula$
    
    kt <- rep(NA,length(t))
    kt <- t - t * log(t)/(par)
    kt[t==1] <- 1
    kt[t==0] <- 0
    return(kt)  
  }
  
  ken_cdf=kendall.Gumbel(gc,F_x1x2)
  
  #cdf
  p_and=(1-(F_x1+F_x2-F_x1x2))
  
  p_or=(1-F_x1x2)
  
  p_cond1= (1-(F_x1+F_x2-F_x1x2))/(1-F_x1)
  
  p_cond2=1-(F_x1x2/F_x1)
  
  p_k=1-ken_cdf
  
  F_all=data.frame(p_and,p_or,p_k, p_cond1, p_cond2 )
  
  # return periods
  rt_ken=1/(1-ken_cdf)
  
  rt_bi_and=1/(1-(F_x1+F_x2-F_x1x2))
  
  rt_bi_or=1/(1-F_x1x2)
  
  rt_cond1=(1-(F_x1+F_x2-F_x1x2))/(1-F_x1)
  
  rt_cond2=1-(F_x1x2/F_x1)
  
  
  
  
  rt=data.frame(rt_x1,rt_x2,rt_bi_and, rt_bi_or, rt_cond1,rt_cond2)
  
  #### Goodness fit of copula function
  ken_emp_fun <- genEmpKenFun(gc,cbind(F_x1,F_x2))
  f_ecdf=ken_emp_fun(F_x1x2)
 
  bi_cdf=data.frame(f_ecdf,ken_cdf)
  
  # AD test
  ad=ad.test(as.matrix(bi_cdf[,1]),as.matrix(bi_cdf[,2]), method="exact", dist=FALSE, Nsim=1000)
  # KS Test
  ks=ks.test(as.matrix(bi_cdf[,1]),as.matrix(bi_cdf[,2]))
  
  test=as.data.frame(c("AD test","KS Test"))
  test["t"]=as.data.frame(c("AD test","KS test"))
  test["p-value"]=c(ad[["ad"]][[5]],ks[["p.value"]])
  test=subset(test, select = c("t","p-value"))
  
  
  
  
  
  
  
  bi_model=list(gc,F_x1x2,ken_cdf,F_all, bi_cdf,rt_ken ,rt,test)
  
  names(bi_model)=c("model","obs_CDF","ken_cdf","cdf","kendall" ,"rt_ken","Rp_bi_obs","Hypothesis Test")
  return(bi_model)
}



