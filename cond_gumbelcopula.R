# this function helps in finding conditional probability of streamflow conditioned 
# upon soil mositure or rainfall
# input : copula parameter, bi_obs_cdf, marginal_cdfs
library(SciViews)
con_prob=function(alpha,bi_cdf,u, v){
   
  temp=matrix(data=NA, nrow=length(bi_cdf),ncol=1)
  for(i in 1:length(bi_cdf)){
    temp[i,]=bi_cdf[i]*(((-ln(u[i]))^(1/alpha))+((-ln(v[i]))^(1/alpha)))^(-1+(1/alpha))* ((-ln(u[i]))^(alpha-1))/u[i]
  }
  
  #kendall_t <- function(alpha, t){
  #  par = alpha
    
   # kt <- rep(NA,length(t))
    #kt <- t - t * log(t)/(par)
    #kt[t==1] <- 1
    #kt[t==0] <- 0
    #return(kt)  
  #}
  
  #k_t=kendall_t(alpha,temp)
  
  rt_cond=1/(1-temp)
  
  d=list(temp,rt_cond)
  
  names(d)=c("Ccdf","C_rt")
  
  return(d)
}

