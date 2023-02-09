##### fitting bicopula function

#family=14 ## survival gumbel copula
bifit=function(uni_data){
  library(VineCopula)
  library(evd)
  library(stats)
  library(evd)
  library(copula)
  source("ecdf_kendall.r")
  u1=uni_data[,1]
  u2=uni_data[,2]
  tau=cor(u1,u2, method="kendall")
  cop=c(3,4,5,6,7,8,9,10)
  #rho <- BiCopTau2Par(family, tau ) "itau procedure"
  #model=BiCop(family, par = rho, par2 = 0)## alternatevely select best among available set
  ## maximum likelihood apporach to select the copula family
  modl=BiCopSelect(u1,u2, familyset = c(3:10), selectioncrit = "AIC", indeptest = FALSE, level = 0.05, weights = NA, rotations = TRUE, se = TRUE, presel = TRUE, method = "mle" )
  ## goodness of model
  
  F_bi=BiCopCDF(u1,u2,modl)
  n <- nrow(uni_data)
  d <- ncol(uni_data) 
  
  ecdf <- vapply(seq_len(n), function(i) sum(colSums(t(uni_data) < uni_data[i, ]) == d)/(n + 1), NA_real_) # empirical distribution
  
  obs_bi_cdf=data.frame(uni_data,ecdf,F_bi)
  colnames(obs_bi_cdf)=c("u1","u2","ecdf","F_bi")
  
  ## hypothesis testing
  
  ## ks test ## AD test
  # AD test
  ad=ad.test(F_bi,ecdf, method="exact", dist=FALSE, Nsim=1000)
  # KS Test
  ks=ks.test(F_bi,ecdf)
  
  test=as.data.frame(c("AD test","KS Test"))
  test["t"]=as.data.frame(c("AD test","KS test"))
  test["p-value"]=c(ad[["ad"]][[5]],ks[["p.value"]])
  test=subset(test, select = c("t","p-value"))
  
  ## graphical plot
  plot(ecdf,F_bi, xlab="Empirical", ylab="Model")
  abline(0,1)
  
  ## independene test
  dep=BiCopIndTest(u1,u2)
  # graphical aid
  t=BiCopKPlot(u1,u2,PLOT=FALSE)
  plot(t$W.in,t$Hi.sort)
  abline(0,1)
  BiCopKPlot(u1,u2,PLOT=TRUE)
  BiCopChiPlot(u1,u2,PLOT=TRUE,xlim = c(-1,1), ylim = c(-1,1) )
  
  
  ### model building is done
  
  
  ## simulation from model
  sim_data=BiCopSim(25000,modl)
  F_cdf_sim=BiCopCDF(sim_data[,1],sim_data[,2],modl)
  t=empCopula(sim_data)
  tt=F.n(sim_data,sim_data)
  ##empirical cdf values
  n_sim <- nrow(sim_data)
  d_sim <- ncol(sim_data) 
  sim=cbind(sim_data[,1],sim_data[,2])
  ecdf_sim <- vapply(seq_len(n_sim), function(i) sum(colSums(t(sim) < sim[i, ]) == d_sim)/(n_sim + 1), NA_real_) # empirical distribution

  ### 
  sim_bi_cdf=data.frame(sim,ecdf_sim,F_cdf_sim)
  colnames(sim_bi_cdf)=c("u1","u2","ecdf_sim","F_cdf_sim")
  
  model_bi=list(tau,modl,obs_bi_cdf,sim_bi_cdf,test,dep)
  names(model_bi)=c("kendall's tau","model","observed cdf","simulated cdf","hypothesis test","dependence test")
  
  return(model_bi)
  
}


 
 
  #fig=plot_ly(x=modeled$X1,y=modeled$X2, z=modeled$F_cdf_sim, type="contour",contours=list(coloring="lines",showlabels=TRUE))%>%
  #  add_trace(x=empirical$X1,y=empirical$X2,z=empirical$ecdf_sim,type="contour", contours=list(coloring="lines",showlabels=TRUE))%>%
   # add_trace(x=as.vector(u1),y=as.vector(u2), type="scatter", mode="markers", color=I("black"))
  #fig
  

  
  
  
  