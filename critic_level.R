## this function helps in getting critical levels at 2,5,10,25,50,100,200 years return levels
pars_critic=function(model){
setwd("F:/Research Work/Synchronization work/Multivariate IDF/R Codes/spcopula-master/R")
source("returnperiods.r")

library(copula)
dim=2
name=model[["familyname"]]
par=model[["par"]]
if(model[["par2"]]!=0){
  par2=model[["par2"]]
}


mod=copula(dim,par,name)
obj=archmCopula(name,par,dim)

cl_rt=criticalLevel(getKendallDistr(obj),KRP=c(2,5,10,25,50,100,200))

u=seq(from=0.0001,to=0.999,length.out=500)

k_RT=list()

for( i in 1:length(cl_rt)){
  kk=criticalPair(obj,cl_rt[i],u,1)
  k_RT=append(k_RT,list(kk))
}
names(k_RT)=c("2","5","10","25","50","100","200")
}

t=data.frame(u,k_RT[1])
tt=rep(2,times=500)
t["rt"]=tt
fig=plot_ly(x=t$u,y=t$X200,z=t$rt,type = "contour")
fig

