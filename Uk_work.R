## work on UK basin
library("readxl")
library("fitur")
library("evd")
library(readr)
library(copula)
library(snpar)
library(ggplot2)
library(Ternary)
library(plotly)
library(ggtern)
library(kSamples)
library(corrplot)
library(VineCopula)
#library(ggplot2)
set.seed(417)
setwd("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Multivarite Analysis/with SM")
data=read_excel("ams_cbridge.xlsx")

colnames(data)={c('Year','SF','Surge','RF','SM')}
tau=cor(data[,2:5], method="kendall")
#Normalize  the data

normalize=function(data){
  data_n=(data-min(data))/(max(data)-min(data))
}
# normalizing the data
data_n=matrix(data=NA, nrow=nrow(data), ncol=(ncol(data)-1)) # frist column is year
for(i in 1:(ncol(data)-1)){
  
  data_n[,i]=normalize(as.matrix(data[,(i+1)]))
}
colnames(data_n)={c('SF','Surge','RF','SM')}

# checking whether dependence is changed bcoz of normalization
s=rowSums(data_n[,2:4])
data_in=data_n[,2:4] # SF is dependent variable
data_ternary=matrix(data=NA, nrow=nrow(data_in), ncol=(ncol(data_in)))
for (j in 1:nrow(data_in)){
  s_temp=s[j]
  data_ternary[j,]=as.matrix(data_in[j,])/s_temp
  
}

## input data for ternany plot should be dataframe
data_ternary=as.data.frame(data_ternary)
colnames(data_ternary)=c('Surge','RF','SM')

### ternary plot
fig_t=plot_ly(data=data_ternary, a=data_ternary[,1],b=data_ternary[,2],c=data_ternary[,3],
        type="scatterternary",
        mode="markers",marker=list(size=15,color=data$SF, colorscale="Blues", showscale=TRUE, reversescale=TRUE),alpha = 1.0) 
fig_t=fig_t%>%  layout(title="",showlegend=F,
         xaxis=list(title="", showgrid=F, zeroline=F, showticklabels=F),
         yaxis=list(title="", showgrid=F, zeroline=F, showticklabels=F),
         sum=1,
         ternary=list(
           aaxis=list(title="Surge", tickformat=".0", tickfont=list(size=15)),
           baxis = list(title = "RainFall", tickformat = ".0", tickfont = list(size = 15)),
           caxis = list(title = "Soil Moisture", tickformat = ".0", tickfont = list(size = 15))))
         #annotations = list(
          # list(xref = "paper", yref = "paper", align = "center",
          #      x = 0.1, y = 1, text = "", ax = 0, ay = 0,
          #      font = list(family = "serif", size = 20, color = "white"),
          #      bgcolor = "#b3b3b3", bordercolor = "black", borderwidth = 2)))

fig_t

## Model fitting
# using source code
setwd("F:/Research Work/Synchronization work/work_R/MultivariateIDFwork")
source("univariate_fit.R")
uni_fit=list()
for(i in 1:ncol(data_n)){
  da=as.matrix(data_n[,i])
  listtemp=marginal_fit(da)
  uni_fit<- append(uni_fit, list(listtemp))
}
names(uni_fit)=c("sf","surge","rf","sm")


# bivariate fitting
# correlation matrix
tau_uk=cor(data_n, method = "kendall")

# correlation test 
k_sfsurge=cor.test(data$Surge,data$SF, method="kendall")
k_sfsm=cor.test(data$SF,data$SM,method="kendall")
k_sfrf=cor.test(data$SF,data$RF,method="kendall")
k_rfsurge=cor.test(data$Surge,data$RF,method="kendall")
k_rfsm=cor.test(data$SM,data$RF,method="kendall")
k_smsurge=cor.test(data$Surge,data$SM,method="kendall")

# all the correlations are statistically significant
# use
# corrplot(tau_uk,method='number')

source("bi_return.R")
bi_data=list(cbind(data_n[,1],data_n[,2]),cbind(data_n[,1],data_n[,3]),cbind(data_n[,1],data_n[,4]))
names(bi_data)=c("sfsurge","sfrf","sfsm")
kendall.Gumbel <- function(copula, t){
  par = copula@parameters
  
  kt <- rep(NA,length(t))
  kt <- t - t * log(t)/(par)
  kt[t==1] <- 1
  kt[t==0] <- 0
  return(kt)  
}

bi_fit=list()
for(i in 1:length(bi_data)){
  data_temp=as.matrix(bi_data[[i]])
  bi_list_temp=bi_return(data_temp)
  bi_fit=append(bi_fit,list(bi_list_temp))
  
}


names(bi_fit)=c("sfsurge","sfrf","sfsm")
### plotting empirical kendall distribution and survival kendall copula function
pt=bi_fit[["sfsurge"]][["kendall"]]

plot(as.vector(pt$f_ecdf),as.vector(pt$ken_cdf), pch=21, col="black",main="Probability Plot",
     xlab="Empirical",ylab="Model")
abline(0,1)
### 
pt_sfrf=bi_fit[["sfrf"]][["kendall"]]
plot(as.vector(pt_sfrf$f_ecdf),as.vector(pt_sfrf$ken_cdf), pch=21, col="black",main="Probability Plot",
     xlab="Empirical",ylab="Model")
abline(0,1)

###
pt_sfsm=bi_fit[["sfsm"]][["kendall"]]
plot(as.vector(pt_sfsm$f_ecdf),as.vector(pt_sfsm$ken_cdf), pch=21, col="black",main="Probability Plot",
     xlab="Empirical",ylab="Model")
abline(0,1)
####
#goodness of fitting
uni_data=cbind(uni_fit[["sf"]][["CDF"]],uni_fit[["surge"]][["CDF"]])

ht_sfrf=BiCopGofTest(uni_fit[["sf"]][["CDF"]],uni_fit[["rf"]][["CDF"]], family = 4,method = "white")
ht_sfsurge=BiCopGofTest(uni_fit[["sf"]][["CDF"]],uni_fit[["surge"]][["CDF"]], family = 4,method = "white")
ht_sfsm=BiCopGofTest(uni_fit[["sf"]][["CDF"]],uni_fit[["sm"]][["CDF"]], family = 4,method = "white")


### scatter
rp_and_obs=cbind(uni_fit[["sf"]][["RP"]],bi_fit[["sfrf"]][["Rp_bi_obs"]][["rt_bi_and"]],bi_fit[["sfsurge"]][["Rp_bi_obs"]][["rt_bi_and"]],bi_fit[["sfsm"]][["Rp_bi_obs"]][["rt_bi_and"]])
data_out=cbind(data[,1:2],rp_and_obs)
colnames(data_out)=c("year","sf","rp_sf","sfrf","sfsurge","sfsm")
write.csv(data_out,file="F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Multivarite Analysis/with SM/rp_sf_obs.csv")

## Plotting
## bivaraite distribution
p=seq(0.005,0.995,length=500)
sf_sim=qgev(p,loc=as.numeric(uni_fit[["sf"]][["parameters"]][1]),scale=as.numeric(uni_fit[["sf"]][["parameters"]][2]),shape=as.numeric(uni_fit[["sf"]][["parameters"]][3]))
surge_sim=qgev(p,loc=as.numeric(uni_fit[["surge"]][["parameters"]][1]),scale=as.numeric(uni_fit[["surge"]][["parameters"]][2]),shape=as.numeric(uni_fit[["surge"]][["parameters"]][3]))
rf_sim=qgev(p,loc=as.numeric(uni_fit[["rf"]][["parameters"]][1]),scale=as.numeric(uni_fit[["rf"]][["parameters"]][2]),shape=as.numeric(uni_fit[["rf"]][["parameters"]][3]))
sm_sim=qgev(p,loc=as.numeric(uni_fit[["sm"]][["parameters"]][1]),scale=as.numeric(uni_fit[["sm"]][["parameters"]][2]),shape=as.numeric(uni_fit[["sm"]][["parameters"]][3]))
## bivariate distribution
# univariate CDF
sf_cdf=pgev(p,loc=as.numeric(uni_fit[["sf"]][["parameters"]][1]),scale=as.numeric(uni_fit[["sf"]][["parameters"]][2]),shape=as.numeric(uni_fit[["sf"]][["parameters"]][3]))
surge_cdf=pgev(p,loc=as.numeric(uni_fit[["surge"]][["parameters"]][1]),scale=as.numeric(uni_fit[["surge"]][["parameters"]][2]),shape=as.numeric(uni_fit[["surge"]][["parameters"]][3]))
rf_cdf=pgev(p,loc=as.numeric(uni_fit[["rf"]][["parameters"]][1]),scale=as.numeric(uni_fit[["rf"]][["parameters"]][2]),shape=as.numeric(uni_fit[["rf"]][["parameters"]][3]))
sm_cdf=pgev(p,loc=as.numeric(uni_fit[["sm"]][["parameters"]][1]),scale=as.numeric(uni_fit[["sm"]][["parameters"]][2]),shape=as.numeric(uni_fit[["sm"]][["parameters"]][3]))
# copula parameter
## sf-surge combination
# gumbel model
gc=bi_fit[["sfsurge"]][["model"]]
# bivariate distribution function
data_plt=matrix(data=NA,nrow=25000,ncol=3)## vector form
dd=expand.grid(SF=sf_sim,Surge=surge_sim)
## sf-surge combination
F_sf_surge=matrix(data=NA, nrow=500,ncol=500)
ken_F_sf_surge=matrix(data=NA, nrow=500,ncol=500)
rp_sf_surge=matrix(data=NA, nrow=500,ncol=500)
F_and=matrix(data=NA, nrow=500,ncol=500)
#x1_temp=matrix(data=NA, nrow=500, ncol=500) ## RF and soil moisture 
#x2_temp=matrix(data=NA, nrow=500, ncol=500) 
for (i in 1:500 ){
 
  for (j in 1:500){
    tp=BiCopCDF(sf_cdf[i],surge_cdf[j],modl)
    F_sf_surge[i,j]=tp
    #and=1-sf_cdf[i]-surge_cdf[j]+tp
    kp=kendall.Gumbel(modl,tp)
    #ken_F_sf_surge[i,j]=kp
    #F_and[i,j]=and
    #rp_sf_surge[i,j]=1/and
   
    #x1_temp[i,j]=U1_sim[i,1]
    #x2_temp[i,j]=U2_sim[j,1]
    
  }
  
}

FF=as.vector(t(ken_F_sf_surge))
dd["CDF"]=FF
### observed bi_Cdf
obs=cbind(uni_fit[["sf"]][["CDF"]],uni_fit[["surge"]][["CDF"]],bi_fit[["sfsurge"]][["cdf"]][["p_and"]])
obs=as.data.frame(obs)
obs[obs<=0.07]=0.007
colnames(obs)=c("SF","Surge","CDF")
## denormalizing the data
x=c(0,0.2,0.4,0.6,0.8,1.0)
data=as.matrix(data)
data_sf=x*(max(data[,2])-min(data[,2]))+min(data[,2])
data_surge=x*(max(data[,3])-min(data[,3]))+min(data[,3])
## bivariate contour plot
fig1=plot_ly(x=sf_cdf,y=surge_cdf,z=as.numeric(FF),type="contour",autocontour="False",ncontours=2, colorscale="RdBu",contours=list(start=0.01,end=0.02,showlabels=TRUE,showlines=TRUE,coloring="lines",line=list(width=2)))%>%
  add_contour(ncontours=2, contours=list(start=0.04,end=0.1, showlabels=TRUE,showlines=TRUE,coloring="lines"))%>%
  add_contour(ncontours=2, contours=list(start=0.2,end=0.5, showlabels=TRUE,showlines=TRUE,coloring="lines"))
fig1=fig1%>%add_trace(x=as.vector(obs$SF),y=as.vector(obs$Surge), type="scatter", mode="markers", color=I("black"))%>%
  hide_guides()
fig1=fig1%>%layout(
  xaxis=list(ticktext=list(134,213,292,371,450,529),
             tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
             tickmode="array", title="StreamFlow(Cumec)",titlefont=list(color="rgb(20,20,20)",size=15),
             tickfont=list(size=15)),
  yaxis=list(ticktext=list(0.833,1.011,1.19,1.368,1.547,1.725),
             tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
             tickmode="array",title="Surge(m)",titlefont=list(color="rgb(20,20,20)",size=15),
             tickfont=list(size=15)))
fig1

#####

#surface plot
fig1s <- plot_ly(z=ken_F_sf_surge , type="surface")

fig1s=fig1s%>% add_surface()
fig1s=fig1s%>%
  layout(
    xaxis=list(ticktext=list(134,213,292,371,450,529),
               tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
               tickmode="array", title="StreamFlow(Cumec)"),
    yaxis=list(ticktext=list(0.833,1.011,1.19,1.368,1.547,1.725),
               tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
               tickmode="array",title="Surge(m)"))
fig1s
###
# contour plot


# for sub plot use start=0.95, end=0.99995, size=0.001 in contours
fig1c <- plot_ly(data=dd, x=~SF, y=~Surge,z=~CDF , type="contour", contours= list(
                                                                                  showlabels=TRUE, line=list(width=5),
                                                                                  labelfont=list(size=12, color="rgb(255, 255, 255)")))%>%
  add_trace(x=as.vector(obs$SF),y=as.vector(obs$Surge), type="scatter", mode="markers", color=I("black"))
fig1c=fig1c%>%
  layout(
    xaxis=list(ticktext=list(134,213,292,371,450,529),
               tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
               tickmode="array", title="StreamFlow(Cumec)",titlefont=list(color="rgb(20,20,20)",size=15),
               tickfont=list(size=15)),
    yaxis=list(ticktext=list(0.833,1.011,1.19,1.368,1.547,1.725),
               tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
               tickmode="array",title="Surge(m)",titlefont=list(color="rgb(20,20,20)",size=15),
               tickfont=list(size=15))
    )
fig1c

#####
dd2=expand.grid(SF=sf_sim,Surge=surge_sim)
dd2["rp"]=as.vector(t(rp_sf_surge))

fig1_rp=plot_ly(x=sf_cdf,y=surge_cdf,z=rp_sf_surge,type="contour",autocontour="False",ncontours=2, colorscale="RdBu",contours=list(start=200,end=100,showlabels=TRUE,showlines=TRUE,coloring="lines",line=list(width=2)))%>%
  add_contour(ncontours=2, contours=list(start=50,end=25, showlabels=TRUE,showlines=TRUE,coloring="lines"))%>%
  add_contour(ncontours=2, contours=list(start=10,end=5, showlabels=TRUE,showlines=TRUE,coloring="lines"))
fig1_rp=fig1_rp%>%add_trace(x=as.vector(obs$SF),y=as.vector(obs$Surge), type="scatter", mode="markers", color=I("black"))%>%
  hide_guides()
fig1_rp=fig1_rp%>%layout(
  xaxis=list(ticktext=list(134,213,292,371,450,529),
             tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
             tickmode="array", title="StreamFlow(Cumec)",titlefont=list(color="rgb(20,20,20)",size=15),
             tickfont=list(size=15)),
  yaxis=list(ticktext=list(0.833,1.011,1.19,1.368,1.547,1.725),
             tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
             tickmode="array",title="Surge(m)",titlefont=list(color="rgb(20,20,20)",size=15),
             tickfont=list(size=15)))

fig1_rp

### sf-soil moisture combination

# gumbel model
gc_sf_sm=bi_fit[["sfsm"]][["model"]]
# bivariate distribution function
d_sf_sm=expand.grid(SF=sf_sim,SM=sm_sim)
F_sf_sm=matrix(data=NA, nrow=500,ncol=500)
ken_F_sf_sm=matrix(data=NA, nrow=500,ncol=500)
rp_sf_sm=matrix(data=NA, nrow=500,ncol=500)
F_and_sf_sm=matrix(data=NA, nrow=500,ncol=500)
#x1_temp=matrix(data=NA, nrow=500, ncol=500) ## RF and soil moisture 
#x2_temp=matrix(data=NA, nrow=500, ncol=500) 
for (i in 1:500 ){
  
  for (j in 1:500){
    tp=pCopula(cbind(sf_cdf[i],sm_cdf[j]),gc)
    F_sf_sm[i,j]=tp
    and=1-sf_cdf[i]-sm_cdf[j]+tp
    kp=kendall.Gumbel(gc,tp)
    ken_F_sf_sm[i,j]=kp
    F_and_sf_sm[i,j]=and
    rp_sf_sm[i,j]=1/and
    
    #x1_temp[i,j]=U1_sim[i,1]
    #x2_temp[i,j]=U2_sim[j,1]
    
  }
  
}

#FF=as.vector(t(ken_F_sf_surge))
#dd["CDF"]=FF
### observed bi_Cdf
obs_sf_sm=cbind(uni_fit[["sf"]][["CDF"]],uni_fit[["sm"]][["CDF"]],bi_fit[["sfsm"]][["cdf"]][["p_and"]])
obs_sf_sm=as.data.frame(obs)

colnames(obs)=c("SF","SM","CDF")
## denormalizing the data
x=c(0,0.2,0.4,0.6,0.8,1.0)
data=as.matrix(data)
data_sf=x*(max(data[,2])-min(data[,2]))+min(data[,2])
data_sm=x*(max(data[,5])-min(data[,5]))+min(data[,5])
## bivariate contour plot
fig2=plot_ly(x=sf_cdf,y=sm_cdf,z=F_and_sf_sm,type="contour",autocontour="False",ncontours=2, colorscale="RdBu",contours=list(start=0.01,end=0.02,showlabels=TRUE,showlines=TRUE,coloring="lines",line=list(width=2)))%>%
  add_contour(ncontours=2, contours=list(start=0.04,end=0.1, showlabels=TRUE,showlines=TRUE,coloring="lines"))%>%
  add_contour(ncontours=2, contours=list(start=0.2,end=0.5, showlabels=TRUE,showlines=TRUE,coloring="lines"))
fig2=fig2%>%add_trace(x=as.vector(obs_sf_sm$SF),y=as.vector(obs_sf_sm$SM), type="scatter", mode="markers", color=I("black"))%>%
  hide_guides()
fig2=fig2%>%layout(
  xaxis=list(ticktext=list(134,213,292,371,450,529),
             tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
             tickmode="array", title="StreamFlow(Cumec)",titlefont=list(color="rgb(20,20,20)",size=15),
             tickfont=list(size=15)),
  yaxis=list(ticktext=list(0.294,0.318,0.341,0.365,0.388,0.411),
             tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
             tickmode="array",title="SoilMoisture(m/m3)",titlefont=list(color="rgb(20,20,20)",size=15),
             tickfont=list(size=15, color="rgb(20,20,20")))
fig2
## 
fig2_rp=plot_ly(x=sf_cdf,y=sm_cdf,z=rp_sf_sm,type="contour",autocontour="False",ncontours=2, colorscale="RdBu",contours=list(start=200,end=100,showlabels=TRUE,showlines=TRUE,coloring="lines",line=list(width=2)))%>%
  add_contour(ncontours=2, contours=list(start=50,end=25, showlabels=TRUE,showlines=TRUE,coloring="lines"))%>%
  add_contour(ncontours=2, contours=list(start=10,end=5, showlabels=TRUE,showlines=TRUE,coloring="lines"))
fig2_rp=fig2_rp%>%add_trace(x=as.vector(obs_sf_sm$SF),y=as.vector(obs_sf_sm$SM), type="scatter", mode="markers", color=I("black"))%>%
  hide_guides()
fig2_rp=fig2_rp%>%layout(
  xaxis=list(ticktext=list(134,213,292,371,450,529),
             tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
             tickmode="array", title="StreamFlow(Cumec)",titlefont=list(color="rgb(20,20,20)",size=15),
             tickfont=list(size=15)),
  yaxis=list(ticktext=list(0.294,0.318,0.341,0.365,0.388,0.411),
             tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
             tickmode="array",title="SoilMoisture(m/m3)",titlefont=list(color="rgb(20,20,20)",size=15),
             tickfont=list(size=15, color="rgb(20,20,20")))
fig2_rp



##SF-Rf combination

gc_sf_rf=bi_fit[["sfrf"]][["model"]]
# bivariate distribution function
d_sf_rf=expand.grid(SF=sf_sim,RF=rf_sim)
F_sf_rf=matrix(data=NA, nrow=500,ncol=500)
ken_F_sf_rf=matrix(data=NA, nrow=500,ncol=500)
rp_sf_rf=matrix(data=NA, nrow=500,ncol=500)
F_and_sf_rf=matrix(data=NA, nrow=500,ncol=500)
#x1_temp=matrix(data=NA, nrow=500, ncol=500) ## RF and soil moisture 
#x2_temp=matrix(data=NA, nrow=500, ncol=500) 
for (i in 1:500 ){
  
  for (j in 1:500){
    tp=pCopula(cbind(sf_cdf[i],rf_cdf[j]),gc)
    F_sf_rf[i,j]=tp
    and=1-sf_cdf[i]-rf_cdf[j]+tp
    F_and_sf_rf[i,j]=and
    kp=kendall.Gumbel(gc,tp)
    ken_F_sf_rf[i,j]=kp
    rp_sf_rf[i,j]=1/and
    
    #x1_temp[i,j]=U1_sim[i,1]
    #x2_temp[i,j]=U2_sim[j,1]
    
  }
  
}
#FF_sf_rf=as.vector(t(ken_F_sf_rf))
#d_sf_rf["CDF"]=FF_sf_rf
### observed bi_Cdf
obs_sf_rf=cbind(uni_fit[["sf"]][["CDF"]],uni_fit[["rf"]][["CDF"]],bi_fit[["sfrf"]][["cdf"]][["p_and"]])
obs_sf_rf=as.data.frame(obs)

colnames(obs)=c("SF","RF","CDF")
## bivariate contour plot
## denormalizing the data
x=c(0,0.2,0.4,0.6,0.8,1.0)
data_sf=x*(max(data[,2])-min(data[,2]))+min(data[,2])
data_rf=x*(max(data[,4])-min(data[,4]))+min(data[,4])


### 
## bivariate contour plot
fig3=plot_ly(x=sf_cdf,y=rf_cdf,z=F_and_sf_rf,type="contour",autocontour="False",ncontours=2, colorscale="RdBu",contours=list(start=0.01,end=0.02,showlabels=TRUE,showlines=TRUE,coloring="lines",line=list(width=2)))%>%
  add_contour(ncontours=2, contours=list(start=0.04,end=0.1, showlabels=TRUE,showlines=TRUE,coloring="lines"))%>%
  add_contour(ncontours=2, contours=list(start=0.2,end=0.5, showlabels=TRUE,showlines=TRUE,coloring="lines"))
fig3=fig3%>%add_trace(x=as.vector(obs_sf_rf$SF),y=as.vector(obs_sf_rf$RF), type="scatter", mode="markers", color=I("black"))%>%
  hide_guides()
fig3=fig3%>%layout(
  xaxis=list(ticktext=list(134,213,292,371,450,529),
             tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
             tickmode="array", title="StreamFlow(Cumec)",titlefont=list(color="rgb(20,20,20)",size=15),
             tickfont=list(size=15)),
  yaxis=list(ticktext=list(17.375,22.765,28.155,33.545,38.935,44.325),
             tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
             tickmode="array",title="Rainfall(mm)",titlefont=list(color="rgb(20,20,20)",size=15),
             tickfont=list(size=15, color="rgb(20,20,20")))
fig3
## 
fig3_rp=plot_ly(x=sf_cdf,y=rf_cdf,z=rp_sf_rf,type="contour",autocontour="False",ncontours=2, colorscale="RdBu",contours=list(start=200,end=100,showlabels=TRUE,showlines=TRUE,coloring="lines",line=list(width=2)))%>%
  add_contour(ncontours=2, contours=list(start=50,end=25, showlabels=TRUE,showlines=TRUE,coloring="lines"))%>%
  add_contour(ncontours=2, contours=list(start=10,end=5, showlabels=TRUE,showlines=TRUE,coloring="lines"))
fig3_rp=fig3_rp%>%add_trace(x=as.vector(obs_sf_rf$SF),y=as.vector(obs_sf_rf$RF), type="scatter", mode="markers", color=I("black"))%>%
  hide_guides()
fig3_rp=fig3_rp%>%layout(
  xaxis=list(ticktext=list(134,213,292,371,450,529),
             tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
             tickmode="array", title="StreamFlow(Cumec)",titlefont=list(color="rgb(20,20,20)",size=15),
             tickfont=list(size=15)),
  yaxis=list(ticktext=list(17.375,22.765,28.155,33.545,38.935,44.325),
             tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
             tickmode="array",title="Rainfall(mm)",titlefont=list(color="rgb(20,20,20)",size=15),
             tickfont=list(size=15, color="rgb(20,20,20")))
fig3_rp



## bivariate contour plot
#surface plot
#fig2s <- plot_ly(z=ken_F_sf_surge , type="surface")

#fig2s=fig2s%>% add_surface()
#fig2s=fig2s%>%
#  layout(
#    xaxis=list(ticktext=list("134","213","292","371","450","529"),
#               tickvals=list(0,100,200,300,400,500),
#               tickmode="array", title="StreamFlow(Cumec)"),
#    yaxis=list(ticktext=list(0.833,1.011,1.19,1.368,1.547,1.725),
#               tickvals=list(0,100,200,300,400,500),
#              tickmode="array",title="Soil Moisture(m/m3"))
#fig2s
# contour plot
## to create a sub plot change autocontour=False start=0.95, end=0.99995, size=0.001
fig3c <- plot_ly(data=d_sf_rf, x=~SF, y=~RF ,z=~CDF , type="contour", contours= list(showlabels=TRUE,
                                                                                     labelfont=list(size=15, color="rgb(255,255,255)")))%>%
  add_trace(x=as.vector(obs_sf_rf$SF),y=as.vector(obs_sf_rf$RF), type="scatter", mode="markers", color=I("black"))
fig3c=fig3c%>%
  layout(
    xaxis=list(ticktext=list(134,213,292,371,450,529),
               tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
               tickmode="array", title="StreamFlow(Cumec)",titlefont=list(color="rgb(20,20,20)",size=15),
               tickfont=list(size=15, color="rgb(20,20,20")),
    yaxis=list(ticktext=list(17.375,22.765,28.155,33.545,38.935,44.325),
               tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
               tickmode="array",title="Rainfall(mm)",titlefont=list(color="rgb(20,20,20)",size=15),
               tickfont=list(size=15, color="rgb(20,20,20")))
fig3c


dd2_sf_rf=expand.grid(SF=sf_sim,RF=rf_sim)
dd2_sf_rf["rp"]=as.vector(t(rp_sf_rf))

fig3rp <- plot_ly(data=dd2_sf_rf, x=~SF, y=~RF ,z=~rp , type="contour", autocontour=F, contours= list(start=2, end=220, size=30,showlabels=TRUE,
                                                                                                      labelfont=list(size=15, color="rgb(255,255,255)")), colorscale="hot")%>%
  add_trace(x=as.vector(obs_sf_rf$SF),y=as.vector(obs_sf_rf$RF), type="scatter", mode="markers", color=I("black"))
fig3rp=fig3rp%>%
  layout(
    xaxis=list(ticktext=list(134,213,292,371,450,529),
               tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
               tickmode="array", title="StreamFlow(Cumec)",titlefont=list(color="rgb(20,20,20)",size=15),
               tickfont=list(size=15, color="rgb(20,20,20")),
    yaxis=list(ticktext=list(17.375,22.765,28.155,33.545,38.935,44.325),
               tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
               tickmode="array",title="Rainfall(mm)",titlefont=list(color="rgb(20,20,20)",size=15),
               tickfont=list(size=15, color="rgb(20,20,20")))
fig3rp=fig3rp%>%colorbar(title="Return \n Period(Yrs)")
fig3rp
