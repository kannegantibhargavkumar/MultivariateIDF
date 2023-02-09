# This file contains procedure to estimate bi/multi return periods. 
# Input data: Rainfall and Streamflow
# 
library(readxl)
library(evd)
library(copula)
library(plotly)

ams=read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Data/Data_at_cbridge/ams_cbridge.xlsx")
colnames(ams)=c("time","SF","RF")

# model fitting
setwd("F:/Research Work/Synchronization work/work_R/MultivariateIDFwork")
source("univariate_fit.R")
uni_fit=list()
for(i in 1:(ncol(ams)-1)){
  da=as.matrix(ams[8:70,i+1])
  listtemp=marginal_fit(da)
  uni_fit<- append(uni_fit, list(listtemp))
}
names(uni_fit)=c("SF","RF")

# removing NAs
ams=na.omit(ams)
tau=cor(ams$SF,ams$RF, method = "kendall")

# bivariate copula model fit
source("bi_return.R")
bi=bi_return(as.matrix(ams))

## calculating the PCP_0.8, SF_0.8
sf_0.8=qgev(0.8,loc=as.numeric(uni_fit[["SF"]][["parameters"]][1]),scale=as.numeric(uni_fit[["SF"]][["parameters"]][2]),shape=as.numeric(uni_fit[["SF"]][["parameters"]][3]))
rf_0.8=qgev(0.8,loc=as.numeric(uni_fit[["RF"]][["parameters"]][1]),scale=as.numeric(uni_fit[["RF"]][["parameters"]][2]),shape=as.numeric(uni_fit[["RF"]][["parameters"]][3]))
# plotting 
axx= list(ticktext=list(100,"SF_0.8" ,600),
          tickvals=list(100, 332,600),
          tickmode="array", title="StreamFlow(Cumec)",linecolor="black",linewidth=2, mirror=T)
axy=list(ticktext=list(15,"RF_0.8",60),
         tickvals=list(15,34, 60),
         tickmode="array", title="Rainfall(mm)",linecolor="black",linewidth=2, mirror=T)
hline=function(y=0, color="black"){
  list(type="line",
       x0=100,
       x1=600,
       xref="paper",
       y0=y,
       y1=y,
       line=list(color=color))
}
vline <- function(x = 0, color = "green") {
  list(
    type = "line",
    y0 = 15,
    y1 = 60,
    yref = "paper",
    x0 = x,
    x1 = x,
    line = list(color = color, dash="dot")
  )
}
fig=plot_ly(data=ams, x=ams$SF, y=ams$RF, type="scatter", mode="markers",marker=list(size=8), color = I("black"))
fig=fig%>% layout(plot_bgcolor="white",xaxis=axx, yaxis=axy,shapes=list(hline(34),vline(332), list(type="rect", fillcolor="grey",opacity=0.4, line=list(color="grey"),x0= 332, x1= 620, y0=34, y1=60)))
fig

## or proabability
fig2=plot_ly(data=ams, x=ams$SF, y=ams$RF, type="scatter", mode="markers",marker=list(size=8), color = I("black"))
fig2=fig2%>% layout(plot_bgcolor="white",xaxis=axx, yaxis=axy,shapes=list(hline(34),vline(332), list(type="rect", fillcolor="grey",opacity=0.3, line=list(color="grey"),x0=95, x1= 620, y0=34, y1=60),
                                                                          list(type="rect", fillcolor="grey",opacity=0.3, line=list(color="grey"),x0=332, x1= 620, y0=15, y1=34)))
fig2
##  condition probability

## fig3 cond1
fig3=plot_ly(data=ams, x=ams$SF, y=ams$RF, type="scatter", mode="markers",marker=list(size=8), color = I("black"))
fig3=fig3%>% layout(plot_bgcolor="white",xaxis=axx, yaxis=axy,shapes=list(hline(34),vline(332), list(type="rect", fillcolor="grey",opacity=0.3, line=list(color="grey"),x0=332, x1= 620, y0=34, y1=60)))
fig3
