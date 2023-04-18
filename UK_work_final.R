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


# export the data 
rt_uni=as.data.frame(cbind(uni_fit$sm$RP,uni_fit$rf$RP,uni_fit$surge$RP,uni_fit$sf$RP))
colnames(rt_uni)={c('sm','rf','Surge','sf')}
write.csv(rt_uni, "F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Data/Data_at_cbridge/return period_univariate.csv")

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


## copula fit
#uni_data=cbind(uni_fit[["sf"]][["CDF"]],uni_fit[["surge"]][["CDF"]])

bi_data=list(cbind(uni_fit[["sf"]][["CDF"]],uni_fit[["rf"]][["CDF"]]),
             cbind(uni_fit[["sf"]][["CDF"]],uni_fit[["surge"]][["CDF"]]),
             cbind(uni_fit[["sf"]][["CDF"]],uni_fit[["sm"]][["CDF"]]))
names(bi_data)=c("sfrf","sfsurge","sfsm")

source("copula_fit.r")
source("kendall_functions.r")


bi_fit=list()
kt_dist=list()
rt=list()
for(i in 1:3){
  data_temp=bi_data[[i]]
  fit_temp=bifit(data_temp)
  model=fit_temp$model
  t=fit_temp[["simulated cdf"]][,4]
  kt_e=kendall_dist(model,t)
  bi_fit=append(bi_fit,list(fit_temp))
  rt_f=1/(1-kt_e)
  kt_dist=append(kt_dist,list(kt_e))
  rt=append(rt,list(rt_f))
}

names(rt)=c("sfrf","sfsurge","sfsm")
names(bi_fit)=c("sfrf","sfsurge","sfsm")
## exporting data
rt_sf_surge=data.frame(bi_fit[["sfsurge"]][["simulated cdf"]][,1:2],rt["sfsurge"])


fig1_rp=plot_ly(x=rt_sf_surge$u1,y=rt_sf_surge$u2,z=rt_sf_surge$sfsurge,type="contour",autocontour="True", colorscale="RdBu",contours=list(showlabels=TRUE,showlines=TRUE,coloring="lines",line=list(width=2)))%>%
  add_contour(ncontours=2, contours=list(start=50,end=25, showlabels=TRUE,showlines=TRUE,coloring="lines"))%>%
  add_contour(ncontours=2, contours=list(start=10,end=5, showlabels=TRUE,showlines=TRUE,coloring="lines"))

source("returnperiods.r")

write.csv(rt_sf_surge,"F:/Research Work/Synchronization work/Multivariate IDF/Draft/Results/Bivariate results/rt_sfsurge.csv")
 


