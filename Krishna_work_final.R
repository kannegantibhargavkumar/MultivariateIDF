#####
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
#library(ggplot2)
set.seed(417)
#####
#Importing the dataset
setwd("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/data_input")
ams_warunji=read_excel("ams_warunji.xlsx") # the annual maxima series
# calculating the dependence between variables

ams=subset(ams_warunji,select = c("sm","rf","sf"))
ams[24,1]=c(0.462)
ams[24,2]=c(124.25268)
ams[24,3]=c(9378.232)

normalize=function(data){
  data_n=(data-min(data))/(max(data)-min(data))
}
# normalizing the data
data_n=matrix(data=NA, nrow=35, ncol=3)
for(i in 1:3){
 
  data_n[,i]=normalize(as.matrix(ams[,i]))
}


#check this later
data_n=as.data.frame(data_n)
colnames(data_n)=c("SoilMoisture","Rainfall", "Streamflow")
## kendall correlation
tau=cor(data_n, method="kendall")

### plots
######
par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(ams$sf, ams$sm, pch = 16, col = 2, xlab="", ylab="")              # Create first plot
par(new = TRUE)                             # Add new plot
plot(ams$sf, ams$rf, pch = 17, col = 3,              # Create second plot without axes
     axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(ams$rf)))      # Add second axis
mtext("Rainfall(mm)", side = 4, line = 3)
mtext("Streamflow(cumec)", side=1, line=3) # Add second axis label
mtext("SoilMoisture(mm)", side=2, line=3)
legend(5530,125, legend=c("SoilMoisture","Rainfall"), col=c("Red","Green"), lty=1:2, cex=0.8)
### scatter with contour 

fig=plot_ly(data=ams, x=~sf,y=~sm, size=~rf, type="scatter", mode="markers")
fig=fig%>% plot_ly(data = ams , x=~sf,y=~sm, z=~rf, type = "contour",contours=list(fillcolor='None'))
fig

plot_ly(data = ams , x=~sf,y=~sm, z=~rf, type = "contour", contours=list(fillcolor='None', showlables= TRUE)) %>%
  add_trace(x=~sf,y=~sm, size=~rf, type = "scatter",mode="markers", color=I("black"))
## 3d scatter plot
axx=list(title="Soil Moisture(mm)")
axy=list(title="Precipitation(mm)")
axz=list(title="Streamflow(Cumec)")

fig=plot_ly(x=as.vector(ams$sm),y=as.vector(ams$rf), z=as.vector(ams$sf),type="scatter3d", mode="markers", color=z)
fig=fig%>% layout(scene=list(xaxis=axx, yaxis=axy, zaxis=axz))
fig

s=rowSums(data_n[,1:2])
data_ternary=matrix(data=NA, nrow=35, ncol=2)
for (j in 1:35){
  s_temp=s[j]
  data_ternary[j,]=as.matrix(data_n[j,1:2])/s_temp
  
}
## input data for ternany plot should be dataframe
data_ternary=as.data.frame(data_ternary)
colnames(data_ternary)=c("SoilMoisture","Rainfall")


# scatter plot
plot_ly(data=data_ternary, x=~SoilMoisture, y=~Rainfall, size=20,color=~as.vector(ams$sf), type="scatter", mode="markers", colors="Blues", alpha=1.0)%>%
  layout(title="", showlegend=F,
         xaxis=list(title="Soil Moisture m/m3",titlefont=list(color="rgb(20,20,20)",size=15),
                    tickfont=list(size=15)),
         yaxis=list(title="Rainfall(mm)",titlefont=list(color="rgb(20,20,20)",size=15),
                   tickfont=list(size=15)))




### ternary plot
plot_ly(data=data_ternary, a=data_ternary[,1],b=data_ternary[,2],c=data_ternary[,3],
        type="scatterternary",
        mode="markers",
        evaluate=T) %>%
  layout(title="",showlegend=F,
         xaxis=list(title="", showgrid=F, zeroline=F, showticklabels=F),
         yaxis=list(title="", showgrid=F, zeroline=F, showticklabels=F),
         sum=1,
         ternary=list(
           aaxis=list(title="SoilMoisture", tickformat=".0", tickfont=list(size=10)),
           baxis = list(title = "Rainfall", tickformat = ".0", tickfont = list(size = 10)),
           caxis = list(title = "Streamflow", tickformat = ".0", tickfont = list(size = 10))),
         annotations = list(
           list(xref = "paper", yref = "paper", align = "center",
           x = 0.1, y = 1, text = "Ternary Plot in R(Markers)", ax = 0, ay = 0,
           font = list(family = "serif", size = 15, color = "white"),
           bgcolor = "#b3b3b3", bordercolor = "black", borderwidth = 2)))

## using ggplot
ggtern(data=data_n, aes(x=Rainfall, y=SoilMoisture, z=Streamflow))+
  geom_point()+
  labs(title="Dependence between Hydrologic variables")+
  theme_rgbw()

###### 
# multiple linear regression modelling 


## model fitting

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
names(uni_fit)=c("sm","rf","sf")


# bivariate fitting
# correlation matrix
tau_uk=cor(data_n, method = "kendall")

# correlation test 
k_sfsm=cor.test(ams$sm,ams$sf,method="kendall")
k_sfrf=cor.test(ams$sf,ams$rf,method="kendall")
k_rfsm=cor.test(ams$sm,ams$rf,method="kendall")

## copula fit
#uni_data=cbind(uni_fit[["sf"]][["CDF"]],uni_fit[["surge"]][["CDF"]])

bi_data=list(cbind(uni_fit[["sf"]][["CDF"]],uni_fit[["rf"]][["CDF"]]),
             cbind(uni_fit[["sf"]][["CDF"]],uni_fit[["sm"]][["CDF"]]))
names(bi_data)=c("sfrf","sfsm")

source("copula_fit.r")
source("kendall_functions.r")


bi_fit=list()
kt_dist=list()
rt=list()
for(i in 1:2){
  data_temp=bi_data[[i]]
  fit_temp=bifit(data_temp)
  model=fit_temp$model
  t=fit_temp[["observed cdf"]][,4]
  kt_e=kendall_dist(model,t)
  bi_fit=append(bi_fit,list(fit_temp))
  rt_f=1/(1-kt_e)
  kt_dist=append(kt_dist,list(kt_e))
  rt=append(rt,list(rt_f))
}

names(rt)=c("sfrf","sfsm")
names(bi_fit)=c("sfrf","sfsm")
## exporting data

write.csv(uni_fit[["sf"]][["RP"]], "F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/RP_sF_UNI.csv")





































# using source code
setwd("F:/Research Work/Synchronization work/work_R/MultivariateIDFwork")
source("univariate_fit.R")
uni_fit=list()
for(i in 1:ncol(data_n)){
  data=data_n[,i]
  listtemp=marginal_fit(data)
  uni_fit<- append(uni_fit, list(listtemp))
}
names(uni_fit)=c("sm","rf","sf")

# bivariate fitting
# correlation matrix
tau=cor(ams,method="kendall")
# correlation test 
#cor.test(as.matrix(ams[,1]),as.matrix(ams[,2]))
# all the correlations are statistically significant
source("bi_return.R")
bi_data=list(cbind(data_n[,1],data_n[,2]),cbind(data_n[,1],data_n[,3]),cbind(data_n[,2],data_n[,3]))
names(bi_data)=c("smrf","smsf","rfsf")
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
names(bi_fit)=c("smrf","smsf","rfsf")

## plotting
## bivariate distribution modelling
p=seq(0.005,0.995,length=500)
p1=seq(0.005,0.995, length=500)
p2=seq(0.005,0.995, length=500)
sf_sim=qgev(p,loc=as.numeric(uni_fit[["sf"]][["parameters"]][1]),scale=as.numeric(uni_fit[["sf"]][["parameters"]][2]),shape=as.numeric(uni_fit[["sf"]][["parameters"]][3]))
rf_sim=qgev(p2,loc=as.numeric(uni_fit[["rf"]][["parameters"]][1]),scale=as.numeric(uni_fit[["rf"]][["parameters"]][2]),shape=as.numeric(uni_fit[["rf"]][["parameters"]][3]))
sm_sim=qgev(p1,loc=as.numeric(uni_fit[["sm"]][["parameters"]][1]),scale=as.numeric(uni_fit[["sm"]][["parameters"]][2]),shape=as.numeric(uni_fit[["sm"]][["parameters"]][3]))
## univariate cdf
sf_cdf=pgev(sf_sim,loc=as.numeric(uni_fit[["sf"]][["parameters"]][1]),scale=as.numeric(uni_fit[["sf"]][["parameters"]][2]),shape=as.numeric(uni_fit[["sf"]][["parameters"]][3]))
rf_cdf=pgev(rf_sim,loc=as.numeric(uni_fit[["rf"]][["parameters"]][1]),scale=as.numeric(uni_fit[["rf"]][["parameters"]][2]),shape=as.numeric(uni_fit[["rf"]][["parameters"]][3]))
sm_cdf=pgev(sm_sim,loc=as.numeric(uni_fit[["sm"]][["parameters"]][1]),scale=as.numeric(uni_fit[["sm"]][["parameters"]][2]),shape=as.numeric(uni_fit[["sm"]][["parameters"]][3]))
# copula parameter
## sf-sm combination
# gumbel model
gc=bi_fit[["smsf"]][["model"]]
# bivariate distribution function
data_plt=matrix(data=NA,nrow=25000,ncol=3)## vector form
dd=expand.grid(SF=sf_sim,SM=sm_sim) ## use normalized values
## sf-surge combination
F_sf_sm=matrix(data=NA, nrow=500,ncol=500)
ken_F_sf_sm=matrix(data=NA, nrow=500,ncol=500)
rp_sf_sm=matrix(data=NA, nrow=500,ncol=500)
F_and=matrix(data=NA, nrow=500,ncol=500)
rp_and=matrix(data=NA, nrow=500,ncol=500)
#x1_temp=matrix(data=NA, nrow=500, ncol=500) ## RF and soil moisture 
#x2_temp=matrix(data=NA, nrow=500, ncol=500) 
for (i in 1:500 ){
  
  for (j in 1:500){
    tp=pCopula(cbind(sf_cdf[i],sm_cdf[j]),gc)
    F_sf_sm[i,j]=tp
    and=1-sf_cdf[i]-sm_cdf[j]+tp;
    kp=kendall.Gumbel(gc,tp)
    ken_F_sf_sm[i,j]=kp
    rp_sf_sm[i,j]=1/(1-kp)
    rp_and[i,j]=1/and
    F_and[i,j]=and
    
    #x1_temp[i,j]=U1_sim[i,1]
    #x2_temp[i,j]=U2_sim[j,1]
    
  }
  
}

FF=as.vector(t(ken_F_sf_sm))
dd["CDF"]=FF
### observed bi_Cdf
obs=cbind(uni_fit[["sm"]][["CDF"]],uni_fit[["sf"]][["CDF"]],bi_fit[["smsf"]][["cdf"]][["p_and"]])
obs=as.data.frame(obs)
#obs[obs<=0.07]=0.007
colnames(obs)=c("SM","SF","CDF")
## denormalizing the data
x=c(0,0.2,0.4,0.6,0.8,0.995)
data_sf=x*(max(ams$sf)-min(ams$sf))+min(ams$sf)
data_sm=x*(max(ams$sm)-min(ams$sm))+min(ams$sm)

# contour plot
fig1=plot_ly(x=sf_cdf,y=sm_cdf,z=F_and,type="contour",autocontour="False",ncontours=2, colorscale="RdBu",contours=list(start=0.005,end=0.01,showlabels=TRUE,showlines=TRUE,coloring="lines",line=list(width=2)))%>%
  add_contour(ncontours=2, contours=list(start=0.04,end=0.1, showlabels=TRUE,showlines=TRUE,coloring="lines"))%>%
  add_contour(ncontours=2, contours=list(start=0.2,end=0.5, showlabels=TRUE,showlines=TRUE,coloring="lines"))
fig1=fig1%>%add_trace(x=as.vector(obs$SF),y=as.vector(obs$SM), type="scatter", mode="markers", color=I("black"))%>%
  hide_guides()
fig1=fig1%>%layout(
  xaxis=list(ticktext=list(774.10, 2494.92, 4215.75, 5936.57, 7657.40, 26832),
             tickvals=list(0,0.2,0.4,0.6,0.8,0.995),
             tickmode="array", title="StreamFlow(Cumec)",titlefont=list(color="rgb(20,20,20)",size=15),
             tickfont=list(size=15)),
  yaxis=list(ticktext=list(0.38, 0.39, 0.41, 0.43, 0.45, 0.470),
             tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
             tickmode="array",title="Soil Moisture(m/m3)",titlefont=list(color="rgb(20,20,20)",size=15),
             tickfont=list(size=15)))

  
# for sub plot use start=0.95, end=0.99995, size=0.001 in contours
fig1c <- plot_ly(data=dd, x=~SF, y=~SM,z=~CDF , type="contour", autocontours="True", contours= list(
  showlabels=TRUE, line=list(width=5),
  labelfont=list(size=12, color="rgb(255, 255, 255)")))%>%
  add_trace(x=as.vector(obs$SF),y=as.vector(obs$SM), type="scatter", mode="markers", color=I("black"))
fig1c=fig1c%>%
  layout(
    xaxis=list(ticktext=list(774.10, 2494.926, 4215.753, 5936.579, 7657.406, 9378.232),
               tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
               tickmode="array", title="StreamFlow(Cumec)",titlefont=list(color="rgb(20,20,20)",size=15),
               tickfont=list(size=15)),
    yaxis=list(ticktext=list(0.38, 0.39, 0.41, 0.43, 0.45, 0.46),
               tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
               tickmode="array",title="Soil Moisture(m/m3)",titlefont=list(color="rgb(20,20,20)",size=15),
               tickfont=list(size=15))
  )
fig1c


## return period

dd2=expand.grid(SF=sf_sim,SM=sm_sim)
dd2["rp"]=as.vector(t(rp_sf_sm))
fig1_rp=plot_ly(x=sf_cdf,y=sm_cdf,z=rp_and,type="contour",autocontour="False",ncontours=2, colorscale="RdBu",contours=list(start=200,end=100,showlabels=TRUE,showlines=TRUE,coloring="lines",line=list(width=2)))%>%
  add_contour(ncontours=2, contours=list(start=25,end=10, showlabels=TRUE,showlines=TRUE,coloring="lines"))%>%
  add_contour(ncontours=2, contours=list(start=5,end=2, showlabels=TRUE,showlines=TRUE,coloring="lines"))
fig1_rp=fig1_rp%>%add_trace(x=as.vector(obs$SF),y=as.vector(obs$SM), type="scatter", mode="markers", color=I("black"))%>%
  hide_guides()
fig1_rp=fig1_rp%>%layout(
  xaxis=list(ticktext=list(774.10, 2494.926, 4215.753, 5936.579, 7657.406, 9378.232),
             tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
             tickmode="array", title="StreamFlow(Cumec)",titlefont=list(color="rgb(20,20,20)",size=15),
             tickfont=list(size=15)),
  yaxis=list(ticktext=list(0.38, 0.39, 0.41, 0.43, 0.45, 0.460),
             tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
             tickmode="array",title="Soil Moisture(m/m3)",titlefont=list(color="rgb(20,20,20)",size=15),
             tickfont=list(size=15)))


#start=2, end=220, size=10,
fig1rp <- plot_ly(data=dd2, x=~SF, y=~SM,z=~rp , type="contour", autocontour=TRUE, contours= list(showlabels=TRUE,
                                                                                                  labelfont=list(size=15, color="rgb(255,255,255)")), colorscale="hot")%>%
  add_trace(x=as.vector(obs$SF),y=as.vector(obs$SM), type="scatter", mode="markers", color=I("black"))
fig1rp=fig1rp%>%
  layout(
    xaxis=list(ticktext=list(774.10, 2494.926, 4215.753, 5936.579, 7657.406, 9378.232),
               tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
               tickmode="array", title="StreamFlow(Cumec)",titlefont=list(color="rgb(20,20,20)",size=15),
               tickfont=list(size=15, color="rgb(20,20,20")),
    yaxis=list(ticktext=list(0.38, 0.39, 0.41, 0.43, 0.45, 0.46),
               tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
               tickmode="array",title="Soil Moisture(m/m3)",titlefont=list(color="rgb(20,20,20)",size=15),
               tickfont=list(size=15, color="rgb(20,20,20")))
fig1rp=fig1rp%>%colorbar(title="Return \n Period(Yrs)")
fig1rp

### SF-RF combination
# gumbel model
gc=bi_fit[["rfsf"]][["model"]]
# bivariate distribution function
data_plt=matrix(data=NA,nrow=25000,ncol=3)## vector form
dd_rf=expand.grid(SF=sf_sim,RF=rf_sim) ## use normalized values
## sf-surge combination
F_sf_rf=matrix(data=NA, nrow=500,ncol=500)
ken_F_sf_rf=matrix(data=NA, nrow=500,ncol=500)
rp_sf_rf=matrix(data=NA, nrow=500,ncol=500)
#x1_temp=matrix(data=NA, nrow=500, ncol=500) ## RF and soil moisture 
#x2_temp=matrix(data=NA, nrow=500, ncol=500) 
for (i in 1:500 ){
  
  for (j in 1:500){
    tp=pCopula(cbind(sf_cdf[i],rf_cdf[j]),gc)
    F_sf_rf[i,j]=tp
    kp=kendall.Gumbel(gc,tp)
    ken_F_sf_rf[i,j]=kp
    rp_sf_rf[i,j]=1/(1-kp)
    
    #x1_temp[i,j]=U1_sim[i,1]
    #x2_temp[i,j]=U2_sim[j,1]
    
  }
  
}

FF_rf=as.vector(t(ken_F_sf_rf))
dd_rf["CDF"]=FF_rf
### observed bi_Cdf
obs_rf=cbind(data_n[,2],data_n[,3],bi_fit[["rfsf"]][["ken_cdf"]])
obs_rf=as.data.frame(obs_rf)
#obs[obs<=0.07]=0.007
colnames(obs_rf)=c("RF","SF","CDF")
## denormalizing the data
x=c(0,0.2,0.4,0.6,0.8,1.0)
data_sf=x*(max(ams$sf)-min(ams$sf))+min(ams$sf)
data_rf=x*(max(ams$rf)-min(ams$rf))+min(ams$rf)

# contour plot
# for sub plot use start=0.95, end=0.99995, size=0.001 in contours
fig2c <- plot_ly(data=dd_rf, x=~SF, y=~RF,z=~CDF , type="contour", autocontours="False", contours= list(start=0.92, end=0.99995, size=0.005,
  showlabels=TRUE, line=list(width=5),
  labelfont=list(size=12, color="rgb(255, 255, 255)")))%>%
  add_trace(x=as.vector(obs_rf$SF),y=as.vector(obs_rf$RF), type="scatter", mode="markers", color=I("black"))
fig2c=fig2c%>%
  layout(
    xaxis=list(ticktext=list(774.10, 2494.926, 4215.753, 5936.579, 7657.406, 9378.232),
               tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
               tickmode="array", title="StreamFlow(Cumec)",titlefont=list(color="rgb(20,20,20)",size=15),
               tickfont=list(size=15)),
    yaxis=list(ticktext=list(33.551, 51.691,69.832,87.972,106.112,124.252),
               tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
               tickmode="array",title="Rainfall(mm)",titlefont=list(color="rgb(20,20,20)",size=15),
               tickfont=list(size=15))
  )
fig2c


## return period

dd2_rf=expand.grid(SF=sf_sim,RF=rf_sim)
dd2_rf["rp"]=as.vector(t(rp_sf_rf))
#start=2, end=220, size=10,
fig2rp <- plot_ly(data=dd2_rf, x=~SF, y=~RF,z=~rp , type="contour", autocontour=TRUE, contours= list(showlabels=TRUE,
                                                                                                  labelfont=list(size=15, color="rgb(255,255,255)")), colorscale="hot")%>%
  add_trace(x=as.vector(obs$SF),y=as.vector(obs$RF), type="scatter", mode="markers", color=I("black"))
fig2rp=fig2rp%>%
  layout(
    xaxis=list(ticktext=list(774.10, 2494.926, 4215.753, 5936.579, 7657.406, 9378.232),
               tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
               tickmode="array", title="StreamFlow(Cumec)",titlefont=list(color="rgb(20,20,20)",size=15),
               tickfont=list(size=15, color="rgb(20,20,20")),
    yaxis=list(ticktext=list(33.551, 51.691,69.832,87.972,106.112,124.252),
               tickvals=list(0,0.2,0.4,0.6,0.8,1.0),
               tickmode="array",title="Rainfall(mm)",titlefont=list(color="rgb(20,20,20)",size=15),
               tickfont=list(size=15, color="rgb(20,20,20")))
fig2rp=fig2rp%>%colorbar(title="Return \n Period(Yrs)")
fig2rp

























## conditional distribution functions
### need to use conditional copula
setwd("F:/Research Work/Synchronization work/work_R/MultivariateIDFwork")
source("cond_gumbelcopula.R")
#conditional distributions
cond_sf_rt=list()
alpha_1=bi_fit[["smsf"]][["model"]]@parameters
bi_obs_cdf_1=bi_fit[["smsf"]][["obs_CDF"]]
u_1=uni_fit[["sm"]][["CDF"]]
v_1=uni_fit[["sf"]][["CDF"]]
cond_sf_sm=con_prob(alpha_1,bi_obs_cdf_1,u_1,v_1)
###
alpha_2=bi_fit[["rfsf"]][["model"]]@parameters
bi_obs_cdf_2=bi_fit[["rfsf"]][["obs_CDF"]]
u_2=uni_fit[["rf"]][["CDF"]]
v_2=uni_fit[["sf"]][["CDF"]]
cond_sf_rf=con_prob(alpha_2,bi_obs_cdf_2,u_2,v_2)
###
#Hypothesis testing
source("hypothesis_test.R")
## Sf and Sf conditioned on SM
options(scipen = 999) # to avoid e^ values
t_sf_sm=ht(as.matrix(cbind(uni_fit[["sf"]][["RP"]],cond_sf_sm[["C_rt"]])))
t_sf_rf=ht(as.matrix(cbind(uni_fit[["sf"]][["RP"]],cond_sf_rf[["C_rt"]])))
write.csv(t_sf_sm,file="F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/ht_sm_sf.csv")
write.csv(t_sf_rf,file="F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/ht_rf_sf.csv")
# writing the output
data_to_test_smsf=as.matrix(cbind(uni_fit[["sf"]][["RP"]],cond_sf_sm[["C_rt"]],as.matrix(ams$sf)))
colnames(data_to_test_smsf)=c("RP","Conditional RP","sf")
write.csv(data_to_test_smsf,file="F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/rps_sm_sf_final.csv")

data_to_test_smrf=as.matrix(cbind(uni_fit[["sf"]][["RP"]],cond_sf_rf[["C_rt"]],as.matrix(ams$sf)))
colnames(data_to_test_smrf)=c("RP","Conditional RP","sf")
write.csv(data_to_test_smrf,file="F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/rps_sm_rf_final.csv")

### 
#plot
plot(ecdf(data_to_test_smrf[,1]), main="", lwd=2, col="red", xlab="Return Period(Years)", ylab="Probability")
lines(ecdf(data_to_test_smsf[,2]),lwd=2, col="green")
legend("bottomright",legend=c("RP","Conditioned RP"), col=c("red","green"), lty=1.5, lwd=2)

## bivariate distribution 
p=seq(0.05,0.9995,length=500)

## streamflow and soil moisture
x1_pseudo=qgev(as.matrix(p),loc=as.numeric(uni_fit[["sm"]][["parameters"]][1]), scale =as.numeric(uni_fit[["sm"]][["parameters"]][2]), shape =as.numeric(uni_fit[["sm"]][["parameters"]][3]))
x2_pseudo=qgev(as.matrix(p),loc=as.numeric(uni_fit[["sf"]][["parameters"]][1]), scale =as.numeric(uni_fit[["sf"]][["parameters"]][2]), shape = as.numeric(uni_fit[["sf"]][["parameters"]][3]))
x3_pseudo=qgev(as.matrix(p),loc=as.numeric(uni_fit[["rf"]][["parameters"]][1]), scale =as.numeric(uni_fit[["rf"]][["parameters"]][2]), shape = as.numeric(uni_fit[["rf"]][["parameters"]][3]))
U1_sim=pgev(as.matrix(x1_pseudo),loc=as.numeric(uni_fit[["sm"]][["parameters"]][1]), scale =as.numeric(uni_fit[["sm"]][["parameters"]][2]), shape =as.numeric(uni_fit[["sm"]][["parameters"]][3]))
U2_sim=pgev(as.matrix(x2_pseudo),loc=as.numeric(uni_fit[["sf"]][["parameters"]][1]), scale =as.numeric(uni_fit[["sf"]][["parameters"]][2]), shape = as.numeric(uni_fit[["sf"]][["parameters"]][3]))
U3_sim=pgev(as.matrix(x3_pseudo),loc=as.numeric(uni_fit[["rf"]][["parameters"]][1]), scale =as.numeric(uni_fit[["rf"]][["parameters"]][2]), shape = as.numeric(uni_fit[["rf"]][["parameters"]][3]))
U_sim=cbind(U1_sim,U2_sim,U3_sim)

data_sim=cbind(x1_pseudo,x2_pseudo)

# kendall correlation
tau_sim=cor(data_sim, method="kendall")
# copula parameter
alpha=iTau(gumbelCopula(), tau=0.361)
# dimension 2
d=2
# gumbel model
gc=gumbelCopula(alpha, dim=d)
# bivariate distribution function
data_plt=matrix(data=NA,nrow=25000,ncol=3)## vector form
dd=expand.grid(SM=x1_pseudo,SF=x2_pseudo)
FF=as.vector(t(F_x1x2_temp))
dd["CDF"]=FF
F_x1x2_temp=matrix(data=NA, nrow=500,ncol=500)
x1_temp=matrix(data=NA, nrow=500, ncol=500) ## RF and soil moisture 
x2_temp=matrix(data=NA, nrow=500, ncol=500) 
for (i in 1:500 ){
  u1_temp=U1_sim[i,1]
  for (j in 1:500){
    F_x1x2_temp[i,j]=pCopula(cbind(u1_temp,U2_sim[j,1]),gc)
    x1_temp[i,j]=U1_sim[i,1]
    x2_temp[i,j]=U2_sim[j,1]
    
  }
  
}

write.csv(F_x1x2_temp, "F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/simulated_bivarite/bi_cdf.csv")
write.csv(x1_temp,"F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/simulated_bivarite/x1_cdf.csv")
write.csv(x2_temp, "F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/simulated_bivarite/x2_cdf.csv")



## plotting surface plot

plt_u=c(0.2,0.4,0.6,0.8,0.9995)
x1_plt=qgev(plt_u,loc=as.numeric(uni_fit[["sm"]][["parameters"]][1]), scale =as.numeric(uni_fit[["sm"]][["parameters"]][2]), shape =as.numeric(uni_fit[["sm"]][["parameters"]][3]))
x2_plt=qgev(plt_u,loc=as.numeric(uni_fit[["sf"]][["parameters"]][1]), scale =as.numeric(uni_fit[["sf"]][["parameters"]][2]), shape = as.numeric(uni_fit[["sf"]][["parameters"]][3]))
x3_plt=qgev(plt_u,loc=as.numeric(uni_fit[["rf"]][["parameters"]][1]), scale =as.numeric(uni_fit[["rf"]][["parameters"]][2]), shape = as.numeric(uni_fit[["rf"]][["parameters"]][3]))

## plotting 
#axx=list(ticktext = list("51.3", "61.9", "72.5", "86.6", "165.5"),tickvals = list(0.2, 0.4, 0.6, 0.8, 1.0),title="Rainfall(mm)")
axx=list(ticktext = list("0.41", "0.42", "0.429", "0.44", "0.473"),tickvals = list(0.2, 0.4, 0.6, 0.8, 1.0),title="Soil Moisture(m/m3")
axy=list(ticktext = list("1312", "1896", "2799", "4909", "337400"),tickvals = list(0.2, 0.4, 0.6, 0.8, 1.0), title="Streamflow(Cumec)")
axz=list(title="CDF")


## observed data
fig <- plot_ly(data = da, x = ~sm, y = ~sf,type="scatter", marker=list(size=5,color="red"))
fig

fig <- plot_ly(x=x1_temp,y=x2_temp, z =F_x1x2_temp , type="surface", contours= list(showlabels=TRUE))
#fig=add_trace(fig,type='scatter',mode="line",x=da$sm,y=da$sf)
fig=fig%>% layout(scene=list(xaxis=axx, yaxis=axy, zaxis=axz))
fig


fig <- plot_ly(data=dd, x=~SM, y=~SF,z=~CDF , type="contour", contours= list(showlabels=TRUE))
#fig=add_trace(fig,data=da, x=da$sm,y=da$sf, type='scatter',mode="line", marker=list(size=5, color="red",symbol=104))
fig=fig%>% layout(scene=list(xaxis=axx, yaxis=axy, zaxis=axz))
fig



x=x1_temp,y=x2_temp, 

%>% plot_ly(da, x=da$sm,y=da$sf,type = "scatter")

  
  
plot_ly(F_x1x2_temp)
plot3D(x = d$x1_pseudo, y = x2_pseudo, z = F_x1x2_temp, colkey=TRUE, bty="b2",
       phi = 40, theta = 30, main="Half of a Torus")
