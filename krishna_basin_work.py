# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 12:28:05 2022

@author: bhargav
"""

#%%
import pandas as pd 
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import seaborn as sns

from distfit import distfit
import math
import ternary
import os
os.chdir("F:/Research Work/Synchronization work/Multivariate IDF/Py_codes")
from xarray_to_df import xr_to_df
from latlon_points import latlon_to_points
from circular_statistics import c_stats
import pymannkendall as mk

#import plotly.express as px
#%%
## Rainfall data 
# finding the nearest upstream precipitation stations

lat_warunji= 17.272
lon_warunji=74.165

file="F:/Research Work/Synchronization work/Multivariate IDF/PCP_IMD_2019g/rain_1951_2020.nc"
ds=xr.open_dataset(file)

ds_f=ds.rain.where((ds.rain.lat==17.25) & (ds.rain.lon==74) , drop=True) 

pcp_p=ds_f.sel(time=slice('1965-01-01','2018-12-31'))

tr=pcp_p.to_dataframe()

# rainfall grid near to warunji stream flow station

rf_warunji=tr.reset_index(level=[1,2]).drop(['spatial_ref'], axis=1) 
#rf_warunji.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/data_input/rf_warunji.xlsx")
#%%
#pre-processing the data (modified on 23/06/2022)
sf_warunji=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/data_input/sf_Data.xlsx",index_col=0, parse_dates=True)
sf_warunji=sf_warunji.rename(columns={"Discharge":"sf"}) # changing variable name
rf_warunji=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/data_input/rf_warunji.xlsx", sheet_name="rf", index_col=0, parse_dates=(True))
rf_warunji=rf_warunji.rename(columns={"rain":"rf"}) # changing variable name
sm_warunji=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/data_input/sm_warunji.xlsx", index_col=0,parse_dates=(True))

#merging of data
data=sf_warunji.join(rf_warunji,how='inner').join(sm_warunji,how="inner")
# ## optional
# ## working with water year (i.e. Hydrologic year)
# data["water_year"]=data.index.to_series().dt.year.where(data.index.to_series().dt.month<6, data.index.to_series().dt.year+1)

# annual_maxima_wateryear=data.groupby('water_year').max()
# annual_maxima_wateryear_dates=data.groupby('water_year').idxmax()


## extracting annual maxima values 
ams=data.groupby(data.index.year).max()
ams.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/data_input/ams_warunji.xlsx")

## date of occurence of annual maxima
ams_date=data.groupby(data.index.year).idxmax()
ams_date.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/data_input/dates_ams.xlsx")


# extractig day of the year from date of occurence 
# day=ams_day.dt.dayofyear
month_lables={1:'J', 2:'F', 3:'M',4:'A', 5:'M',6:'J',7:'J',8:'A',9:'S',10:'O',11:'N',12:'D'}
## circular statistics without soil moisture
for i in range(ams.shape[1]):
    temp=data.iloc[:,i]
    temp=pd.DataFrame(temp)
    k=''.join(list(temp.columns))
    s= c_stats(temp)
    dy=ams_date.iloc[:,i].dt.dayofyear
    mn=ams_date.iloc[:,i].dt.month
    mn_str=mn.apply(lambda x:month_lables[x])
    globals()[k+ "_circular"]={"stats":s, "day":dy, "month":mn_str, "mn":mn}


#day of occurrence of ams of all variables
d=pd.concat([sf_circular["day"],rf_circular["day"],sm_circular["day"]], axis=1)
mn_warunji=pd.concat([sf_circular["mn"],rf_circular["mn"],sm_circular["mn"]], axis=1)

mn_warunji.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/data_input/months_warunji.xlsx")


## annual maxima values and dates
ams=data.groupby(data.index.year).max()
ams.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Multivarite Analysis/ams_cbridge.xlsx")


#%%
# data rearrangement
# extracting AMS and date of occurence

sf_warunji=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/data_input/sf_Data.xlsx",index_col=0, parse_dates=True)
sf_warunji=sf_warunji.rename(columns={"Discharge":"sf"}) # changing variable name

ams_warunji=sf_warunji.groupby(sf_warunji.index.year).max()
# date of occurence of AMS
ams_sf_rday=sf_warunji.groupby(sf_warunji.index.year).idxmax(axis=0)

ams_warunji.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/ams_warunji.xlsx")
ams_sf_rday.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/ams_sfday_warunji.xlsx")
# extracting the day of occurence of AMS 
day_sf=ams_sf_rday["Discharge"].dt.dayofyear



ams_rf_w=rf_warunji.groupby(rf_warunji.index.year).max()
# date of occurence of AMS
ams_rf_rday=rf_warunji.groupby(rf_warunji.index.year).idxmax(axis=0)
# extracting day of the year 
day_rf=ams_rf_rday["rain"].dt.dayofyear
ams_rf_rday.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/ams_rfday_warunji.xlsx")

# combine streamflow data pcp data near to warunji
# annual maximum data
rf_sf_warunji= ams_warunji.join(ams_rf_w['rain'], how='inner')
corr_rf_sf_warunji=rf_sf_warunji.corr(method='kendall')

rf_sf_warunji.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/ams_rf_sf_warunji.xlsx")

ams_rf_rday.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/ams_sfday_warunji.xlsx")
# daily data merge
rf_sf_d_warunji=sf_warunji.join(rf_warunji['rain'])
corr_daily=rf_sf_d_warunji.corr(method='kendall')
## fitting marginal distribution to return levels

rt_dis= pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/dis_rt.xlsx", sheet_name='rt', index_col=0, parse_dates=(True))




#%%
#Input data
dis_cbridge=pd.read_excel('F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/sf_Data.xlsx',index_col=0, parse_dates=True)
dis_cbridge=dis_cbridge.rename(columns={"runoff_mean":"sf"}) # changing variable name
surge_newport=pd.read_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/data/storm surge_Reconstructed_north altanta/excel_files/newport_p057_uk.xlsx', sheet_name='surge', index_col=0, parse_dates=(True))
surge_newport=surge_newport.rename(columns={"surge_reconsturcted":"surge"}) # changing variable name
rf=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Data/Data_at_cbrdge/rainfall_cbridge.xlsx", index_col=0, parse_dates=(True))
sm=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Data/Data_at_cbrdge/soil_moisture_cbridge.xlsx", index_col=0, parse_dates=(True))
temp=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Data/Data_at_cbrdge/temp_cbridge.xlsx", index_col=0, parse_dates=(True))
##merging of data
data=dis_cbridge.join(surge_newport,how='inner').join(rf,how="inner")
data_sm=dis_cbridge.join(surge_newport,how='inner').join(rf,how="inner").join(sm,how="inner")
## working with water year (i.e. Hydrologic year)
data_sm["water_year"]=data_sm.index.to_series().dt.year.where(data_sm.index.to_series().dt.month<10, data_sm.index.to_series().dt.year+1)

annual_maxima_wateryear=data_sm.groupby('water_year').max()
annual_maxima_wateryear_dates=data_sm.groupby('water_year').idxmax()

## annual maxima values and dates
ams=data.groupby(data.index.year).max()
ams.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Multivarite Analysis/ams_cbridge.xlsx")

ams_sm=data_sm.groupby(data_sm.index.year).max()

## date of occurence of annual maxima
ams_date=data.groupby(data.index.year).idxmax()

ams_sm_date=data_sm.groupby(data_sm.index.year).idxmax()

##day of occurence
# year series
# year=np.array(ams.index.year)
# ams_day=ams_date.squeeze()

# extractig day of the year from date of occurence 
# day=ams_day.dt.dayofyear
month_lables={1:'J', 2:'F', 3:'M',4:'A', 5:'M',6:'J',7:'J',8:'A',9:'S',10:'O',11:'N',12:'D'}
## circular statistics without soil moisture
for i in range(ams.shape[1]):
    temp=data.iloc[:,i]
    temp=pd.DataFrame(temp)
    k=''.join(list(temp.columns))
    s= c_stats(temp)
    dy=ams_date.iloc[:,i].dt.dayofyear
    mn=ams_date.iloc[:,i].dt.month
    mn_str=mn.apply(lambda x:month_lables[x])
    globals()[k+ "_circular"]={"stats":s, "day":dy, "month":mn_str, "mn":mn}


#day of occurrence of ams of all variables
d=pd.concat([sf_circular["day"],rf_circular["day"],surge_circular["day"]], axis=1)
d_smooth=d.rolling(3).mean()
d_smooth.plot()
d.plot()

d.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Multivarite Analysis/dates_cbridge.xlsx")

###along with soil moisture

for i in range(ams_sm.shape[1]):
    temp_sm=data_sm.iloc[:,i]
    temp_sm=pd.DataFrame(temp_sm)
    k=''.join(list(temp_sm.columns))
    s_sm= c_stats(temp_sm)
    dy_sm=ams_sm_date.iloc[:,i].dt.dayofyear
    mn_sm=ams_sm_date.iloc[:,i].dt.month
    globals()[k+ "_sm_circular"]={"stats":s_sm, "day":dy_sm, "month":mn_sm}
    
d_s=pd.concat([sf_sm_circular["day"],rf_sm_circular["day"],surge_sm_circular["day"], sm_sm_circular["day"]], axis=1)
d_s.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Multivarite Analysis/with SM/dates_cbridge.xlsx")
d_smooth_s=d_s.rolling(3).mean()
d_smooth_s.plot()
d_s.plot()
## month series
ams_sm_mn=pd.concat([sf_sm_circular["month"],rf_sm_circular["month"],surge_sm_circular["month"], sm_sm_circular["month"]], axis=1)

ams_sm_m_wy=ams_sm_mn.replace([1,2,3,4,5,6,7,8,9,10,11,12],[5,6,7,8,9,10,11,12,1,2,3,4])

ams_sm_m_wy.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Multivarite Analysis/with SM/months_modified_cbridge.xlsx")

ams_sm.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Multivarite Analysis/with SM/ams_cbridge.xlsx")

#%%
#writing function for circular statistics


def c_stats(data):
    # This function calculates the circular statistics such as 
    # Mean date of occurence, peristance and variablity in date of 
    #occurence of extremes
    # input: (i) ams dataframe contain two columns
    # dayofyear, variable and year series 
    #data["year"]=data.index.values
    temp=data.to_numpy()
    dy=temp[:,0]
    j=temp[:,2] # year 
    sf=temp[:,1] # variable
    x=np.zeros(len(dy))
    y=np.zeros(len(dy))
    
    for i in range(len(dy)):
         # if leap year 
        if((j[i]%400==0)or(j[i]%100 != 0)and(j[i]%4 ==0)):
          len_yr=366
          theta= (dy[i]*(2*180))/len_yr
          x[i]=sf[i]*math.cos(theta*math.pi/180)
          y[i]=sf[i]*math.sin(theta*math.pi/180)
        else:
          len_yr=365
          theta= (dy[i]*(2*180))/len_yr
          x[i]=sf[i]*math.cos(theta*math.pi/180)
          y[i]=sf[i]*math.sin(theta*math.pi/180)
    
    # coordinates of extremes occurences
    
    x_bar=sum(x)/sum(sf)
    y_bar=sum(y)/sum(sf)
    
    # intializing the value
    ma=[]
    # converting into polar coordinates
    if (x_bar>=0 and y_bar>=0):
        ma=math.degrees(math.atan(y_bar/x_bar));
    elif (x_bar<0 and y_bar>0):
        ma=180-abs(math.degrees(math.atan(y_bar/x_bar)));
    elif (x_bar<0 and y_bar<0):
        ma=180+abs(math.degrees(math.atan(y_bar/x_bar)));
    elif (x_bar>=0 and y_bar<=0):
        ma=360-abs(math.degrees(math.atan(y_bar/x_bar)));
        
    #  # ma= Mean Angle , md= Mean date of occurence   
    md=ma*(365/(2*180))
    # r= persistance, sigma_sq= Variance
    r_bar=np.power((np.square(x_bar)+np.square(y_bar)),0.5)
    sigma_sq=-2*math.log(r_bar)
    
    return x_bar,y_bar,ma,md,r_bar,sigma_sq
  
     
     
    


#%%
#computing seasonal statistics
#streamflow
ams=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/ams_rf_sf_warunji.xlsx",index_col=0,parse_dates=(True))
day=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/ams_sfday_warunji.xlsx",index_col=0,parse_dates=(True))

#input data
#sf
dy_sf=day["Discharge"].squeeze()
day_occ=dy_sf.dt.dayofyear
ams_sf=pd.DataFrame()
ams_sf["day"]=day_occ
ams_sf["sf"]=ams["Discharge"]
ams_sf["year"]=ams.index.year
ams_sf.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/ams_streamflow.xlsx")

# using function here
[x_bar,y_bar,ma,md,r_bar,sigma_sq] = c_stats(ams_sf)

# rainfall

dy_rf=day["RF"].squeeze()
day_occ_rf=dy_rf.dt.dayofyear
ams_rf=pd.DataFrame()
ams_rf["day"]=day_occ_rf
ams_rf["rf"]=ams["rain"]
ams_rf["year"]=ams.index.year
ams_rf.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/ams_rainfall.xlsx")
# using function here
[x_bar,y_bar,ma,md,r_bar,sigma_sq] = c_stats(ams_rf)



#%%
#Importing ams data 
ams_sm=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/ams_soilmoisture.xlsx", index_col=0, parse_dates=(True))
ams_rf=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/ams_rainfall.xlsx", index_col=0, parse_dates=(True))
ams_sf=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/ams_streamflow.xlsx", index_col=0, parse_dates=(True))
# merging all variables at warunji
ams_warunji=pd.merge(pd.merge(ams_sm,ams_rf, on="year"), ams_sf, on="year")

ams_warunji.index=ams_warunji["year"]
ams_warunji=ams_warunji.drop("year", axis=1)
# seperating data and dates
ams_values=pd.DataFrame(index=ams_warunji.index )
ams_values["sm"]=ams_warunji["sm"]
ams_values["rf"]=ams_warunji["rf"]
ams_values["sf"]=ams_warunji["sf"]

ams_values.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/ams_warunji.xlsx")
# dependence among variables
corr_warunji=ams_values.corr(method="spearman")
corr_warunji_2=ams_values.corr(method="pearson")







#%%

# creating a dataframe with lags 
def df_derived_by_shift(df,lag=0,NON_DER=[]):
    df = df.copy()
    if not lag:
        return df
    cols ={}
    for i in range(1,lag+1):
        for x in list(df.columns):
            if x not in NON_DER:
                if not x in cols:
                    cols[x] = ['{}_{}'.format(x, i)]
                else:
                    cols[x].append('{}_{}'.format(x, i))
    for k,v in cols.items():
        columns = v
        dfn = pd.DataFrame(data=None, columns=columns, index=df.index)    
        i = 1
        for c in columns:
            dfn[c] = df[k].shift(periods=i)
            i+=1
        df = pd.concat([df, dfn], axis=1)
        df=df.reindex(df.index)
    return df

#%%

## lagged correlation 
NON_DER = ['index',]
df_new = df_derived_by_shift(ams_values, 3, NON_DER)
#df_new=df_new.drop(['Discharge','Tide'],axis=1)
df_new = df_new.dropna()
corr_lag=df_new.corr()
# heat map of correlation
colormap = plt.cm.RdBu
plt.figure(figsize=(15,10))
plt.title(u'Cross and Auto Correlation', y=1.05, size=20)

mask = np.zeros_like(df_new.corr(method='kendall')) # default method is pearson
mask[np.triu_indices_from(mask)] = True

svm = sns.heatmap(df_new.corr(method='kendall'), mask=mask, linewidths=0.1,vmax=1.0, 
            square=True, cmap=colormap, linecolor='white', annot=True)



#%%
# The dataframe contains dateofyear of each variable
ams_days=pd.DataFrame(index=ams_warunji.index)
ams_days["sm_dy"]=ams_warunji["day_x"]
ams_days["rf_dy"]=ams_warunji["day_y"]
ams_days["sf_dy"]=ams_warunji["day"]

#%%



# extracting the soil moisture corresponding to AMS of streamflow

dates=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/ams_sfday_warunji.xlsx",index_col=0,parse_dates=(True))
sm_warunji=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/sm_warunji.xlsx", index_col=0, parse_dates=True)
ams_sf=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/ams_streamflow.xlsx", index_col=0, parse_dates=(True))
dates_sf=dates["Discharge"]

dy=pd.DataFrame()
dy.index=dates_sf
d=pd.merge(dy, sm_warunji, left_index=True, right_index=True)


# intializing the empty variables
sm_sf_day_0=pd.DataFrame()
sm_sf_day_1=pd.DataFrame()
sm_sf_day_2=pd.DataFrame()
sm_sf_day_3=pd.DataFrame()

for i in range(len(dy)):
    t=sm_warunji.index.values==dy.index.values[i];
    k=[j for j, x in enumerate(t) if x]
    if len(k)==0:
        continue
    else:
        sm_sf_day_0=sm_sf_day_0.append(sm_warunji.iloc[k]);
        sm_sf_day_1=sm_sf_day_1.append(sm_warunji.iloc[[k[0]-1]]);
        sm_sf_day_2=sm_sf_day_2.append(sm_warunji.iloc[[k[0]-2]]);
        sm_sf_day_3=sm_sf_day_3.append(sm_warunji.iloc[[k[0]-3]]);
        
# correlation between soil moisture on the day of AMS streamflow

ams_s=ams_sf.iloc[19:54,1]
d_2=pd.DataFrame(data=[ams_s.values,sm_sf_day_2.values]).transpose()

du=pd.concat([ams_s, sm_sf_day_2["sm"]], axis=1)

    
  



#%%
# taken from other file
# heat map of correlation
colormap = plt.cm.RdBu
plt.figure(figsize=(15,10))
plt.title(u'Cross and Auto Correlation', y=1.05, size=16)

mask = np.zeros_like(df_new.corr(method='kendall')) # default method is pearson
mask[np.triu_indices_from(mask)] = True

svm = sns.heatmap(df_new.corr(method='kendall'), mask=mask, linewidths=0.1,vmax=1.0, 
            square=True, cmap=colormap, linecolor='white', annot=True)

# plot of correlation line 
x=[-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6]
# these values are extracted from correlation matrix computed eariler
corr=[0.12,0.14,0.15,0.18,0.19,0.21,0.22,0.24,0.23,0.22,0.22,0.22,0.23]

plt.plot(x,corr,linewidth=2, markersize=12)
plt.scatter(1,0.24, color='r')
plt.ylim([0.1,0.26])
plt.xlabel('Time(Days)')
plt.ylabel('Kendall Tau')

#%%
# Plotting the results
data_sim=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/data_sim.xlsx", index_col=0)
cdf_sim=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/cdf_sim.xlsx", index_col=0)
rt_sim=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/rt_sim.xlsx", index_col=0)
cdf_sim_array=cdf_sim.to_numpy()
rt_sim_array=rt_sim.to_numpy()


## contour plot
x_plt=data_sim["sm_sim"].values
y_plt=data_sim["sf_sim"].values

[X, Y] = np.meshgrid(x_plt, y_plt)
  
fig, ax = plt.subplots(1, 1)
  
Z = rt_sim_array
  
# plots filled contour plot
ax.contourf(X, Y, Z, 20, cmap="RdGy")
fig.colorbar(Z) 
ax.set_title('Filled Contour Plot')
ax.set_xlabel('Soil Moisture')
ax.set_ylabel('Streamflow')

plt.show()

# 3d plot
plt.contourf(X, Y, Z, 20, cmap='RdGy')
plt.colorbar(Z);


# plotting histogram
sns.set_style('white')
sns.set_context("paper", font_scale=2)
sns.displot(data=rt_dis, x="u_rt", kind="hist", bins=54, aspect=1.5 )

# data preparation

rt_uni=rt_dis["u_rt"].values

# fitting distributions
dist=distfit()
dist.fit_transform(rt_uni)

print(dist.summary)

dist.plot()
# fitting distribution for bi_return period
rt_bi=rt_dis["b_rt"].values
dist2=distfit(distr="powerlaw")
dist2.fit_transform(rt_bi)

print(dist2.summary)
dist2.plot()

# bivariate return periods distribution parameters
# shape=-1.46316 loc=1.41995 scale=0.630905
# univarite return periods distribution parameters
# shape=-1.45582 loc=1.40233 scale=0.611259





