# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 10:35:36 2022

@author: bhargav
"""

#%%
import pandas as pd
import netCDF4 as nc
import xarray as xr
import numpy as np
import seaborn as snp
import matplotlib.pyplot as plt
from fitter import Fitter, get_common_distributions, get_distributions
import math
import os
os.chdir("F:/Research Work/Synchronization work/Multivariate IDF/Py_codes")
from xarray_to_df import xr_to_df
from latlon_points import latlon_to_points
from circular_statistics import c_stats
import pymannkendall as mk

#%%
# ### don't run
# rf=xr.open_dataarray("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/data/rr_ens_mean_0.25deg_reg_v25.0e.nc")
# temp=xr.open_dataarray("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/data/tg_ens_mean_0.25deg_reg_v25.0e.nc")
# sm=xr.open_dataarray("F:/Research Work/Synchronization work/Multivariate IDF/soil moisture/1984_2020sm.nc")
# lat=rf.latitude.values
# lon=rf.longitude.values

# grids_pts=latlon_to_points(lat, lon)

# # required points
# points=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Data/station details/nearest_points.xlsx", index_col=0)
# points=points.drop("S_No", axis=1)
# lat=points["lat"].tolist()
# lon=points["lon"].tolist()

# data_rf=xr_to_df(rf, lat, lon)
# data_temp=xr_to_df(temp, lat, lon)
# data_sm=xr_to_df(sm, lat, lon)

# rf_f=data_rf.mean(axis=1)
# temp_f=data_temp.mean(axis=1)
# sm_f=data_sm.mean(axis=1)
# rf_f.to_excel("")


#%%
#Input data
dis_cbridge=pd.read_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/data/northatlantic streamflowdata/sf_6608300.xlsx', sheet_name='sf',index_col=0, parse_dates=True)
dis_cbridge=dis_cbridge.rename(columns={"runoff_mean":"sf"}) # changing variable name
surge_newport=pd.read_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/data/storm surge_Reconstructed_north altanta/excel_files/newport_p057_uk.xlsx', sheet_name='surge', index_col=0, parse_dates=(True))
surge_newport=surge_newport.rename(columns={"surge_reconsturcted":"surge"}) # changing variable name
rf=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Data/Data_at_cbrdge/rainfall_cbridge.xlsx", index_col=0, parse_dates=(True))
sm=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Data/Data_at_cbrdge/soil_moisture_cbridge.xlsx", index_col=0, parse_dates=(True))
temp=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Data/Data_at_cbrdge/temp_cbridge.xlsx", index_col=0, parse_dates=(True))
##merging of data
data=dis_cbridge.join(surge_newport,how='inner')

## subset of extremes
data_sub=data[(data["surge"]>=0.8) & (data["sf"]>=300)]

data_sub2=data[(data["surge"]>=0.8)]


## plotting data
ax=data.plot(kind="scatter", x="surge", y="sf", 
             s=10,color="DarkBlue")
data_sub.plot(kind="scatter", x="surge",y="sf",
              s=10, color="red", ax=ax)
ax.axhline(y=296, color="black")
ax.axvline(x=0.8, color="black")
plt.xlabel("Daily maximum storm surge"  r'$[m]$')
plt.ylabel("Daily Streamflow" r'$[m^3/s]$')


# conditional plot

ax2=data.plot(kind="scatter", x="surge", y="sf", 
             s=10,color="DarkBlue")
data_sub2.plot(kind="scatter", x="surge",y="sf",
              s=10, color="red", ax=ax2)
ax2.axvline(x=0.78, color="black")
plt.xlabel("Daily maximum storm surge"  r'$[m]$')
plt.ylabel("Daily Streamflow" r'$[m^3/s]$')





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


## months of occurrence
# ams_mn=pd.concat([sf_circular["month"],rf_circular["month"],surge_circular["month"]], axis=1)



# ams_m=pd.concat([sf_circular["mn"],rf_circular["mn"],surge_circular["mn"]], axis=1)
# ams_m.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Multivarite Analysis/months_cbridge.xlsx")
# ams_m.plot()

# ams_mn.plot()
# df=ams_mn["sf"].apply(lambda x:month_lables[x])

# ams_m_wy=ams_m.replace([1,2,3,4,5,6,7,8,9,10,11,12],[5,6,7,8,9,10,11,12,1,2,3,4])

# ams_m_wy.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Multivarite Analysis/months_wy_cbridge.xlsx")


# ams["sf"].join(ams_w_wy["sf"])
# ams_m_smooth=ams_m.rolling(3).mean()
ax=ams_m_smooth.plot()
### plotting the annual maxima series
# fig, ax=plt.subplots()
# ax.plot(ams_m.index.values, ams_m["sf"],'bo' ,lw=2.0, label="SF", linestyle="solid")
# ax.plot(ams_m.index.values, ams_m["rf"],'rD', lw=2.0, label="RF",linestyle="solid")
# ax.plot(ams_m.index.values, ams_m["surge"],'g^',lw=2.0, label="Surge",linestyle="solid")
# ax.legend()




#%%

# dis=xr.open_dataset("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/GRDC-Daily.nc")
# dis_river=dis.river_name
# dis_r=dis_river.to_dataframe()
# dis_s=dis.station_name.to_dataframe()
# dis_f=dis.to_dataframe()
# # multilevel index to time index

# dis_nao=dis_f.reset_index()


# #dis_nao.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_data.xlsx")
# dis_ss=list(dis_nao['id'].unique())

# for i in range(23):
#     globals()[f'sf_{dis_ss[i]}']=dis_nao.loc[(dis_nao['id']== dis_ss[i])] # []- output is dataframe
    
# # writing the data into excel files
# # sf_6608100.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608100.xlsx') 
# # sf_6608101.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608101.xlsx') 
# # sf_6608110.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608110.xlsx') 
# # sf_6608150.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608150.xlsx') 
# # sf_6608160.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608160.xlsx') 
# # sf_6608170.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608170.xlsx') 
# # sf_6608190.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608190.xlsx') 
# # sf_6608200.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608200.xlsx') 
# # sf_6608210.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608210.xlsx') 
# # sf_6608220.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608220.xlsx') 
# # sf_6608230.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608230.xlsx') 
# # sf_6608235.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608235.xlsx') 
# # sf_6608240.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608240.xlsx') 
# # sf_6608250.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608250.xlsx') 
# # sf_6608300.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608300.xlsx') 
# # sf_6608310.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608310.xlsx') 
# # sf_6608500.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608500.xlsx') 
# # sf_6608501.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608501.xlsx') 
# # sf_6608502.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608502.xlsx') 
# # sf_6608510.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608510.xlsx') 
# # sf_6608520.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608520.xlsx') 
# # sf_6608530.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608530.xlsx') 
# # sf_6608600.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608600.xlsx') 

#%%

# # importing the reconstructed storm surge date at newport station
# dis_cbridge=pd.read_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/northatlantic streamflowdata/sf_6608300.xlsx', sheet_name='sf',index_col=0, parse_dates=True)
# surge_newport=pd.read_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/storm surge_Reconstructed_north altanta/excel_files/newport_p057_uk.xlsx', sheet_name='surge', index_col=0, parse_dates=(True))

# dis_cbride_surge_newport=surge_newport.join(dis_cbridge, how='inner')
# dis_cbride_surge_newport.corr(method='kendall')

# dis_cbride_surge_newport.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/dis_cbridge_surge_newport.xlsx')
# # extracting the annual maxima values
# ams=dis_cbride_surge_newport.groupby(dis_cbride_surge_newport.index.year).max()
# ams.corr(method='kendall')

# # extraacting the date of occurence

# ams_dis_cbridge_day=dis_cbride_surge_newport["runoff_mean"].groupby(dis_cbride_surge_newport.index.year).idxmax(axis=0)

# ams_surge_newport_day=dis_cbride_surge_newport["surge_reconsturcted"].groupby(dis_cbride_surge_newport.index.year).idxmax(axis=0)
# # stroing the data into excel files

# ams.to_excel('F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/ams_dis_cbridge_surge_newport.xlsx')
# ams_surge_newport_day.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/ams_surge_newport_day.xlsx")
# ams_dis_cbridge_day.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/ams_dis_cbridge_day.xlsx")


#%%

#streamflow circular statistics

ams=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Multivarite Analysis/ams_dis_cbridge_surge_newport.xlsx", index_col=0, parse_dates=True)

sf_day_ams=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Multivarite Analysis/AMS/ams_dis_cbridge_day.xlsx", index_col=0, parse_dates=(True))

# year series
year=np.array(ams.index.year)
sf_day=sf_day_ams.squeeze()

# extractig day of the year from date of occurence 
day_sf=sf_day.dt.dayofyear
#series to pandas dataframe
ams_dis=pd.DataFrame()
ams_dis["day"]=pd.DataFrame(data=day_sf)
ams_dis["sf"]=ams["runoff_mean"]
ams_dis["year"]=year
# use circular statistics function

#[x_bar,y_bar,ma,md,r_bar,sigma_sq]= c_stats(ams_dis)
season_dis=c_stats(ams_dis)
    
#%%
#storm surge circular statistics
surge_day_ams=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Multivarite Analysis/AMS/ams_surge_newport_day.xlsx", index_col=0, parse_dates=(True))


# extracting day of the year from date of occurence
storm_day=surge_day_ams.squeeze()
day_storm=storm_day.dt.dayofyear

# input dataframe 
ams_surge=pd.DataFrame() 
ams_surge["day"]=pd.DataFrame(data=day_storm) 
ams_surge["surge"]=ams["surge_reconsturcted"]           
ams_surge["year"]=year             

#[x_bar,y_bar,ma,md,r_bar,sigma_sq] = c_stats(ams_surge)

season_surge=c_stats(ams_surge)


#%%%
#nc.file, input_points
rf=xr.open_dataarray("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/data/rr_ens_mean_0.25deg_reg_v25.0e.nc")
temp=xr.open_dataarray("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/data/tg_ens_mean_0.25deg_reg_v25.0e.nc")
sm=xr.open_dataarray("F:/Research Work/Synchronization work/Multivariate IDF/soil moisture/1984_2020sm.nc")
lat=rf.latitude.values
lon=rf.longitude.values

grids_pts=latlon_to_points(lat, lon)

# required points
points=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Data/station details/nearest_points.xlsx", index_col=0)
points=points.drop("S_No", axis=1)
lat=points["lat"].tolist()
lon=points["lon"].tolist()

data_rf=xr_to_df(rf, lat, lon)
data_temp=xr_to_df(temp, lat, lon)
data_sm=xr_to_df(sm, lat, lon)

rf_f=data_rf.mean(axis=1)
temp_f=data_temp.mean(axis=1)
sm_f=data_sm.mean(axis=1)


# annual maxima series 
# rainfall
ams_rf=rf_f.groupby(rf_f.index.year).max()
ams_rf_day=rf_f.groupby(rf_f.index.year).idxmax()
rf_day=ams_rf_day.squeeze()
day_rf=rf_day.dt.dayofyear

# input dataframe 
rf_ams=pd.DataFrame() 
rf_ams["day"]=pd.DataFrame(data=day_rf) 
rf_ams["rf"]=ams_rf          
rf_ams["year"]=rf_day.dt.year 
rf_ams=rf_ams.set_index(ams_rf_day)
rf_stat=c_stats(rf_ams)

# soil moisture
ams_sm=sm_f.groupby(sm_f.index.year).max()
ams_sm_day=sm_f.groupby(sm_f.index.year).idxmax()
sm_day=ams_sm_day.squeeze()
day_sm=sm_day.dt.dayofyear

# input dataframe 
sm_ams=pd.DataFrame() 
sm_ams["day"]=pd.DataFrame(data=day_sm) 
sm_ams["rf"]=ams_sm          
sm_ams["year"]=sm_day.dt.year 
sm_ams=sm_ams.set_index(ams_sm_day)
sm_stat=c_stats(sm_ams)




#merging the seasonality statistics 
#seasonal_stats=pd.merge(season_dis,season_surge, how="left")
#                        ,sm_stat,rf_stat)



#%%
ams_temp=temp_f.groupby(temp_f.index.year).max()
ams_sm=sm_f.groupby(sm_f.index.year).max()


# first extarct days also 
# then look the sf and follow the same procedure.

## circular statistics
rf_stat=c_stats(rf_f)
temp_stat=c_stats(temp_f)
sm_stat=c_stats(sm_f)


#%%
# writing the data excel files
rf_f.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Data/rainfall_cbridge.xlsx")
temp_f.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Data/temp_cbridge.xlsx")
sm_f.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Data/soil_moisture_cbridge.xlsx")
