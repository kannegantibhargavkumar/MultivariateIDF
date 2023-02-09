# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 09:54:47 2022

@author: bhargav
"""

#%%
# reading multiple text files of IMD data
import pandas as pd
import os
import glob
import pylab as py
import seaborn as sns 
import numpy as np
from pylab import plot, show, hist, figure,title
import matplotlib.pyplot as plt
import netCDF4 as nc
import xarray as xr
#%%
fil="F:/Research Work/Synchronization work/Multivariate IDF/PCP_IMD_2019g/rain_1951_2020.nc"
ds=xr.open_dataset(fil)
# choosing the desired grid points
ds_f=ds.rain.where((ds.rain.lat>=9) & (ds.rain.lat<=10.25) & (ds.rain.lon>=76) & (ds.rain.lon<=77.5), drop=True) 

pcp_p=ds_f.sel(time=slice('1985-01-01','2018-12-31'))

rf_pcp=pcp_p.to_dataframe()

rf_periyar=rf_pcp.reset_index(level=[1,2]).drop(['spatial_ref'], axis=1) # final rainfall dataset

## sepearte the time series at each grid point
#lat=rf_periyar.lat.unique()
#lon=rf_periyar.lon.unique()
time=rf_periyar.index.unique()


# importing the grid points using excel file
grid_pts=pd.read_excel('F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_pcp/pcp_points_periyar.xlsx', index_col=0, parse_dates=True)

## rearranging the data 
rf=pd.DataFrame(index=time, columns=range(9))
rf.index=time

## extracting all the grid points of data
for i in range(9):
    col_id=[rf.columns[i]]
    a=grid_pts.lat[i]
    b=grid_pts.lon[i]
    pcp=rf_periyar.loc[(rf_periyar.lat==a) & (rf_periyar.lon==b)]
    rf[col_id]=pcp.drop(['lat','lon'], axis=1)



## extracting 2018 year data only

rf.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Kerala Floods Data_to Bhargav/Rainfall/IMD_rf_1985_2018.xlsx")

rf_95=rf.quantile(0.95, axis=0)
rf_99=rf.quantile(0.99, axis=0)
rf_995=rf.quantile(0.995, axis=0)

# example of grid 8
t=pd.date_range('2018-06-01','2018-08-17')
rf_2018_monsoon=rf.loc[t]

# grid 7

rf_2018_monsoon[7].plot()



threshold=42.0899  # 95 % percentile at grid 7
values=rf_2018_monsoon[7]
x=range(len(values))
# split it up
above_threshold=np.maximum(values-threshold,0)
below_threshold=np.minimum(values,threshold)
# and plotting it
fig, ax=plt.subplots()
ax.bar(x, below_threshold, 0.35, color="g")
ax.bar(x, above_threshold, 0.35, color="r", bottom=below_threshold)
## horizontal line indicating the threshold
ax.plot([], [])





## annual maximas 1d, 3d, 4d, 5d
ams_pcp_1day=rf.groupby(rf.index.year).max()
ams_pcp_1day_rday=rf.groupby(rf.index.year).idxmax(axis=0)

# 3day maximum
rf_3d=rf.rolling(3).sum()

ams_pcp_3day=rf_3d.groupby(rf_3d.index.year).max()
ams_pcp_3day_rday=rf_3d.groupby(rf_3d.index.year).idxmax(axis=0)

# 4day maxium
rf_4d=rf.rolling(4).sum()
ams_pcp_4day=rf_4d.groupby(rf_4d.index.year).max()
ams_pcp_4day_rday=rf_4d.groupby(rf_4d.index.year).idxmax(axis=0)

# 5day maxium
rf_5d=rf.rolling(5).sum()

ams_pcp_5day=rf_5d.groupby(rf_5d.index.year).max()
ams_pcp_5day_rday=rf_5d.groupby(rf_5d.index.year).idxmax(axis=0)

## writing into excel files

#1 day excel file
ams_pcp_1day.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_pcp/New/ams_pcp_1day.xlsx")
ams_pcp_1day_rday.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_pcp/New/ams_pcp_1day_rday.xlsx")
# 3 day excel file

rf_3d.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_pcp/New/rf_3day.xlsx")

ams_pcp_3day.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_pcp/New/ams_pcp_3day.xlsx")
ams_pcp_3day_rday.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_pcp/New/ams_pcp_3day_rday.xlsx")
# 5 day excel file
rf_5d.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_pcp/New/rf_5day.xlsx")

ams_pcp_5day.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_pcp/New/ams_pcp_5day.xlsx"rday)
ams_pcp_5day_rday.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_pcp/New/ams_pcp_5day_rday.xlsx")
## correlation

# run the discharge code file
ams_dis_1day=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_dis_new/ams_dis_1day.xlsx", index_col=0,parse_dates=True )

ams_1day_pcp_dis=ams_dis_1day.join(ams_pcp_1day, how='inner')
ams_1day_pcp_dis=ams_1day_pcp_dis.drop('day', axis=1)

cor_1day=ams_1day_pcp_dis.corr(method='kendall')

## 3DAY 
ams_dis_3day=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_dis_new/ams_dis_3day.xlsx", index_col=0,parse_dates=True )
ams_dis_3day.index=ams_pcp_3day.index

ams_3day_pcp_dis=ams_dis_3day.join(ams_pcp_3day, how='inner')
ams_3day_pcp_dis=ams_3day_pcp_dis.drop('day', axis=1)

cor_3day=ams_3day_pcp_dis.corr(method='kendall')

## 5day
ams_dis_5day=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_dis_new/ams_dis_5day.xlsx", index_col=0,parse_dates=True )
ams_dis_5day.index=ams_pcp_5day.index

ams_5day_pcp_dis=ams_dis_5day.join(ams_pcp_5day, how='inner')
ams_5day_pcp_dis=ams_5day_pcp_dis.drop('day', axis=1)

cor_5day=ams_5day_pcp_dis.corr(method='kendall')

# writing the correlation to excel files
cor_1day.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/Kendall_correlation/cor_1day.xlsx")
cor_3day.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/Kendall_correlation/cor_3day.xlsx")
cor_5day.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/Kendall_correlation/cor_5day.xlsx")
