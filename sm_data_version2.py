# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 21:15:21 2022

@author: bhargav
"""
#%%
import xarray as xr
import pandas as pd
import numpy as np
import math
import os
os.chdir("F:/Research Work/Synchronization work/Multivariate IDF/Py_codes")
from circular_statistics import c_stats

#%%

sm=xr.open_dataset("F:/Research Work/Synchronization work/soil moisture data/1984_2020sm_to_IMD_grid/1984_2020sm_to_IMD_grid.nc")

# renaming the data variable
sm['sm'] = sm['__xarray_dataarray_variable__']

# dropping the duplicate variable
sm= sm.drop(['__xarray_dataarray_variable__'])

temp=sm.sm.where((sm.sm.lat==17.25) & (sm.sm.lon==74) , drop=True)

t=temp.to_dataframe()
# reindexing the dataframe 


sm_warunji=t.reset_index(level=[1,2]).drop(['lat','lon'], axis=1)
# daily soil moisture data at warunji
sm_warunji.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/sm_warunji.xlsx")
##
#Extracting the annaul maxima series
ams_sm=sm_warunji.groupby(sm_warunji.index.year).max()
sm_ams_day=sm_warunji.groupby(sm_warunji.index.year).idxmax()

#%%
# calculating the circular statistics
day=sm_ams_day.squeeze()
dy_sm=day.dt.dayofyear
yr=ams_sm.index.values
# input data for circular statistics function
ams=pd.DataFrame()
ams["day"]=dy_sm
ams["sm"]=ams_sm
ams["year"]=yr

#%%
# using function here'
# file is saved as circular statistics
[x_bar,y_bar,ma,md,r_bar,sigma_sq] = c_stats(ams)     
    

#%%
#exporting the soil moisture data
sm_warunji.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/soilmoisture.xlsx")
ams.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/ams_soilmoisture.xlsx")
