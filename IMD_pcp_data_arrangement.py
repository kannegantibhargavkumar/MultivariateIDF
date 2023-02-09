# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 11:43:31 2021

@author: Bhargav
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

#%%

# define the location of the directory
file_dir2=(r'F:\Research Work\Synchronization work\Multivariate IDF\PCP_IMD_2019g')
# different file directory
#file_dir=(r'D:\Research Work\Synchronization work\IMD RF PT25 (1901-2017)\MASTER DATA PT25 (1901-2017)\RF_DATA_PT25\pt25_data')
# change the directory
#os.chdir(file_dir)
os.chdir(file_dir2)
# !pip install imdlib
#import imdlib as im



#%%
#filenames=[i for i in glob.glob("*.grd")]
# reading all the file names in the directory
filenames2=[i for i in glob.glob("*.txt")] 
# taking file names from 1985-2017(bcoz of streamflow data)
d=filenames2[84:119]
# merging all years data into time series
frames=[pd.read_csv(i,header=None, sep="\t") for i in d]
result=pd.concat(frames)
# generating a date series for the data
result['date']=pd.date_range(start='01-Jan-1985', periods=len(result), freq='D')

# indexing the dataframe by date series
# convert the column(it's a string) to datetime type
datetime_series=pd.to_datetime(result['date'])
# create datetime index passing the datetime series
datetime_index=pd.DatetimeIndex(datetime_series.values)
pcp_1985_2019=result.set_index(datetime_index)
# we don't need the date column anymore
pcp_1985_2019.drop('date',axis=1, inplace=True)

# writing pcp data into nc file 
fn='D:\Research Work\Synchronization work\Multivariate IDF\PCP_IMD_2019\pcp_1985_2019.nc'
ds=nc.Dataset(fn,'w',format='NETCDF4')
# creating a dimensions
time=ds.createDimension('time', None)
lat=ds.createDimension('lat', None)
long=ds.createDimension('lon', None)
# add NetCDF variable
times=ds.createVariable('time', 'f4', ('time',))
lats=ds.createVariable('lat','f4',('lat',))
lons=ds.createVariable('lon','f4',('lon',))
value=ds.createVariable('Value','f4',('time','lat','lon',))
# assign latitude and longitude values
lats[:]=np.arange(8.25,37.25,0.25)
lons[:]=np.arange(68,97.25,0.25)
#assign NetCDF data values
# importing the IMD grid points
lat_long_ind=pd.read_excel(r'D:\Research Work\Synchronization work\Multivariate IDF\PCP_IMD_2019\seasonality.xlsx')
lat_ind=lat_long_ind['Latitude'].dropna()
long_ind=lat_long_ind['Longitude']
grid_points=np.transpose(np.reshape(np.meshgrid(lat_ind,long_ind), (2,13806))) # square boundary
# actual grid points
imd_points=pd.read_excel(r'D:\Research Work\Synchronization work\Multivariate IDF\PCP_IMD_2019\seasonality_original.xlsx')





for i in range(12783):
    



# extract grid points of interest
cols=[40,53,64,65,75,76,77]
pcp_1985_2019_trimed=pcp_1985_2019[pcp_1985_2019.columns[cols]]

# plots
sns.set(rc={'figure.figsize':(20,10)})
td=pcp_1985_2019_trimed[40].plot(linewidth=0.75, color='b')
td.set_xlabel('Date')
td.set_ylabel('Precipitation(mm)')



# average precipitation computation
weights=pd.Series([0.168,0.203,0.179,0.107,0.105,0.136,0.104])
#pcp_1985_2019_trimed.multiply(weights,axis="columns", level=None, fill_value=None)

pcp_1985_2019_trimed[40]=0.168*pcp_1985_2019_trimed[40]
pcp_1985_2019_trimed[53]=0.168*pcp_1985_2019_trimed[53]
pcp_1985_2019_trimed[64]=0.168*pcp_1985_2019_trimed[64]
pcp_1985_2019_trimed[65]=0.168*pcp_1985_2019_trimed[65]
pcp_1985_2019_trimed[75]=0.168*pcp_1985_2019_trimed[75]
pcp_1985_2019_trimed[76]=0.168*pcp_1985_2019_trimed[76]
pcp_1985_2019_trimed[77]=0.168*pcp_1985_2019_trimed[77]

# average precipitation in the catchment
pcp_final=pcp_1985_2019_trimed.sum(axis=1, skipna=True)

pcp_pr=pcp_final.to_frame()
pcp_pr.columns=['pcp']
pcp_pr['year']=pcp_pr.index.year

pcp_ams=pcp_pr.groupby(pcp_pr['year']).max()
pcp_ams.sort_values(by='pcp', ascending=False, inplace=True)
pcp_ams.to_excel(r'D:\Research Work\Synchronization work\Multivariate IDF\PCP_IMD_2019\pcp_ams.xlsx')

sns.set(rc={'figure.figsize':(20,10)})
td=pcp_final.plot(linewidth=0.75, color='b')
td.set_xlabel('Date')
td.set_ylabel('Precipitation(mm)')

# exporting pandas dataframe to excel
pcp_final.to_excel(r'D:\Research Work\Synchronization work\Multivariate IDF\PCP_IMD_2019\pcp_final.xlsx')

pcp_final['year']=pcp_final.index.year


# Exporting IMD grids in periyar basin
pcp_1985_2019_trimed.to_excel(r'D:\Research Work\Synchronization work\Multivariate IDF\PCP_IMD_2019\pcp_1985_2019_periyar.xlsx')

#%%%
# extracting the annual maxima series

pcp_periyar=pd.read_excel('F:/Research Work/Synchronization work/Multivariate IDF/PCP_IMD_2019g/pcp_1985_2019_periyar.xlsx', sheet_name='pcp_periyar', index_col=0, parse_dates=True)

pcp_periyar['year']=pcp_periyar.index.year

ams_1=pcp_periyar.groupby(pcp_periyar['year']).max()
ams_1days=pcp_periyar.groupby(pcp_periyar['year']).idxmax(axis=0)

# 3day rolling total
pcp_3d=pcp_periyar.rolling(3).sum()
pcp_3d['year']=pcp_3d.index.year
ams_3days=pcp_3d.groupby(pcp_3d['year']).idxmax(axis=0)

# ams 3day values
ams_3=pcp_3d.groupby(pcp_3d['year']).max()


############################################################
# 5 day discharge 
pcp_5d=pcp_periyar.rolling(5).sum()
pcp_5d['year']=pcp_5d.index.year

# ams 3day index values
ams_5days=pcp_5d.groupby(pcp_5d['year']).idxmax(axis=0)

# ams 3day values
ams_5=pcp_5d.groupby(pcp_5d['year']).max()

## writing into the excel files
ams_1.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_pcp/ams_1_pcp.xlsx")

ams_1days.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_pcp/ams_1day_pcp.xlsx")

ams_3.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_pcp/ams_3_pcp.xlsx")

ams_3days.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_pcp/ams_3day_pcp.xlsx")

ams_5.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_pcp/ams_5_pcp.xlsx")

ams_5days.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_pcp/ams_5day_pcp.xlsx")

pcp_3d.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_pcp/pcp_3day.xlsx")

pcp_5d.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_pcp/pcp_5day.xlsx")
