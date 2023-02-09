# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 11:00:18 2021

@author: Bhargav
"""
#!pip install netcdf4
#!pip install xarray
#%%
import netCDF4 as nc
from netCDF4 import Dataset
import pandas as pd
import os
#filedir=("D:\Research Work\Synchronization work\Multivariate IDF\soil moisture')
#os.chdir(filedir)
import glob
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# import geopandas 
# import xarray
# import rasterio
# import glob
# =============================================================================

#%%
# importing the soil moisture data for india 

ds=xr.open_mfdataset("F:/Research Work/Synchronization work/soil moisture data/1984_2020sm_to_IMD_grid/1984_2020sm_to_IMD_grid.nc")
latt=ds.lat.values
long=ds.lon.values













#%%%

# # india boundaries
# lat_ind=latt[248:329]
# long_ind=long[991:1112]


# grid_points=np.transpose(np.reshape(np.meshgrid(lat_ind,long_ind), (2,9801)))

# grid_ind_sm=pd.DataFrame(grid_points)
# grid_ind_sm.to_excel(r'D:\Research Work\Synchronization work\Multivariate IDF\soil moisture\grids_ind_sm.xlsx')

# ds2=xr.open_mfdataset(r'D:\Research Work\Synchronization work\Multivariate IDF\PCP_IMD_2019\rain_1951_2020.nc')

# ds3=ds.regrids(ds2,'conservative')


# ds1=cf.read('1984_2019sm.nc')[0]

# # =============================================================================
# 
# # periyarbasin boundaries
# k=latt[318:324]
# j=long[1023:1030]
# 
# # clipping the soil moisture data of periyar basin
# sm_periyar=ds.sm.sel(lat= k, lon=j)
# 
# gd=pd.read_excel('lat_long_details.xlsx',sheet_name='grids')
# 
#   
# sm_periyar.to_dataframe()
# 
# gd=pd.read_excel('lat_long_details.xlsx',sheet_name='grids')
# 
# =============================================================================

 # handling the arcgis functions here
 # periyar=sp.Reader(r'D:\Research Work\Synchronization work\Multivariate IDF\Kerala Floods Data_to Bhargav\Shapefiles\periyar basin boundary\Periyar basin boundary\FINAL_BOUNDARYKFP.shp')