# -*- coding: utf-8 -*-
"""
Created on Tue May 31 10:05:21 2022

@author: bhargav
"""
#%%
import pandas as pd
import xarray as xr
import os
path=os.chdir()
from glob import glob

#%%
# reading the nc file names in the folder and sub folder
filenames=glob("F:\Research Work\Camles dataset\*.nc")
# converting the list to string
k=''.join(filenames)
# length of string is 50 




data=xr.open_dataset("camels_01022500.nc")
sf_data=data["streamflow"]
sf=sf_data.to_dataframe()
