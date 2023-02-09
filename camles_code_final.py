# -*- coding: utf-8 -*-
"""
Created on Tue May 31 10:05:21 2022

@author: bhargav
"""
#%%
import pandas as pd
import xarray as xr
import os
path=os.chdir("F:\\Research Work\\Camles dataset")
from glob import glob
from scipy.stats import genextreme as gev
import csv
import pickle
#%%
# reading the nc file names in the folder and sub folder
filenames=glob("F:\Research Work\Camles dataset\*.nc")
# converting the list to string
k=''.join(filenames)
data=xr.open_dataset("camels_01022500.nc")
sf_data=data["streamflow"]
sf=sf_data.to_dataframe()
ams=sf.groupby(sf.index.year).max()
ams=sf.groupby(sf.index.year).max().dropna()
f=gev.fit(ams)
cdf=gev.cdf(ams, f[0],f[1],f[2])
rt=[0.5,0.9, 0.98, 0.99]
sf_qunatile=gev.ppf(rt,f[0],f[1],f[2])
# the file name is sf_model_ catchmentid
fg={"data":sf,"ams":ams,"model":"GEV","param":f,"qunatile":sf_qunatile}
#globals()['sf_model_'+k[-11:-3]]={"data":sf,"ams":ams,"model":"GEV","param":f,"qunatile":sf_qunatile}
# writing the dictionary into pickle file
n='sf_model_'+k[-11:-3]+'.pkl'
# writing the dictionary into pickle format
f=open( n,"wb")
pickle.dump(fg, f)
f.close()


