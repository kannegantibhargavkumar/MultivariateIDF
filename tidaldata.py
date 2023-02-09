# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 12:07:03 2021

@author: Bhargav
"""

# new data from INCOIS 
#%%
import os
filedir=(r'D:\Research Work\Synchronization work\Multivariate IDF\Tidal data\Cochin')
os.chdir(filedir)
import glob
import pandas as pd
import numpy as np
#%%
filenames2=[i for i in glob.glob("*.xlsx")] 
# merging all years data into time series
frames=[pd.read_excel(i, index_col=1, parse_dates=True) for i in filenames2]
result=pd.concat(frames)
tidal_data=result.drop(['station','Station'],axis=1)
tidal_data.replace(9999.99,np.nan)



tidal_data['day']=tidal_data.index.day

re_tidal=tidal_data.resample('D').max() # minute data to daily data

re_tidal=re_tidal.drop('day',axis=1)
re_tidal.to_excel(r'D:\Research Work\Synchronization work\Multivariate IDF\Tidal data\dailymaxtide.xlsx')
