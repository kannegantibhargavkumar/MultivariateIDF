# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 11:29:34 2022

@author: bhargav
"""

#%%
import pandas as pd
import xarray as xr
#%%

dt=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/ganga_work/sf_data.xlsx", index_col=0, parse_dates=(True))
dt_sf=dt.groupby(dt.index.year).max()
dt_sf.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/ganga_work/ams_joshimath.xlsx")
