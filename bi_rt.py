# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 21:19:13 2022

@author: bhargav
"""

#this function extracts the annual maxima of precipitation and streamflow
#%%
import pandas as pd

#%%
sf=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Data/northatlantic streamflowdata/sf_6608300.xlsx", index_col=2, parse_dates=(True))

rf=pd.read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Data/Data_at_cbridge/rainfall_cbridge.xlsx", index_col=0, parse_dates=True)

ams_sf=sf["runoff_mean"].groupby(sf.index.year).max()
ams_rf=rf.groupby(rf.index.year).max()
ams_sf=pd.DataFrame(ams_sf)
ams=ams_sf.join(ams_rf, how="inner")
ams.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Data/Data_at_cbridge/ams_cbridge.xlsx")
