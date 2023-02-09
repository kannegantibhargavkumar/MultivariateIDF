# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 21:15:21 2022

@author: bhargav
"""
#%%
import xarray as xr
import pandas as pd
#%%

sm=xr.open_dataset("F:/Research Work/Synchronization work/soil moisture data/1984_2020sm_to_IMD_grid/1984_2020sm_to_IMD_grid.nc")

# renaming the data variable
sm['sm'] = sm['__xarray_dataarray_variable__']

# dropping the duplicate variable
sm= sm.drop(['__xarray_dataarray_variable__'])

temp=sm.sm.where((sm.sm.lat==17.25) & (sm.sm.lon==74) , drop=True)

t=temp.to_dataframe()
