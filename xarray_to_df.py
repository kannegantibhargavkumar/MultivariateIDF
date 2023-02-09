# -*- coding: utf-8 -*-
"""
Created on Thu May 26 12:44:56 2022

@author: bhargav
"""
#xarray_to_df.py>
# function file for extracting data from xarray datasets
# input: 1) xarray.core.dataset.Dataset 2) list of lon and lat values
# output: pandas dataframe with datetime index , each column represents
#         data at each grid point
def xr_to_df(df, lat, lon):
    import pandas as pd
    import xarray as xr
    import numpy as np
    h=list(df.dims)
    df=df.rename({h[1]:"lat",h[2]:"lon"})
    data=pd.DataFrame()
    data.index=df.time
    n=np.arange(4)
    for i in range(len(lat)):
        data[str(n[i])]=df.sel(lat=lat[i], lon=lon[i])
        
    return(data)

        
        
        
        