# -*- coding: utf-8 -*-
"""
Created on Fri May 27 05:58:58 2022

@author: bhargav
"""
# This function creates grid points given lat and lon series
# Input: lat and lon np series
# output: grid points dataframe
def latlon_to_points(lat,lon):
    import numpy as np
    import pandas as pd
    ## initializing the array
    t_lat=lat[0]
    tt_lat=np.array([t_lat]*len(lon))
    pts=np.column_stack((tt_lat,lon))

    for i in range(1,len(lat)):
        t_lt=lat[i]
        tt_lt=np.array([t_lt]*len(lon))
        k=np.column_stack((tt_lt,lon))
        pts=np.vstack((pts,k))
    grid_points=pd.DataFrame(pts)
    grid_points.columns=["lat","lon"]
    
    return(grid_points)

    
