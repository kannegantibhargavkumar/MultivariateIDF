# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 23:48:44 2022

@author: bhargav
"""

import pandas as pd
import numpy as np
import os
os.chdir("F:\Research Work\Synchronization work\Streamflow Data\Swapan\WRIS Data")

import glob as glob
ff=pd.ExcelFile("Stations _all_d.xlsx")
### useful command for merging all sheets 
data=pd.concat(pd.read_excel("Stations _all_d.xlsx", sheet_name=None), ignore_index=True)

data.to_csv("stations.csv")


#%%
# Reading all xlsx files names
files_list=[]

for root, directories, files in os.walk("F:\Research Work\Synchronization work\Streamflow Data\Swapan\WRIS Data"):
    for name in files:
        files_list.append(os.path.join(root, name))
        

print(files_list)

##alternatively
file_list2=glob.glob("F:/Research Work/Synchronization work/Streamflow Data/Swapan/WRIS Data\*.csv")
## finding string in list
gauge_files=[match for match in files_list if "Gauge" in match]

#%%
## search stations




stations=pd.read_csv("F:/Research Work/Synchronization work/Streamflow Data/stations_new.csv")
st=stations["Station"]
dateindex=pd.date_range(start="01/01/1965", end="31/12/2018")
sf_data=pd.DataFrame()
data_extremes=pd.DataFrame()
data_extremes.index=range(len(gauge_files))
#sf_extremes=pd.DataFrame()
sf_data.index=dateindex
for i in range(len(gauge_files)):
    file=gauge_files[i]
    temp=pd.read_csv(file, sep='delimiter', header=None)
    temp=pd.read_csv(file, skiprows=2, index_col=0, parse_dates=True)
    tp=temp["Discharge (cumecs)"]
    #te=pd.merge(stations.iloc[i,1], tp.max())
    #sf_extremes=pd.concat([sf_extremes,])
    t=tp.max()
    st=stations.iloc[i,1:4]
    st.loc[4,]=t
    df=pd.DataFrame(list(st))
    df=df.transpose()
    data_extremes=pd.concat([data_extremes,df])
    #sf_data=pd.concat([sf_data,(tp.max(), stations.iloc[i,1])], axis=0)
    
## remove the nan rows 
data_extremes=data_extremes.dropna(axis=0)
data_extremes.columns=["Station","Lat","Long","Max_SF"]
data_extremes.index=range(329)
data_extremes.to_csv("F:\Research Work\Synchronization work\Regional envelope curve apporach\data\sf_extremes.csv")
