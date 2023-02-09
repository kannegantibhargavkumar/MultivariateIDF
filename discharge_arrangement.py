# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 14:44:48 2021

@author: Bhargav
"""
#%%
import pandas as pd
#import pyextremes as py
import matplotlib.pyplot as plt
import os


import seaborn as sns
import numpy as np
#import scipy 
#import scipy.stats
from scipy.stats import genextreme
from numpy import linspace
#from pylab import plot, show, hist, figure,title
import pylab



#%%
filedir=(r'F:\Research Work\Synchronization work\Multivariate IDF\Kerala Floods Data_to Bhargav\Discharge\SF_data')
os.chdir(filedir)
discharge=pd.read_excel('1985_2021_final.xlsx', sheet_name='1985_2021', index_col=0, parse_dates=True)
# indexing the dataframe by date series
# extracting 1day annual maxima

discharge['year']=discharge.index.year
discharge['month']=discharge.index.month
discharge['day']=discharge.index.day

ams_1_days=discharge.groupby(discharge['year']).idxmax(axis=0)
ams_1_days=ams_1_days.drop(['month','day'], axis=1)
ams_1_days.rename(columns={'discharge':'day'}, inplace = True)
# ams values
ams_1=discharge.groupby(discharge['year']).max()
ams_1=ams_1.drop(['month','day'], axis=1)
# combine two dataframes
ams_dis_1day=pd.concat([ams_1,ams_1_days], axis=1, join='inner')

############################################################
# 3 day discharge
dis_3d=discharge['discharge'].rolling(3).sum()
dis_3d=dis_3d.to_frame()
dis_3d['year']=dis_3d.index.year
dis_3d['month']=dis_3d.index.month
dis_3d['day']=dis_3d.index.day
# ams 3day index values
ams_3days=dis_3d.groupby(dis_3d['year']).idxmax(axis=0)
ams_3days=ams_3days.drop(['month','day'], axis=1)
ams_3days.rename(columns={'discharge':'day'}, inplace = True)
# ams 3day values
ams_3=dis_3d.groupby(dis_3d['year']).max()
ams_3=ams_3.drop(['month','day'], axis=1)
# combine two dataframes
ams_dis_3day=pd.concat([ams_3,ams_3days], axis=1, join='inner')

############################################################
# 5 day discharge 
dis_5d=discharge['discharge'].rolling(5).sum()
dis_5d=dis_5d.to_frame()
dis_5d['year']=dis_5d.index.year
dis_5d['month']=dis_5d.index.month
dis_5d['day']=dis_5d.index.day
# ams 3day index values
ams_5days=dis_5d.groupby(dis_5d['year']).idxmax(axis=0)
ams_5days=ams_5days.drop(['month','day'], axis=1)
ams_5days.rename(columns={'discharge':'day'}, inplace = True)
# ams 3day values
ams_5=dis_5d.groupby(dis_5d['year']).max()
ams_5=ams_5.drop(['month','day'], axis=1)
# combine two dataframes
ams_dis_5day=pd.concat([ams_5,ams_5days], axis=1, join='inner')

##### saving excel files

#1 day excel file
ams_dis_1day.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_new/ams_dis_1day.xlsx")
# 3 day excel file

dis_3d.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_new/dis_3day.xlsx")

ams_dis_3day.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_new/ams_dis_3day.xlsx")

# 5 day excel file
dis_5d.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_new/dis_5day.xlsx")

ams_dis_5day.to_excel("F:/Research Work/Synchronization work/Multivariate IDF/AMS/ams_new/ams_dis_5day.xlsx")




# convert the column(it's a string) to datetime type
datetime_series=pd.to_datetime(discharge['Date'])
# create datetime index passing the datetime series
datetime_index=pd.DatetimeIndex(datetime_series.values)
discharge_1985_2018=discharge.set_index(datetime_index)
# we don't need the date column anymore
discharge_1985_2018.drop('Date',axis=1, inplace=True)
pot_discharge=discharge_1985_2018.loc[discharge_1985_2018['Discharge']>1500,:]
