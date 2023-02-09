# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 12:10:57 2022

@author: bhargav
"""

import pandas as pd
#%%
data=pd.read_csv("F:/naman/namandata.csv", index_col=3,parse_dates=(True))

#%%annual sum
data_a=data["Mean precipitation (mm)"].groupby(data.index.year).sum()

data_m=data["Mean precipitation (mm)"].resample("M").sum()

dd=data_m.values.reshape(56,12)
dd=pd.DataFrame(dd)
data_monthlymean=dd.mean(axis=0)



data_mm=data["Mean precipitation (mm)"].groupby(data.index.month).mean()

dd.plot()
