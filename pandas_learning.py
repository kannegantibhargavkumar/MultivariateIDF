# -*- coding: utf-8 -*-
"""
Created on Tue May 17 18:41:18 2022

@author: bhargav
"""
#%%
import pandas as pd
import matplotlib.pyplot as plt

#%%
# creating a dataframe
df=pd.DataFrame()
# checking type of variable
type(df)

# random data
name=["aish","bhargav","ppm"]
age=[25, 29,65]
occupation=["software","student","professor"]

df["name"]=name
df["age"]=age
df["occupation"]=occupation
# alternative way
#df_a=pd.DataFrame({
 #   "name":["aish","bhargav","ppm"],
  #  "age":[]})

#%%
# finding min and maximum in series
max_age=df["age"].max()
min_age=df["age"].min()

#%%
# reading tabular data
prime_titles=pd.read_excel("F:/Aishu_job/titles.csv/titles.xlsx")
# type of variables
prime_titles["runtime"].dtypes

# exporting table
prime_titles.to_excel()
