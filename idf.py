# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 09:17:33 2021

@author: Bhargav
"""
# univariate analysis of streamflow series

#%%
#   !pip install pyextremes# # line for install for package
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
filedir=(r'D:\Research Work\Synchronization work\Multivariate IDF\Kerala Floods Data_to Bhargav\Discharge\SF_data')
os.chdir(filedir)


#import scipy 
#import scipy.stats
from scipy.stats import genextreme
from numpy import linspace
#from pylab import plot, show, hist, figure,title
import pylab
import pyextremes as py


#%%

# =============================================================================
# #dis_neel=pd.read_excel( 'Neeleswaram 1985-2018_SUMQH - Copy.xls')
# 
# #discharge=dis_neel.iloc[:,8] # 9th column is of interest
# #date=dis_neel.iloc[:,3]
# 
# #result=pd.concat([date,discharge],axis=1,join='inner')
# 
# #result.sort_values(by=['Discharge'], ascending=False, inplace=True)
# 
# #tr_data=result.iloc[1:100,]
# 
# # 
# #my_range=range(1,len(result.index)+1)
# #plt.stem(result['Discharge'])
# #plt.xticks(my_range,result['Date'], rotation='vertical')
# 
# =============================================================================

####
#%%
# =============================================================================
# trial_discharge=pd.read_excel(r'D:\Research Work\Synchronization work\Multivariate IDF\Kerala Floods Data_to Bhargav\Discharge\Neeleswarm_1985-2018.xlsx', index_col=3, parse_dates=True)
# trial_discharge['year']=trial_discharge.index.year
# trial_discharge['month']=trial_discharge.index.month
# trial_discharge['day']=trial_discharge.index.day
# 
# # plotting
# # line plot
# sns.set(rc={'figure.figsize':(20,10)})
# ts_d=trial_discharge['Discharge'].plot(linewidth=0.75);
# ts_d.set_ylabel('Daily Discharge (Cumec)')
# 
# # seasonality
# fig, axes=plt.subplots(1,1, figsize=(16,10),sharex=True)
# sns.boxplot(data=trial_discharge, x='month', y='Discharge', ax=axes)
# axes.set_ylabel('Daily Discharge (Cumec)')
# axes.set_title('Seasonality')
#   
# # sorted line plot
# trial_discharge.sort_values(by='Discharge', ascending=False, inplace=True)
# trial_discharge['Discharge'].plot(linewidth=0.75)
# =============================================================================
#%%

# the analysis starts from here
# import the data from excel sheet using pandas
# indexing the data with dates 
Q_neel=pd.read_excel('1985_2021_final.xlsx', sheet_name='1985_2021',  index_col=0, parse_dates=True)

Q_neel['year']=Q_neel.index.year
Q_neel['month']=Q_neel.index.month
Q_neel['day']=Q_neel.index.day

# plotting the time series of discharge using the sns package
sns.set(rc={'figure.figsize':(20,10)})
ts_d=Q_neel['discharge'].plot(linewidth=0.75)
ts_d.set_ylabel('Daily Discharge (Cumec)')


k=g.sort_values('discharge')

# annual maxima series
g=Q_neel.groupby(Q_neel['year']).max()
g['yr']=g.index.values
g.sort_values(by='discharge', ascending=False, inplace=True)
g['discharge'].stem(linewidth=0.75)

ams_sort=g.drop(['month', 'day'], axis=1)


ams_sort.to_excel(r'D:\Research Work\Synchronization work\Multivariate IDF\Kerala Floods Data_to Bhargav\Discharge\SF_data\ams_neel.xlsx')


g.yr.astype(str)

plt.figure(figsize=(10,6))
# bar plot with matplotlib
plt.bar('yr','discharge',data=g, width=0.1)
plt.xlabel("Year", size=15)
plt.ylabel("Discharge ", size=15)




plt.stem(g.index.values,g['discharge'].values, use_line_collection=True)


# =============================================================================
# fig, axes=plt.subplots(1,1, figsize=(12,10),sharex=True)
# hist(g['Discharge'])
# axes.set_xlabel('Daily Discharge (Cumec)')
# axes.set_ylabel('Frequency')
# axes.set_title('Block Maxima')
# 
# =============================================================================

param=genextreme.fit(g['Discharge'])

# using pyextremes 
# Initialize the EVA object



## 

 dis_data=(pd.read_excel("1985_2021_final.xlsx",sheet_name='1985_2021', index_col=0, parse_dates=True, squeeze=True)
                     .sort_index(ascending=True)
                     .astype(float)
                     .dropna())
 
# =============================================================================
# # =============================================================================
# dis_data.head()
# dis_data.plot()
# =============================================================================

# intialize the model

model=py.EVA(data=dis_data)
# Extract the extreme values
model.get_extremes(method='BM',extremes_type="high",
                   block_size="365.24245D",
                   errors="ignore")
model.plot_extremes()
model.fit_model()
model.plot_diagnostic(alpha=0.95)
model.plot_return_values(return_period=np.logspace(0.01,2,1000),
                         return_period_size="365.2425D",
                         alpha=0.95,
                         )
summary=model.get_summary(
    return_period=[2,5,10,25,50,100,250,500,1000],
    alpha=0.95,
    n_samples=1000,
    )

# =============================================================================
# plot(summary['return value'])
# plot.xlabel('Return Period (Years)')
# plot.ylabel('Discharge (Cumec)')
# 
# =============================================================================
AMS=model.extremes

id=Q_neel.index[Q_neel['discharge']>=1500].tolist()


# alternatively POT apporach
# thershold selection
# mean residual plot
ax=py.plot_mean_residual_life(ts=dis_data)
fig=ax.get_figure()
fig.savefig("Mean-residual-life-high.png",dpi=96,bbox_inches="tight")
# parameter stability plot
ax_shape,ax_scale=py.plot_parameter_stability(ts=dis_data,alpha=None, progress=True )
fig=ax_shape.get_figure()
fig.savefig("parameter-stability.png",dpi=96, bbox_inches="tight")

#return value stability plot
ax=py.plot_return_value_stability(ts=dis_data,
                                       return_period=100,
                                       return_period_size="365.24245D",
                                       thresholds=np.linspace(1250,1750,50),
                                       r="24H",
                                       extremes_type="high",
                                       distributions=["genpareto","expon"],
                                       alpha=0.95,
                                       n_samples=100,
                                       progress=True,
                                       )

fig=ax.get_figure()
fig.savefig("return-value-stability.png",dpi=96,bbox_inches="tight")

# plot thrshold stability
axes=py.plot_threshold_stability(dis_data,return_period=100,thresholds=np.linspace(1250,2000,50),progress=True)
fig=axes[0].get_figure()
fig.savefig("threshold-stability.png",dpi=96,bbox_inches="tight")
# use selected threshold
model2=py.EVA(dis_data)
model2.get_extremes("POT",threshold=1500)
#ax=model2.plot_extremes()
#fig=ax.get_figure()
#fig.savefig("POT series", dpi=96,bbox_inches="tight")

model2.fit_model()
fig, ax=model2.plot_diagnostic(alpha=0.95)
fig.savefig("Selected-threshold-diagnostic.png",dpi=96,bbox_inches="tight")

fig, ax=model2.plot_return_values(return_period=np.logspace(0.01,2,100),
                         return_period_size="365.2425D",
                         alpha=0.95,
                         )
fig.savefig("IDF of 100 years.png",dpi=96,bbox_inches="tight")


# customized return levels
summary=model2.get_summary(
    return_period=list(linspace(1,500,50)),
    alpha=0.95,
    n_samples=1000,
    )

pylab.plot(linspace(1,500,50),summary["return value"],'bo')
pylab.xlabel('Return Period(yrs)')
pylab.ylabel("Discharge (Cumec)")


id2=dis_data.index[dis_data>=1500].tolist()


