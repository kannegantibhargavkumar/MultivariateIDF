# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 10:43:36 2021

@author: Bhargav
"""

# In this files co-occurence of extreme sea level and extreme streamflow levels
# will be analysed.

#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
from pandas.plotting import autocorrelation_plot
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
from statsmodels.tsa.stattools import adfuller
import seaborn as sns

#%% 

tidal=pd.read_excel(r'D:\Research Work\Synchronization work\Multivariate IDF\Tidal data/daily_corrected.xlsx',index_col=0, parse_dates=True)
discharge=pd.read_excel(r'D:\Research Work\Synchronization work\Multivariate IDF\Kerala Floods Data_to Bhargav\Discharge\SF_data\1985_2021_final.xlsx',sheet_name='1985_2021',index_col=0, parse_dates=True)
skewsurge=pd.read_excel(r'D:\Research Work\Synchronization work\Multivariate IDF\Tidal data\cochin_skew_surge.xlsx',sheet_name='skewsurge',index_col=0, parse_dates=True )
pcp_data=pd.read_excel(r'D:\Research Work\Synchronization work\Multivariate IDF\PCP_IMD_2019\pcp_final.xlsx',sheet_name='rf',index_col=0, parse_dates=True)
tidal_in=pd.read_excel(r'D:\Research Work\Synchronization work\Multivariate IDF\Tidal data\dailymaxtide.xlsx', sheet_name='Data_clean',index_col=0,parse_dates=True)

sea_level=tidal.drop(['Year','Month','Day'], axis=1)

tide_pcp_gesla=pd.concat([sea_level,pcp_data],axis=1, join='inner') # global dataset
tide_pcp=pd.concat([tidal_in,pcp_data],axis=1,join='inner')  # incois dataset
tide_dis=pd.concat([tidal_in,discharge],axis=1,join='inner') # incois dataset

tide_dis.corr(method='kendall')

# correlation in monsoon period
list_months=[6,7,8,9]
list_nms=[1,2,3,4,10,11,12]
tide_dis_jjas=tide_dis[tide_dis.index.to_series().dt.month.isin(list_months)]
tide_dis_nms=tide_dis[tide_dis.index.to_series().dt.month.isin(list_nms)]

tide_dis_jjas.corr(method='kendall')
tide_dis_nms.corr(method='kendall')

#tide_dis_jjas.to_excel(r'D:\Research Work\Synchronization work\Multivariate IDF\Tidal data\Monsoon data\tide_dis_jjas.xlsx')
#%%
# reloading the cleaned data
tide_dis_jjas_clean=pd.read_excel(r'D:\Research Work\Synchronization work\Multivariate IDF\Tidal data\Monsoon data\tide_dis_jjas.xlsx',sheet_name='data_f',index_col=0, parse_dates=True)
tide_dis_jjas_clean.corr(method='kendall')
#checking the correlation with precipitation also
# =============================================================================
# # extracting monsoon data of tide and pcp
# tide_pcp_jjas=tide_pcp[tide_pcp.index.to_series().dt.month.isin(list_months)]
# tide_pcp_jjas_clean=pd.concat([tide_dis_jjas_clean['Tide'],tide_pcp_jjas['pcp']],axis=1,join='inner')
# tide_pcp_jjas_clean.corr(method='kendall')
# 
# =============================================================================
# Filling the nan values with rolling mean
tide_dis_jjas_final= tide_dis_jjas_clean.fillna(tide_dis_jjas_clean.rolling(6,min_periods=1).mean())
#tide_dis_jjas_final.to_excel(r'D:\Research Work\Synchronization work\Multivariate IDF\Tidal data\Monsoon data\tide_dis_jjas_final.xlsx')

#%%
# creating a dataframe with lags 
def df_derived_by_shift(df,lag=0,NON_DER=[]):
    df = df.copy()
    if not lag:
        return df
    cols ={}
    for i in range(1,lag+1):
        for x in list(df.columns):
            if x not in NON_DER:
                if not x in cols:
                    cols[x] = ['{}_{}'.format(x, i)]
                else:
                    cols[x].append('{}_{}'.format(x, i))
    for k,v in cols.items():
        columns = v
        dfn = pd.DataFrame(data=None, columns=columns, index=df.index)    
        i = 1
        for c in columns:
            dfn[c] = df[k].shift(periods=i)
            i+=1
        df = pd.concat([df, dfn], axis=1)
        df=df.reindex(df.index)
    return df

tide_dis_jjas_final['index']=tide_dis_jjas_final.index
NON_DER = ['index',]
df_new = df_derived_by_shift(tide_dis_jjas_final, 6, NON_DER)
#df_new=df_new.drop(['Discharge','Tide'],axis=1)
df_new = df_new.dropna()
corr_lag=df_new.corr()

# heat map of correlation
colormap = plt.cm.RdBu
plt.figure(figsize=(15,10))
plt.title(u'Cross and Auto Correlation', y=1.05, size=16)

mask = np.zeros_like(df_new.corr(method='kendall')) # default method is pearson
mask[np.triu_indices_from(mask)] = True

svm = sns.heatmap(df_new.corr(method='kendall'), mask=mask, linewidths=0.1,vmax=1.0, 
            square=True, cmap=colormap, linecolor='white', annot=True)

# plot of correlation line 
x=[-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6]
# these values are extracted from correlation matrix computed eariler
corr=[0.12,0.14,0.15,0.18,0.19,0.21,0.22,0.24,0.23,0.22,0.22,0.22,0.23]

plt.plot(x,corr,linewidth=2, markersize=12)
plt.scatter(1,0.24, color='r')
plt.ylim([0.1,0.26])
plt.xlabel('Time(Days)')
plt.ylabel('Kendall Tau')


q_1=tide_dis_jjas_final.quantile(q=0.95,axis=0)
q_2=tide_dis_jjas_final.quantile(q=0.99,axis=0)


plt.scatter(tide_dis_jjas_final['Tide'],tide_dis_jjas_final['Discharge'], color='b',s=25)
plt.axvline(x=1.24275, color='g', linestyle='--')
plt.axhline(y=1212.13,color='g', linestyle='--')
plt.axvline(x=1.36115, color='r', linestyle='--')
plt.axhline(y=2441.36,color='r', linestyle='--')
plt.xlim([0.5,1.8])
plt.ylim([0,7000])
plt.xlabel('Daily maximum Tide Level(m)')
plt.ylabel('Daily Discharge (Cumec)')
plt.show()



# autocorrelation plots
plt.figure(figsize=(20,10))
plot_acf(tide_dis_jjas_final['Discharge'],lags=48)
plt.show()

plt.figure(figsize=(20,10))
plot_acf(tide_dis_jjas_final['Tide'],lags=48)
plt.show()

plt.figure(figsize=(20,10))
plot_pacf(tide_dis_jjas_final['Discharge'],lags=48)
plt.show()

plt.figure(figsize=(20,10))
plot_pacf(tide_dis_jjas_final['Tide'],lags=48)
plt.show()



#%% NON monsoon season




#%%
# lagged correlation during monsoon  
# 
# def crosscorr(datax,datay,lag=0):
#     
#     #lag n cross correlation 
#     # parameters-lag: int;
#     #default 0
#     #datax,datay: pandas.Series objects of equal length 
#     #retuns
#     # crosscorr:float
#     return datax.corr(datay.shift(lag))
# lag1_corr=crosscorr(tide_dis_jjas_clean['Tide'],tide_dis_jjas_clean['Discharge'],lag=1)
# 
# correlation=[crosscorr(tide_dis_jjas_clean['Tide'],tide_dis_jjas_clean['Discharge'],lag=i) for i in range (5)]
# 
# 
# 
# 
# # autocorrelation
# ax1=autocorrelation_plot(tide_dis_jjas_clean['Discharge'])
# ax1.set.xlim([0,40])
# ax2=autocorrelation_plot(tide_dis_jjas_clean['Tide'])
# 
# plot_acf(tide_dis_jjas_clean['Discharge'], lags=40);
# plot_pacf(tide_dis_jjas_clean['Discharge'], lags=40)
# 
# plot_acf(tide_pcp_jjas_clean['Tide'], lags=40)
# 
# filled_dataset = tide_dis_jjas_clean.fillna(tide_dis_jjas_clean.rolling(6,min_periods=1).mean())
# 
# 
# 
# 
# 
# plot_acf(filled_dataset['Tide'],lags=40)
# plot_pacf(filled_dataset['Tide'],lags=40)
# # non stationarity test
# result_dis = adfuller(filled_dataset['discharge'])
# print('ADF Statistic: %f' % result_dis[0])
# print('p-value: %f' % result_dis[1])
# print('Critical Values:')
# for key, value in result_dis[4].items():
# 	print('\t%s: %.3f' % (key, value))
#     
# result_tide = adfuller(filled_dataset['Tide'].to_numpy())
# print('ADF Statistic: %f' % result_tide[0])
# print('p-value: %f' % result_tide[1])
# print('Critical Values:')
# for key, value in result_tide[4].items():
# 	print('\t%s: %.3f' % (key, value))