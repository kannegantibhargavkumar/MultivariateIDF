# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 22:09:15 2021

@author: Bhargav
"""

#%%
import pandas as pd
import pyextremes as py
import os
filedir=(r'F:\Research Work\Synchronization work\Multivariate IDF\Kerala Floods Data_to Bhargav\Discharge\SF_data')
os.chdir(filedir)
import numpy as np
import pylab as pyl


#%%

temp=pd.read_excel("F:/Research Work/Synchronization work/Work Related/Data and work/codes and results/SF data_Peninsular/Godavari/temp.xlsx", index_col=0,parse_dates=True,squeeze=(True))
temp=temp[(temp!='-')]
temp=temp.dropna()


model=py.EVA(data=temp)

model.get_extremes(method="BM",extremes_type="high", block_size="365.2425D", errors="raise")

model.plot_extremes()

model.fit_model()

model.plot_diagnostic(alpha=0.95)

#%%
# peak over thershold apporach 

Q_data=(pd.read_excel("1985_2021_final.xlsx",sheet_name='1985_2021', index_col=0, parse_dates=True, squeeze=True)
                     .sort_index(ascending=True)
                     .astype(float)
                     .dropna())

Q_data=Q_data.drop(["Unnamed:2"], axis=1)
Q_data=Q_data.sort_index(ascending=True)
Q_data=Q_data.astype(float)
 # now we need to choose thershold parameter
 # mean residual plot
 ax=py.plot_mean_residual_life(Q_data)
 fig=ax.get_figure()
 fig.savefig("mean-residual-life-high.png",dpi=96, bbox_inches="tight")
 # parameter stability
ax_shape,ax_scale=py.plot_parameter_stability(Q_data,alpha=None, progress=True )
fig=ax_shape.get_figure()
fig.savefig("parameter-stability.png",dpi=96, bbox_inches="tight")
#
# return value stability
#return value stability plot
ax=py.plot_return_value_stability(Q_data,
                                       return_period=100,
                                       return_period_size="365.24245D",
                                       thresholds=np.linspace(1000,1500,50),
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
axes=py.plot_threshold_stability(Q_data,return_period=100,thresholds=np.linspace(1000,1500,50),progress=True)
fig=axes[0].get_figure()
fig.savefig("threshold-stability.png",dpi=96,bbox_inches="tight")
# model fitting
model2=py.EVA(Q_data)
model2.get_extremes("POT",threshold=1500)
fig=ax.get_figure()
ax=model2.plot_extremes()
fig.savefig("POT series", dpi=96,bbox_inches="tight")

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
    return_period=list(np.linspace(1,500,50)),
    alpha=0.95,
    n_samples=1000,
    )

pyl.plot(np.linspace(1,500,50),summary["return value"],'bo')
pyl.xlabel('Return Period(yrs)')
pyl.ylabel("Discharge (Cumec)")



