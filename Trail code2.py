# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 10:09:01 2021

@author: Bhargav
"""

import matplotlib.pyplot as plt
import scipy 
import scipy.stats

size=20000
x=scipy.arange(size)
# creating the dummy sample (using beta distribution)

y=scipy.int_(scipy.round_(scipy.stats.beta.rvs(6,2,size=size)*47))

# creating the histogram 

h=plt.hist(y,bins=range(48))
dist_names=['alpha','beta','arcsine','weibull_min','weibull_max','rayleigh']

for dist_name in dist_names:
    dist=getattr(scipy.stats, dist_name)
    params=dist.fit(y)
    arg=params[:-2]
    loc=params[-2]
    scale=params[-1]
    if arg:
        pdf_fitted=dist.pdf(x,*arg, loc=loc,scale=scale)*size
    else:
        pdf_fitted=dist.pdf(x,loc=loc,scale=loc)*size
    plt.plot(pdf_fitted, label=dist_name)
    plt.xlim(0,47)
plt.legend(loc='upper left')
plt.show()

