# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 11:15:44 2021

@author: Bhargav
"""

from scipy.stats import genpareto
from numpy import linspace
from pylab import plot, show, hist, figure,title

# picking 150 of from a normal distribution
# with mean 0 and standard deviatation 1
samp=genpareto.rvs(0.1,size=100)

param=genpareto.fit(samp) # distribution fitting


# now, param[0] and param [1] are the mean and the 
# standard deviation of the fitted distribution

x=linspace(0,5,100)

# fitted distribution
pdf_fitted=genpareto.pdf(x, param[0],param[1], param[2])

# orginal distribution
pdf=genpareto.pdf(x,0.1)

title('Normal distribution')
plot(x,pdf_fitted,'r-', x, pdf,'b-')
hist(samp,normed=1, alpha=0.3)
show()
