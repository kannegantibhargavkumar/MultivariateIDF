# -*- coding: utf-8 -*-
"""
Created on Fri May 27 06:40:08 2022

@author: bhargav
"""
    # This function calculates the circular statistics such as 
    # Mean date of occurence, peristance and variablity in date of 
    #occurence of extremes
    # input: time series of data with datetime index
    
    # 
def c_stats(data):
    
    
    import numpy as np
    import pandas as pd
    import math
    date=data.groupby(data.index.year).idxmax() # return date of maximum in a year
    #t=data.groupby(data.index.year).max().to_numpy() # returns the annual maxima
    t=data.iloc[:,1]
    #d=date.squeeze()
    #dy=d.dt. dayofyear.to_numpy() # dayofyear
    dy=data["day"]
    j=date.index.values # year
    
# intializing the for loop
    x=np.zeros(len(dy))
    y=np.zeros(len(dy))
    
    for i in range(len(dy)):
         # if leap year 
        if((j[i]%400==0)or(j[i]%100 != 0)and(j[i]%4 ==0)):
          len_yr=366
          theta= (dy[i]*(2*180))/len_yr
          x[i]=t[i]*math.cos(theta*math.pi/180)
          y[i]=t[i]*math.sin(theta*math.pi/180)
        else:
          len_yr=365
          theta= (dy[i]*(2*180))/len_yr
          x[i]=t[i]*math.cos(theta*math.pi/180)
          y[i]=t[i]*math.sin(theta*math.pi/180)
    
    # coordinates of extremes occurences
    
    x_bar=sum(x)/sum(t)
    y_bar=sum(y)/sum(t)
    
    # intializing the value
    ma=[]
    # converting into polar coordinates
    if (x_bar>=0 and y_bar>=0):
        ma=math.degrees(math.atan(y_bar/x_bar));
    elif (x_bar<0 and y_bar>0):
        ma=180-abs(math.degrees(math.atan(y_bar/x_bar)));
    elif (x_bar<0 and y_bar<0):
        ma=180+abs(math.degrees(math.atan(y_bar/x_bar)));
    elif (x_bar>=0 and y_bar<=0):
        ma=360-abs(math.degrees(math.atan(y_bar/x_bar)));
        
    #  # ma= Mean Angle , md= Mean date of occurence   
    md=ma*(365/(2*180))
    # r= persistance, sigma_sq= Variance
    r_bar=np.power((np.square(x_bar)+np.square(y_bar)),0.5)
    sigma_sq=-2*math.log(r_bar)
    
    s=pd.DataFrame()
    s["statistic"]=["x_bar","y_bar","ma","md","r_bar","sigma_sq"]
    
    s["value"]=list([x_bar, y_bar, ma, md, r_bar,sigma_sq])

    return s


