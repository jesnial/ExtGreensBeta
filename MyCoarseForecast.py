# -*- coding: utf-8 -*-
"""
Created on Tue May 21 12:53:41 2019

@author: Jessie
"""


from netCDF4 import Dataset
import numpy as np
import os
import math


def coarse_forecast(SimRes, Bath, lowlim=-500, highlim=-50,hE=100):
    '''
    wave height from the coarse grid simulation for all points : SimRes
    2 dim for (x,y) and 10 dim(SimRes)=(number of deg_east, number of deg_north, 9)
    bathymetry table : Bath  dim(bath)=(number of deg_east, number of deg_north, 1)
    epsilon : deep water bathymetry limit
    hE :a parameter value for when the extended green’s law is chosen to be used from (approx. 100m)
    '''
    east,north,N =np.shape(SimRes)
   
    Computed_heights = np.copy(SimRes)
    #copy results for bathy>500m
    
    for i in range (east):
        for j in range (north):
            #use Green's law between 500 and 50m
            
                if lowlim<Bath[i,j]<highlim:
                    #find nearest point at greater depth 
                    k=-1
                    while ( all(np.isnan(a) for a in SimRes[i,np.mod((j+k),north)]) or Bath[i,j]<Bath[i,np.mod((j+k),north)]) and k<north-1:
                        k=k+1
                        if k==0:
                            k=k+1
                    for n in range (N):
                        #use Green's law from this nearest point
                        Computed_heights[i,j,n]= Computed_heights[i,np.mod((j+k),north),n]*(Bath[i,np.mod((j+k),north)]/Bath[i,j])**(1/4)   
                               
                if 0>Bath[i,j]>highlim: #0 if bathy<50m
                    for n in range (N):
                        Computed_heights[i,j,n]=0

                #print(Computed_heights[i,j,:], "computed")
                #print(SimRes[i,j,:], "simres")
                #print(i,j)
    return Computed_heights        
                    
    

