# -*- coding: utf-8 -*-
"""
Created on Tue May 21 12:10:14 2019

@author: Jessie#
"""

from netCDF4 import Dataset
import numpy as np
import os
import scipy.interpolate as scp
#import matplotlib.pyplot as plt
import scipy
from time import time
#plt.show()
#plt.ion()

from graddescents import *
from MyFineForecastNicePromtest import *
from MyCoarseForecast import *


#%% make all simulations from grid files into a numpy array

def make_array(place, tsu_list):
    N=len(tsu_list)
    bathy=Dataset(place,'r')
    var=list(bathy.variables.keys())
    x=bathy.variables[var[0]]
    x=np.asarray(x)
    y=bathy.variables[var[1]]
    y=np.asarray(y)
    z=bathy.variables[var[2]]
    Bath=np.asarray(z)
    (east,north)=np.shape(Bath)
    print(east,north)
    SimRes=np.zeros((east,north,N))
    for i in range (N):
        tsu=Dataset(tsu_list[i],'r')
        var=list(tsu.variables.keys())
        sim=tsu.variables[var[2]]
        SimRes[:,:,i]=np.asarray(sim)
    return x,y,Bath,SimRes

#%% 3h hmax for Bandol, load all simulations
place='Jessie/grids_extGreenLaw/bathy/West-Mediterranean.grd'
tsu_list=['Jessie/grids_extGreenLaw/HMAXnoClip_grd2min/cropped/hmax.gr00.3h_413-75s.grd','Jessie/grids_extGreenLaw/HMAXnoClip_grd2min/cropped/hmax.gr00.3h_413-65.grd','Jessie/grids_extGreenLaw/HMAXnoClip_grd2min/cropped/hmax.gr00.3h_413-70.grd','Jessie/grids_extGreenLaw/HMAXnoClip_grd2min/cropped/hmax.gr00.3h_413-75.grd','Jessie/grids_extGreenLaw/HMAXnoClip_grd2min/cropped/hmax.gr00.3h_J7.1.grd','Jessie/grids_extGreenLaw/HMAXnoClip_grd2min/cropped/hmax.gr00.3h_J7.8.grd','Jessie/grids_extGreenLaw/HMAXnoClip_grd2min/cropped/hmax.gr00.3h_L.grd','Jessie/grids_extGreenLaw/HMAXnoClip_grd2min/cropped/hmax.gr00.3h_S.grd','Jessie/grids_extGreenLaw/HMAXnoClip_grd2min/cropped/hmax.gr00.3h_Y.grd']

lat,lon,Bath,SimRes=make_array(place,tsu_list)

BathBandol_part1='Jessie/grids_extGreenLaw/bathy/FR_Bandol_part1_10m.grd'
listNice=['Jessie/grids_extGreenLaw/Resultat_Calc_multigrd/Casst413_Mw7.5s/hmax.gr03.4.00000_Bandol_part1.grd','Jessie/grids_extGreenLaw/Resultat_Calc_multigrd/Casst413_Mw7.0/hmax.gr03.4.00000_Bandol_part1.grd','Jessie/grids_extGreenLaw/Resultat_Calc_multigrd/Casst413_Mw7.5/hmax.gr03.4.00000_Bandol_part1.grd','Jessie/grids_extGreenLaw/Resultat_Calc_multigrd/Casst413_Mw6.5/hmax.gr03.4.00000_Bandol_part1.grd','Jessie/grids_extGreenLaw/Resultat_Calc_multigrd/Jijel7.1/hmax.gr03.4.00000_Bandol_part1.grd','Jessie/grids_extGreenLaw/Resultat_Calc_multigrd/Jijel7.8/hmax.gr03.4.00000_Bandol_part1.grd','Jessie/grids_extGreenLaw/Resultat_Calc_multigrd/Ligure1887/hmax.gr03.4.00000_Bandol_part1.grd','Jessie/grids_extGreenLaw/Resultat_Calc_multigrd/Boum2003_SemmaneMw7.1/hmax.gr03.4.00000_Bandol_part1.grd','Jessie/grids_extGreenLaw/Resultat_Calc_multigrd/Boum2003_YellesMw7.5/hmax.gr03.4.00000_Bandol_part1.grd']

finelat,finelon,FineBath,FineSimRes=make_array(BathBandol_part1,listNice)
np.save('Dataset/bathy/Bandol1_bathy.npy',FineBath)


#%% compute coarse forecast over the whole Mediterranean

tic=time()
Computed_heights=coarse_forecast(SimRes, Bath, lowlim=-500, highlim=-50,hE=100)
toc=time()
print(toc-tic)

##%% save output
np.save('Jessie/grids_extGreenLaw/MyRes/CoarseHeights',Computed_heights)

#%% load output if needed
Computed_heights = np.load('Jessie/grids_extGreenLaw/MyRes/CoarseHeights.npy','r')
 

#%% compute list of  of indices corresponding to the positions in the grid of each point ranked by bathymetry value
tic=time()
idx_list=idx_table(FineBath)
toc=time()
print(toc-tic)

##%% save output
np.save('Jessie/grids_extGreenLaw/MyRes/idx_Bandol_part1', idx_list)

#%% load output if needed

idx_list=np.load('Jessie/grids_extGreenLaw/MyRes/idx_Bandol_part1.npy','r')


#%% compute forecast and betas over the fine grid

tic=time()
Betas, test =fine_forecast(Computed_heights, Bath, FineSimRes, FineBath, lon, lat, finelon, finelat, idx_list)
toc=time()
print(toc-tic)

#%% save outputs
np.save('Jessie/grids_extGreenLaw/MyRes/Bandol_part1Forecast',test)

np.save('Jessie/grids_extGreenLaw/MyRes/Bandol_part1Betas',Betas)
