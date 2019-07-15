from netCDF4 import Dataset
import numpy as np
import os
import scipy.interpolate as scp
import scipy
from time import time

from graddescents import *


def idx_table(Bathy_array):
    '''
    Bathy_array : table of bathymetry values at the corresponding points in the grid
    returns a list of indices corresponding to the positions in the grid of each point ranked by bathymetry value
    '''
    east,north=np.shape(Bathy_array)
    ordBath=list(np.reshape(Bathy_array.copy(),(east*north,1)))
    ordBath.sort()
    # put bathymetry table in a list and sort it
    
    idx_list=np.zeros((east*north,2))
    
    for i in range (len(ordBath)): # find corresponding indices
        a,b=np.where(np.abs(Bathy_array-ordBath[i])==np.amin(np.abs(Bathy_array-ordBath[i])))
        for u,v in zip(a,b):
                idx_list[i,0]=u
                idx_list[i,1]=v

    return idx_list            



#%% coarse grid to fine grid

def interpolate(x,y,xfine,yfine,Z,jmin,imin):
    '''
    interpolates coarse grid onto fine grid on one point
    imin, jmin : coordinates of the point of interest
    x,y : small coarse grids of latitude and longitude containing (imin,jmin)
    z : small coarse grid of simulation results containing (imin, jmin)
    xfine, yfine : small fine grids of latitude and longitude containing (imin,jmin)

    '''
    # Regularly-spaced, coarse grid
    X=x.copy()
    Y=y.copy()
    print(x,y)
    interp_spline = scp.RectBivariateSpline(x,y, Z,kx=2,ky=2)

    # Regularly-spaced, fine grid
    dx2, dy2 = np.abs(xfine[1]-xfine[0]),np.abs(yfine[1]-yfine[0])
    x2 = np.arange(np.amin(X), np.amax(X), dx2)
    y2 = np.arange(np.amin(Y), np.amax(Y), dy2)

    X2, Y2 = np.meshgrid(x2,y2)
    Z2 = interp_spline(x2,y2)

    #find closest point in the real fine grid
    umin=np.where(np.abs(x2-xfine[imin])==np.amin(np.abs(x2-xfine[imin])))[0]
    vmin=np.where(np.abs(y2-yfine[jmin])==np.amin(np.abs(y2-yfine[jmin])))[0]
    print(umin,vmin)
    print(Z2[umin,vmin])
    return Z2[umin,vmin]
    
def coarse_to_fine(SimRes, Bath, FineSimRes, FineBath, Lat, Lon, FineLat, FineLon, lowlim=-50, highlim=0, hE=50):

     
    '''
    calculates values around 50m depth, from the coarse grid to the fine one
    
    '''
   
    (east, north)=np.shape(Bath)
    (fine_east, fine_north,N)=np.shape(FineSimRes)
    Fine_computed_heights2=np.zeros((fine_east,fine_north,N))

    '''
    Calculates the height at the minimum
    Interpolation with adjacent coarse grid cells
    '''

    imin,jmin=np.where(FineBath==np.amin(FineBath))
    imin,jmin=imin[0],jmin[0]
    print(imin, jmin)

    idx=np.where(FineBath<=0)
    Fine_computed_heights2[idx[0],idx[1],:]=np.nan
        
    '''
    get closest big cell to small cell (imin, jmin)
    '''
    u,v= np.argmin(np.abs(Lat-FineLat[imin])),  np.argmin(np.abs(Lon-FineLon[jmin]))
    
    w=2
    q=2
    z=2
    while q<north and np.count_nonzero(np.isnan(SimRes[u,np.mod(v-q,north)]))>0 :
        q=q+1
    while w<east and np.count_nonzero(np.isnan(SimRes[np.mod(u+w,east),v]))>0 :
        w=w+1
    while z<east and np.count_nonzero(np.isnan(SimRes[np.mod(u-z,east),v]))>0 :
        z=z+1
    

    numb=max(q,z,w)
    '''
    small coarse latitude and longitude tables over deepest point
    '''
    X=Lat[min(np.mod(u-numb,east),np.mod(u+numb,east)):max(np.mod(u-numb,east),np.mod(u+numb,east))]
    Y=Lon[min(np.mod(v-numb,north),np.mod(v+numb,north)):max(np.mod(v-numb,north),np.mod(v+numb,north))]   

    for n in range(N):
        
        '''
        small grid of simulation results from precalculated data
        '''
        Z=SimRes[min(np.mod(u-numb,east),np.mod(u+numb,east)):max(np.mod(u-numb,east),np.mod(u+numb,east)),min(np.mod(v-numb,north),np.mod(v+numb,north)):max(np.mod(v-numb,north),np.mod(v+numb,north)),n]
        Fine_computed_heights2[imin, jmin, n]=interpolate(X,Y, FineLon,FineLat,Z,imin,jmin)
        
    print(Fine_computed_heights2[imin,jmin,:],'interpolated')
    print(FineSimRes[imin,jmin,:],'finesimres')
    print(SimRes[np.mod(u-w,east),v,:],'coarse res')
    
    return imin, jmin, Fine_computed_heights2
    
#%% forecast on fine grid


def fine_forecast(SimRes, Bath, FineSimRes, FineBath, Lat, Lon, FineLat, FineLon, idx_list, lowlim=-50, highlim=-5, hE=50):
    
     '''
     forecast all values below 50m in the fine grid using Green's law
     SimRes, FineSimRes : simulation results for coarse and fine grids
     Bath, FineBath : bathymetry tables for coarse and fine grids
     Lat, Lon, FineLat, FineLon : latitude and longitude tables
     idx_list : list of  of indices corresponding to the positions in the grid of each point ranked by bathymetry value from idx_table(FineBath)
     lowlim, highlim : range for Extended Green's Law
     hE : typical value 
     '''
     
     #create all the tables
     
     (east, north)=np.shape(Bath) #i is actually north and j east, it's inverted
     (fine_east, fine_north,N)=np.shape(FineSimRes)
     Betas = np.zeros((fine_east, fine_north))
     beta_zero=0.5
     print(fine_east, fine_north)
     print(np.shape(FineLat),np.shape(FineLon))
     
     imin, jmin, Fine_computed_heights= coarse_to_fine(SimRes, Bath, FineSimRes, FineBath, Lat, Lon, FineLat, FineLon, lowlim, highlim, hE)
     print(Fine_computed_heights[imin,jmin,:],FineSimRes[imin,jmin,:])
   
     print(fine_north*fine_east,'size')

     x,y=0,0
     u,v=np.where(FineBath==np.amin(FineBath))
     u,v=u[0],v[0]
     for i in range(1, fine_east*fine_north):
         '''
         go over the whole grid
         '''

             x,y=u,v
             u,v=int(idx_list[i,0]),int(idx_list[i,1])
             
             if FineBath[u,v]<-50:
                 '''
                 Use Geen's law
                 '''
                 for n in range (N):
                    Fine_computed_heights[u,v,n]= Fine_computed_heights[x,y,n]*(FineBath[x,y]/FineBath[u,v])**(1/4)
    
             
             if 0>=FineBath[u,v]>=-50: 
                     '''
                     use extended Green's law (calculate beta using a gradient descent) while bathymetry is strictly positive (before beta becomes too big,
                      typical beta values rang from -0.5 to 1.5)
                     '''
                     betatest=computeB(hE, FineBath[x,y], FineBath[u,v], Fine_computed_heights[x,y,:], FineSimRes[u,v,:], descent=fixedstep, epsilon = 10**-4, beta0=beta_zero, h=10**-2)[-1]
                     if np.isnan(betatest)==False and np.abs(betatest)<3:
                         Betas[u,v]=betatest
                         beta_zero=Betas[u,v]
                         for n in range (N):
                            Fine_computed_heights[u,v,n]= Fine_computed_heights[x,y,n]*(1+Betas[u,v]*((hE-FineBath[u,v])/hE))*(FineBath[x,y]/FineBath[u,v])**(1/4)
                        
             
                         
             if np.mod(i,10000)==0:
                 print(i)
                 print(Betas[u,v], "beta")
                 print(FineBath[u,v], 'bathy', FineBath[x,y], 'bathy source')
                 print(Fine_computed_heights[u,v,:], "computed")
                 print(Fine_computed_heights[x,y,:], 'computed source')
                 print(FineSimRes[u,v,:], "simres")
                 np.save('Jessie/grids_extGreenLaw/MyRes/NiceAirForecast',Fine_computed_heights)
                 np.save('Jessie/grids_extGreenLaw/MyRes/NiceAirBetas',Betas)
     
     
                 
     return Betas, Fine_computed_heights


    
#%% beta parameter

def computeB(hE, h1, h2, eta1, eta2, descent = fixedstep, epsilon = 10**-3, beta0=0.5, h=10**-2): 
    '''
    Compute one beta using a fixed step gradient descent 
    Cost function is MSE
    eta 2 : simulation result for the wave height
    etap : wave height to be predicted
    eta 1 : last predicted wave height
    '''
    #tic=time()
    eta1=np.asarray(eta1)
    eta2=np.asarray(eta2)
    
    #find how many tsunami simulations are inputted
    try :
        (N,)=np.shape(eta1)
    except :
        N=1
        
    # compute cost function and gradient of cost function
    if h2!=hE:    
        etap = lambda beta : (h1/h2)**(1/4)*(1+beta*(hE-h2)/hE)*eta1
        Cost=lambda beta : np.vdot(eta2-etap(beta), eta2-etap(beta))*(1/N)
        Cost_grad = lambda beta : -2/N*((hE-h2)/hE)*(h1/h2)**(1/4)*np.vdot(eta1, eta2-etap(beta))
    else : 
        hE=150
        etap = lambda beta : (h1/h2)**(1/4)*(1+beta*(hE-h2)/hE)*eta1
        Cost=lambda beta : np.vdot(eta2-etap(beta), eta2-etap(beta))*(1/N)
        Cost_grad = lambda beta : -2/N*((hE-h2)/hE)*(h1/h2)**(1/4)*np.vdot(eta1, eta2-etap(beta))
        
    #use a gradient descent from graddescents
    beta = descent(beta0 ,Cost, Cost_grad,h, epsilon)
    return beta   
    
