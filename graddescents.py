# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 12:28:08 2019

@author: Jessie
"""
import numpy as np
import numpy.random as rnd
import matplotlib.pyplot as plt


def fixedstep(x0,f,gradf,h,epsilon):
    '''
    a basic fixed-step gradient descent
    x0 starting point
    f and gradf : function and derivative of the function
    h : step
    epsilon : error
    '''
    x0=np.asarray(x0)
    X=[x0]
    n=0
    while np.linalg.norm(gradf(x0))>epsilon and n<10**5:
        
        x0=x0-h*gradf(x0)
        n=n+1
        X.append(x0)
        
    return X   


def gradnorm(x0,f,gradf,h,epsilon):
    '''
    normalized gradient descent
    x0 starting point
    f and gradf : function and derivative of the function
    h : step
    epsilon : error
    '''
    x0=np.asarray(x0)
    X=[x0]
    n=0
    while np.linalg.norm(f(x0))>epsilon and n<10**4:
        a=np.linalg.norm(gradf(x0))
        if a>0:
            x0=x0-h*gradf(x0)/a
        n=n+1
        X.append(x0)
    return X    

def wolfe(q,gradq,m1,m2,a,t0):
    '''
    Wolfe and Armijo conditions

    '''
    t=t0
    tg=0
    td=0
    while q(t)>q(0)+t*m1*gradq(0) or gradq(t)<m2*gradq(0):
        if q(t)>q(0)+t*m1*gradq(0):
            td=t
        else:
            if gradq(t)<m2*gradq(0):
                tg=t
        if td==0:
            t=a*t
        else :
            t=(tg+td)/2
    return t    

def gradwolfe(x0,f,gradf,epsilon,m1,m2,a,t0):
    '''
    gradient descent with step determined by Wolfe and Armijo conditions
    x0 starting point
    f and gradf : function and derivative of the function
    epsilon : error
    '''
    x0=np.asarray(x0)
    X=[x0]
    n=0
    while np.linalg.norm(gradf(x0))>epsilon and n<10**5:
        dx=-gradf(x0)
        q=lambda t: np.linalg.norm(f(x0+t*dx))
        dq= lambda t: np.sum(dx*gradf(x0+t*dx))
        h=wolfe(q,dq, m1,m2,a,t0)
        x0=x0+h*dx
        n=n+1
        X.append(x0)
    return X    

