# -*- coding: utf-8 -*-
"""
Created on Sat May  2 17:37:19 2020

@author: Hubert
"""
from generatory import *
from estymatory import *
import random
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import signal
from scipy import special
import scipy.stats as st
from statsmodels.distributions.empirical_distribution import ECDF



def est_gest_hist(x_d,Xn):
    emp=[]
    alfa=0.5
    n=len(Xn)
    N=np.sqrt(n)
    h=math.pow(N,-alfa)
    for xi in x_d:
        suma=0
        for XN in Xn:
            if(XN>xi and XN<=xi+h):
                suma+=1
            else:
                pass
        q=suma/(N*h)
        emp.append(q)
    return emp


def gaussian_kernel(x):
    return 1/(np.sqrt(2*np.pi))*math.exp((-1*math.pow(x,2))/2)

def epanechnikov_kernel(x):
    if(np.abs(x)<=1):
        return 3/4*(1-math.pow(x,2))
    else:
        return 0

def boxcar_kernel(x):
    if(np.abs(x)<=1):
        return 0.5
    else:
        return 0



def est_jadrowy(x_d,Xn,kernel_func):
    emp=[]
    alfa=0.5
    n=len(Xn)
    N=np.sqrt(n)
    h=math.pow(N,-alfa)
    for xi in x_d:
        suma=0
        for XN in Xn:
            suma+=kernel_func((XN-xi)/h)
        q=suma/(n*h)
        emp.append(q)
    return emp
    


if __name__=='__main__': 
    X_pocz=[]
    Xn=[]
    X_test=[]
    x_d=np.linspace(-5, 5, 2000)
    gen_sawtooth(X_pocz,2000,88000)
    for x in X_pocz:
        y=roz_normalny(x,0,1)
        Xn.append(y)
        
    for x in range(0,1000):
        X_test.append(random.uniform(0,1))
        
    plt.plot(x_d,est_jadrowy(x_d,Xn,epanechnikov_kernel))
    plt.title('Estymator gestosci')
    plt.show()