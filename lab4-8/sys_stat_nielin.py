# -*- coding: utf-8 -*-
"""
Created on Sun May  3 15:30:29 2020

@author: Hubert
"""

from generatory import *
from estymatory import *
from est_gestosc import *
import random
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import signal
from scipy import special
import scipy.stats as st
from statsmodels.distributions.empirical_distribution import ECDF


def est_nadaray_watson(x_d,Xn,Yn,kernel_func):
    result=[]
    alfa=0.5
    n=len(Xn)
    N=np.sqrt(n)
    #h=math.pow(N,-alfa)
    h=0.2
    for xi in x_d:
        sumaYK=0
        sumaK=0
        for xn,yn in zip(Xn,Yn):
            sumaYK+=yn*kernel_func((xn-xi)/h)
            sumaK+=kernel_func((xn-xi)/h)
        q=sumaYK/sumaK
        result.append(q)
    return result
        
    


if __name__ =='__main__':
    X_pocz=[]
    Xn=[]
    Zn=[]
    Yn=[]
    alfa=2
    x_d=np.linspace(-2, 2, 1000)
    gen_sawtooth(X_pocz, 1000, 88000)
    for x in X_pocz:
        tmp=random.uniform(0,1)
        if(tmp>=0.5):
            y=2*x
            Xn.append(y)
        else:
            y=-2*x
            Xn.append(y)            
   
    for x in Xn:
        z=np.random.normal(0, 1)
        yn=math.atan(2*x)+z
        Yn.append(yn)
        
    plt.plot(x_d,est_nadaray_watson(x_d,Xn,Yn,epanechnikov_kernel))
    plt.title("Gaussian kernel")
    plt.show()
    plt.scatter(Xn,Yn)
    plt.title("Estymacja systemu nieliniowego")
    plt.show()
    

    