# -*- coding: utf-8 -*-
"""
Created on Mon May  4 13:50:24 2020

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

def calculate_fun(x,alfa):
    if(np.abs(x)>=0 and np.abs(x)<1):
        return alfa*math.pow(x,2)
    elif(np.abs(x)>=1 and np.abs(x)<2):
       return 1
    elif(np.abs(x)>=2):
       return 0


def calculate_fun_tab(Xn,alfa):
    mx=[]
    for x in Xn:
        if(np.abs(x)>=0 and np.abs(x)<1):
            yn=alfa*math.pow(x,2)
            mx.append(yn)
        elif(np.abs(x)>=1 and np.abs(x)<2):
            mx.append(1)
        elif(np.abs(x)>=2):
            mx.append(0)
    return mx
            
            
if __name__=='__main__':
    X_pocz=[]
    Xn=[]
    Yn=[]
    gen_sawtooth(X_pocz, 1000, 88000)
    x_d=np.linspace(-3,3,1000)
    for x in X_pocz:
        tmp=random.uniform(0,1)
        if(tmp>=0.5):
            y=np.pi*x
            Xn.append(y)
        else:
            y=-np.pi*x
            Xn.append(y)  
    for x in Xn:
        rng=np.random.normal(0,1)
        yn=calculate_fun(x,2)+rng
        Yn.append(yn)
    plt.plot(x_d,calculate_fun_tab(x_d,2),color='green')
    plt.scatter(Xn,Yn)
    plt.title("Chmura nieliniowego systemu")
    plt.show()