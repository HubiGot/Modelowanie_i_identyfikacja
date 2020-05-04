# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 13:45:08 2020

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


def read_data(filepath):
    data=[]
    f=open(filepath,"r")
    for line in f.readlines():
        data.append(float(line))
    return data


def dyst_estymator(x_d,Xn):
    emp=[]
    for xi in x_d:
        suma=0
        for XN in Xn:
            if(XN<=xi):
                suma+=1
            else:
                pass
        q=suma/len(Xn)
        emp.append(q)
    return emp


if __name__=='__main__':
    filepath="ModLab5Dat1.txt"
    dat=read_data(filepath)
    Xn=[]
    X_pocz=[]
    x_d=np.linspace(-5, 5, 10000)
    gen_sawtooth(X_pocz,10000,88000)
    for x in X_pocz:
        y=roz_normalny(x,1,5)
        Xn.append(y)

    fig1=plt.plot(x_d,dyst_estymator(x_d,dat),label='dane')
    fig2=plt.plot(x_d,dyst_estymator(x_d,Xn),label='cauchy')
    plt.legend(loc="upper left")
    plt.title("Dystrybuanty empiryczne")
    plt.show()    


    
