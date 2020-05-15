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

def base_cos(x,i):
    if(i==0):
        return np.sqrt(1/(2*np.pi))
    else:
        return np.sqrt(1/np.pi)*math.cos(i*x)


def alfa_i(Xn,i):
    fi_sum=0
    for x in Xn:
        fi_sum+=base_cos(x,i)
    return fi_sum/len(Xn)

def beta_i(Xn,Yn,i):
    fi_sum=0
    for x,y in zip(Xn,Yn):
        fi_sum+=y*base_cos(x,i)
    return fi_sum/len(Xn)
    
def est_mx(Xn,Yn,L):
    result_tab=[]
    for xd in x_d:
        LSum=0
        MSum=0
        for i in range(0,L):
            LSum+=beta_i(Xn,Yn,i)*base_cos(xd,i)
            MSum+=alfa_i(Xn,i)*base_cos(xd,i)
        mx=LSum/MSum
        result_tab.append(mx)
    return result_tab

def MSE_err(N,L,x_d):
    mse_errtab=[]
    X_pocz=[]
    Xn=[]
    Yn=[]
    real_funtab=[]
    est_funtab=[]
    for l in L:
        mse_sum=0
        gen_sawtooth(X_pocz, N, 88000)
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
        est_funtab=est_mx(Xn,Yn,l)
        real_funtab=calculate_fun_tab(Xn,2)
        for j in range(0,len(real_funtab)):
            mse_sum+=(est_funtab[j]-real_funtab[j])**2
        err=mse_sum/N
        mse_errtab.append(err)
        Xn.clear()
        Yn.clear()
        X_pocz.clear()
        real_funtab.clear()
        est_funtab.clear()
    return mse_errtab
            
            
            
if __name__=='__main__':
    """
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
    
    plt.plot(x_d,est_mx(Xn,Yn,5)) 
    plt.title('Estymacja charakterystyki systemu')
    plt.show()
    """

    L=[]
    for i in range(5,50,5):
        L.append(i)
    N=1000
    x_d=np.linspace(-3,3,N)
    plt.plot(L,MSE_err(N,L,x_d))
    plt.title("Błąd empiryczny MSE")
    plt.xlabel("liczba funkcji bazowych L")
    plt.ylabel("wartość błędu")
    plt.show()
    
    """
    plt.plot(x_d,calculate_fun_tab(x_d,2),label='nieliniowa charakterystyka systemu',color='red')
    plt.scatter(Xn,Yn, label='chmura pomiarów')
    plt.legend(loc="upper left")
    plt.title("Nieliniowa charakterystyka z chmurą pomiarów")
    plt.show()
    """