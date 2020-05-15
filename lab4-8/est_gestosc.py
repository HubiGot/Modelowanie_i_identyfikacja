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
  
def real_gest_fn(x_d):
    dystr=[]
    for x in x_d:
        #y=0.5*(1+special.erf((x-1)/(1*np.sqrt(2))))
        y=1/(1*np.sqrt(2*np.pi))*np.exp(-(x-0)**2/2)
        dystr.append(y)
    return dystr

def MSE_err(N,L,x_d):
    mse_errtab=[]
    uni_tab=[] #tablica z roz. jednostajnym do tworzenia kolejnych rozkladow
    testtab=[] #tablica z rozkładem do estymacji 
    emp_gest=[]
    real_gest=[]
    for n in N:
        mse_sum=0
        for i in range(0,L):
            gen_sawtooth(uni_tab,n,88000)
            for x in uni_tab:
                y=roz_normalny(x,0,1)
                testtab.append(y)
            emp_gest=est_gest_hist(x_d,testtab)
            real_gest=real_gest_fn(x_d)
            for j in range(0,len(real_gest)):
                mse_sum+=(emp_gest[j]-real_gest[j])**2
            uni_tab.clear()
            testtab.clear()
            emp_gest.clear()
            real_gest.clear()
        err=mse_sum/L
        mse_errtab.append(err)
    return mse_errtab

def MSE_err_jadr(N,L,x_d):
    mse_errtab=[]
    uni_tab=[] #tablica z roz. jednostajnym do tworzenia kolejnych rozkladow
    testtab=[] #tablica z rozkładem do estymacji 
    emp_dist=[]
    real_dist=[]
    for n in N:
        mse_sum=0
        for i in range(0,L):
            gen_sawtooth(uni_tab,n,88000)
            for x in uni_tab:
                y=roz_normalny(x,0,1)
                testtab.append(y)
            emp_dist=est_gest_hist(x_d,testtab)
            real_dist=real_gest_fn(x_d)
            for j in range(0,len(real_dist)):
                mse_sum+=(emp_dist[j]-real_dist[j])**2
            uni_tab.clear()
            testtab.clear()
            emp_dist.clear()
            real_dist.clear()
        err=mse_sum/L
        mse_errtab.append(err)
    return mse_errtab

if __name__=='__main__': 
    """
    X_pocz=[]
    Xn=[]
    X_test=[]
    x_d=np.linspace(-2, 2, 3000)
    gen_sawtooth(X_pocz,3000,88000)
    for x in X_pocz:
        #y=roz_normalny(x,0,1)
        y=roz_trojkatny2(x)
        Xn.append(y)
    
    plt.plot(x_d,est_gest_hist(x_d,Xn))
    plt.title("Estymacja gęstości w oparciu o histogram")
    plt.show()
    """
    x_d=np.linspace(-5, 5, 2000)
    N=[]
    for i in range(10,500,40):
        N.append(i)
    L=10
    #res=MSE_err(N,L,x_d)
    plt.plot(N,MSE_err_jadr(N,L,x_d))
    plt.title("Błąd empiryczny MSE przy L=10")
    plt.xlabel("liczba próbek N")
    plt.ylabel("wartość błędu")
    plt.show()
    
    """
    plt.plot(x_d,est_jadrowy(x_d,Xn,boxcar_kernel))
    plt.title("Estymator jądrowy gęstości - Boxcar kernel")
    plt.show()
    
      
    for x in range(0,1000):
        X_test.append(random.uniform(0,1))      
    plt.plot(x_d,est_jadrowy(x_d,Xn,epanechnikov_kernel))
    plt.title('Estymator gestosci')
    plt.show()
    """