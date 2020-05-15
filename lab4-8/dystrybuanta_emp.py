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

def real_dyst(x_d):
    dystr=[]
    for x in x_d:
        #y=0.5*(1+special.erf((x-1)/(1*np.sqrt(2))))
        y=1/np.pi*np.arctan(x-1)+0.5
        dystr.append(y)
    return dystr

def real_dyst_uni(x_d):
    dystr=[]
    for x in x_d:
        if(x<0):
            y=0
        elif(x>=0 and x<=1):
            y=x
        elif(x>1):
            y=1
        dystr.append(y)
    return dystr
        
def MSE_err(N,L,x_d):
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
            emp_dist=dyst_estymator(x_d,testtab)
            real_dist=real_dyst(x_d)
            for j in range(0,len(real_dist)):
                mse_sum+=(emp_dist[j]-real_dist[j])**2
            uni_tab.clear()
            testtab.clear()
            emp_dist.clear()
            real_dist.clear()
        err=mse_sum/L
        mse_errtab.append(err)
    return mse_errtab

def MSe_err(L,N,x_d):
    mse_errtab=[]
    uni_tab=[]
    testtab=[]
    emp_dist=[]
    real_dist=[]
    for l in L:
        mse_sum=0
        for i in range(0,l): 
            gen_sawtooth(uni_tab,N,88000)
            for x in uni_tab:
                y=roz_normalny(x,0,1)
                testtab.append(y)
            emp_dist=dyst_estymator(x_d,testtab)
            real_dist=real_dyst(x_d)
            for j in range(0,len(real_dist)):
                mse_sum+=(emp_dist[j]-real_dist[j])**2                
            uni_tab.clear()
            testtab.clear()
            emp_dist.clear()
            real_dist.clear()
        err=mse_sum/l
        mse_errtab.append(err)
    return mse_errtab

def MSE_err_uni(N,L,x_d):
    mse_errtab=[]
    uni_tab=[] #tablica z roz. jednostajnym do tworzenia kolejnych rozkladow
    testtab=[] #tablica z rozkładem do estymacji 
    emp_dist=[]
    real_dist=[]
    for n in N:
        mse_sum=0
        for i in range(0,L):
            gen_sawtooth(testtab,n,88000)
            emp_dist=dyst_estymator(x_d,testtab)
            real_dist=real_dyst_uni(x_d)
            for j in range(0,len(real_dist)):
                mse_sum+=(emp_dist[j]-real_dist[j])**2
            testtab.clear()
            emp_dist.clear()
            real_dist.clear()
        err=mse_sum/L
        mse_errtab.append(err)
    return mse_errtab

if __name__=='__main__': 
    filepath="ModLab5Dat2.txt"
    dat=read_data(filepath)
    Xn=[]
    X_pocz=[]
    x_d=np.linspace(-5, 5, 10000)
    X_pocz2=[]
    gen_sawtooth(X_pocz,10000,88000)
    for x in X_pocz:
        y=roz_normalny(x,0,1)
        Xn.append(y)
   
    

    fig1=plt.plot(x_d,dyst_estymator(x_d,dat),label='Dystrybuanta empiryczna')
    fig2=plt.plot(x_d,real_dyst(x_d),label='Rozklad Cauchyego x0=1, gamma=1')
    plt.legend(loc="upper left")
    plt.title("Dystrybuanta empiryczna ModLab5Dat2.txt")
    plt.show()    
    """
    N=[]
    for i in range(10,1000,40):
        N.append(i)
    L=10
    res=MSE_err_uni(N,L,x_d)
    plt.plot(N,res)
    plt.title("Błąd empiryczny MSE przy L=10")
    plt.xlabel("liczba próbek N")
    plt.ylabel("wartość błędu")
    plt.show()
    """
