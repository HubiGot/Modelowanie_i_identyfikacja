# -*- coding: utf-8 -*-
"""
Created on Fri May 29 09:32:54 2020

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
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D 

L=10
K=20
N=100
var=1
delta_x=1

def mnk_est(Xn,Yn):
    an=np.dot(np.dot(np.linalg.inv(np.dot(Xn.transpose(),Xn)),Xn.transpose()),Yn)
    return an
    
def cov_an(Xn,delta_z):
    cov_est=np.dot(delta_z**2,np.linalg.inv(np.dot(Xn.transpose(),Xn)))
    return cov_est

def err(Xn,tab_N,a):
    err_tab=[]
    Zn=np.zeros((K,N),dtype=float)
    for k in range(K):
        for j in range(N):
            Zn[k][j]=np.random.normal(0,2)
    for n in tab_N:
        x_tmp=Xn[:n,:]
        err=0
        for k in range(K):
            z_tmp=Zn[k][:n]
            y_tmp=np.add(np.dot(x_tmp,a),z_tmp)
            a_est=mnk_est(x_tmp,y_tmp)
            err+= 1/K*pow(np.linalg.norm(a_est-a),2)
        err_tab.append(err)
    return err_tab

if __name__=="__main__":
    a=[2,2,2,2,2,2,2,2,2,2]
    #macierz Xn
    Xn=np.zeros((N,L),dtype=float)
    I=np.eye(L)
    Sigma=I*delta_x
    mi=np.zeros(L)
    for i in range(N):
        for j in range(L):
            Xn[i][j]=np.random.normal(0,1)            
    #zaklocenia Z
    Zn=np.zeros(N)
    for i in range(N):
        Zn[i]=np.random.normal(0,np.sqrt(var))
    #Yn
    Yn=np.add(np.dot(Xn,a),Zn)
    a_est=mnk_est(Xn,Yn)
    print(a_est)
    print("")
    cov_a=cov_an(Xn,np.sqrt(var))
    
    nx, ny = 10, 10
    x = range(nx)
    y = range(ny)

    data = cov_a
    hf = plt.figure()
    ha = hf.add_subplot(111, projection='3d')

    X, Y = np.meshgrid(x, y)
    ha.plot_surface(X, Y, data,rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
    ha.set_title("Macierz kowariancji estymatora")
    plt.show()
    
    tab_N=np.linspace(15,100,50).astype(int)
    error_tab=err(Xn,tab_N,a)
    plt.plot(tab_N,error_tab)
    plt.title('Zależność wartości błędu od liczby próbek N')
    plt.ylabel('wartosc bledu')
    plt.xlabel('liczba probek N')    
    
    
    