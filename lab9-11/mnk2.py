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
s=0.5

def mnk_est(Xn,Yn):
    an=np.dot(np.dot(np.linalg.inv(np.dot(Xn.transpose(),Xn)),Xn.transpose()),Yn)
    return an
    
def cov_an(Xn,delta_z):
    cov_est=np.dot(delta_z**2,np.linalg.inv(np.dot(Xn.transpose(),Xn)))
    return cov_est

def err(Xn,tab_N,a):
    err_tab=[]
    Zn=np.zeros((K,N),dtype=float)
    Ksi_n=np.zeros((K,N),dtype=float)
    for k in range(K):
        for j in range(N):
            Ksi_n[k][j]=np.random.normal(0,5)
    for k in range(K):
        Zn[k][0]=Ksi_n[k][0];
        for j in range(1,N):
            Zn[k][j]=Ksi_n[k][j]+s*Ksi_n[k][j-1]
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
    #wektor a
    a=[1,2,3,4,5,6,7,8,9,10]
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
    Ksi_n=np.zeros(N)
    for i in range(N):
        Ksi_n[i]=np.random.normal(0,np.sqrt(var))
    Zn[0]=Ksi_n[0]
    for i in range(1,N):
        Zn[i]=Ksi_n[i]+s*Ksi_n[i-1]
    #Yn
    Yn=np.add(np.dot(Xn,a),Zn)
    a_est=mnk_est(Xn,Yn)
    print(a_est)
    print("")
    
    #macierz kowariancji R wektora Zn
    R=np.zeros((N,N),dtype=float)
    for i in range(N):
        for j in range(N):
            if i==j:
                R[i][j] = (1+pow(s,2))*var
            elif np.absolute(i-j)==1:
                R[i][j] = s*var
    
    #kowariancja estymatora
    cov_estymatora = np.dot(np.dot(np.dot(np.dot(np.linalg.inv(np.dot(Xn.transpose(),Xn)),Xn.transpose()), R),Xn),np.linalg.inv(np.dot(Xn.transpose(),Xn)))           
    """
    cov_a=cov_an(Xn,np.sqrt(var))
    """
    # Set up grid and data
    nx, ny = 10, 10
    x = range(nx)
    y = range(ny)

    data = cov_estymatora
    hf = plt.figure()
    ha = hf.add_subplot(111, projection='3d')

    X, Y = np.meshgrid(x, y)
    ha.plot_surface(X, Y, data,rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
    ha.set_title("Macierz kowariancji estymatora an")
    plt.show()
    
    tab_N=np.linspace(15,100,50).astype(int)
    error_tab=err(Xn,tab_N,a)
    plt.plot(tab_N,error_tab)
    plt.title('Zależność wartości błędu od liczby próbek N')
    plt.ylabel('wartosc bledu')
    plt.xlabel('liczba probek N')    
    
    
    

