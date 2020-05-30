# -*- coding: utf-8 -*-
"""
Created on Sat May 30 17:11:29 2020

@author: Hubert
"""
import itertools
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

N=100
La=10
Lb=10
k=1

def fi_vec(i_t,un,yn):
    y_temp=[]
    fi=[]
    k=2
    for i in range(Lb):
        fi.append(un[i_t-k-i])
        y_temp.append(-yn[i_t-i])
    for y in y_temp:
        fi.append(y)
    return fi



if __name__=="__main__":
    Theta=[1,2,3,4,5,6,7,8,9,10,21,22,23,24,25,26,27,28,29,30]
    #Un
    un=np.zeros(N)
    for i in range(N):
        un[i]=np.random.normal(0,1)
    #en
    en=np.zeros(N) 
    for i in range(N):
        en[i]=np.random.normal(0,2)
    #vn
    v_temp=[]
    a_temp=[]
    b_temp=[]
    u_temp=[]
    vn=np.zeros(N)
    
    for i in range(10):
        vn[i]=np.random.normal(0,1)
    for i in range(Lb,Lb+La):
        a_temp.append(Theta[i])
    for i in range(Lb):
        b_temp.append(Theta[i])         
    for i in range(La,N):
        for j in range(i-1,i-11,-1):
            v_temp.append(vn[j])
        for k in range(i,i-10,-1):
            u_temp.append(vn[j])
        for l in range(10):
            vn[i]+=a_temp[l]*v_temp[l]+b_temp[l]*u_temp[l]
        v_temp.clear()
        u_temp.clear()
    #yn
    yn=np.zeros(N)
    for i in range(N):
        yn[i]=vn[i]+en[i]
    
    #fi
    Fi=[]
    for i in range(1,N):
        fi_temp=fi_vec(i,un,yn)
        Fi.append(fi_temp)
    #Ri
    Ri=[]
    alfa=500
    lambd=2
    Ri.append(np.eye(10)*alfa)
    for i in range(1,N):
        Ri_temp=Ri[i-1]-((np.dot(np.dot(np.dot(Ri[i-1],Fi[i-1]),Fi[i-1].transpose()),Ri[i-1]))/(np.add(lambd,np.dot(np.dot(Fi[i-1].transpose(),Ri[i-1]),Fi[i-1])) ))
        print(1/lambd*Ri_temp)
        Ri.append(Ri_temp)
    
    
 