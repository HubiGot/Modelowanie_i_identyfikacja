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
import subprocess as sp



def est_nadaray_watson(x_d,Xn,Yn,kernel_func):
    result=[]
    alfa=0.5
    #n=len(Xn)
    #N=np.sqrt(n)
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
        
def real_char(x_d):
    atan_tab=[]
    for x in x_d:
        atan_tab.append(math.atan(2*x))
    return atan_tab




def MSE_err2(N,L,x_d):
    mse_errtab=[]
    uni_tab=[] #tablica z roz. jednostajnym do tworzenia kolejnych rozkladow
    emp_dist=[]
    real_dist=[]
    X3=[]
    Y3=[]
    for n in N:
        sumL=0
        for i in range(0,L):
            mse_sum=0
            gen_sawtooth(uni_tab,n,88000)
            for x in uni_tab:
                tmp=np.random.uniform(0,1)
                if(tmp>=0.5):
                    y=2*x
                    X3.append(y)
                else:
                    y=-2*x
                    X3.append(y)
            for x in X3:
                z=np.random.normal(0,1)
                yn=math.atan(2*x)+z
                Y3.append(yn)
            emp_dist=est_nadaray_watson(x_d,X3,Y3,epanechnikov_kernel)
            real_dist=real_char(x_d)
            for j in range(0,len(real_dist)):
                mse_sum+=(emp_dist[j]-real_dist[j])**2
            mse_sum=mse_sum/n
            sumL+=mse_sum
            uni_tab.clear()
            emp_dist.clear()
            real_dist.clear()
            X3.clear()
            Y3.clear()
        err=sumL/L
        mse_errtab.append(err)
    return mse_errtab;


if __name__ =='__main__':
    
    X_pocz=[]
    Xn=[]
    Zn=[]
    Yn=[]
    atan_tab=[]
    alfa=[-2,-1,1,2]
    h_tab=[0.1, 0.5, 0.8, 1]
    x_d=np.linspace(-2, 2, 500)
    gen_sawtooth(X_pocz, 500, 88000)
    for x in X_pocz:
        tmp=random.uniform(0,1)
        if(tmp>=0.5):
            y=2*x
            Xn.append(y)
        else:
            y=-2*x
            Xn.append(y)            
    
    
    
    N=[]
    for i in range(100,1000,50):
        N.append(i)
    L=10

    plt.plot(N,MSE_err2(N,L,x_d))
    plt.title("Błąd empiryczny MSE przy L=10")
    plt.xlabel("liczba próbek N")
    plt.ylabel("wartość błędu")
    plt.show()
    
    
    
    """
    for x in Xn:
        z=np.random.normal(0, 1)
        yn=math.atan(2*x)+z
        Yn.append(yn)
     
    
    for x in x_d:
        atan_tab.append(math.atan(2*x))
    
    plt.plot(x_d,est_nadaray_watson(x_d,Xn,Yn,epanechnikov_kernel))
    plt.title("Estymacja charakterystyki systemu - Epanechnikov kernel")
    plt.show()
    """
    """
    for a in alfa:
        for x in Xn:
            z=np.random.normal(0, 1)
            yn=math.atan(a*x)+z
            Yn.append(yn)
        plt.plot(x_d,est_nadaray_watson(x_d,Xn,Yn,epanechnikov_kernel),label=' $a = %s $' % a)
        Yn.clear()
    plt.title("Kształt nieliniowosci w systemie" )
    plt.legend(loc="upper left")
    plt.show()
    
    
    for h in h_tab:
        for x in Xn:
            z=np.random.normal(0, 1)
            yn=math.atan(2*x)+z
            Yn.append(yn)
        plt.plot(x_d,est_nadaray_watson(x_d,Xn,Yn,epanechnikov_kernel,h),label=' $h = %s $' % h)
        Yn.clear()
    plt.title("Estymacja charakterystyki systemu z różnym parametrem h" )
    plt.legend(loc="upper left")
    plt.show()
    
    
    plt.scatter(Xn,Yn,label='chmura pomiarów')
    plt.plot(x_d,atan_tab,color='red',label='nieliniowa charakterysytyka systemu')
    plt.title("Nieliniowa charakterystyka systemu z chmurą pomiarów")
    plt.legend(loc="upper left")
    plt.show()
    """

    