# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 14:42:49 2020

@author: Hubert
"""

import random
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import signal
from scipy import special
import scipy.stats as st
from statsmodels.distributions.empirical_distribution import ECDF


#generator o rozkladzie jednostajnym na przedziale [0,1]
def gen_sin(sintab,n,w):
    sintab.append(random.uniform(0,1))
    for x in range(1,n):
        sintab.append((1+np.sin(w*sintab[x-1]))/2)
    plt.hist(sintab, bins=500)
    plt.title('Histogram generatora sinusoidalnego')
    plt.ylabel('liczba wystąpień próbki')
    plt.xlabel("zmienna losowa")
    plt.show()
    
def gen_sawtooth(sawtab,n,w):
    sawtab.append(random.uniform(0,1))
    for x in range(1,n):
        sawtab.append((1+signal.sawtooth(w*sawtab[x-1]))/2)
    #plt.hist(sawtab, bins=500)
    #plt.title('Histogram generatora piłokształtnego')
    #plt.ylabel('liczba wystąpień próbki')
    #plt.xlabel("zmienna losowa")
    #plt.show()

def gen_fibo(fibotab,n,fibo_p,fibo_q, fibo_mod):
    for x in range(0, fibo_p+1):
        fibotab.append(random.uniform(0,1))
    for x in range(fibo_p+1,n):
        fibotab.append((fibotab[x-fibo_p] - fibotab[x-fibo_q])%fibo_mod)
    plt.hist(fibotab, bins=500)
    plt.title('Histogram generatora Fibonaciego')
    plt.ylabel('liczba wystąpień próbki')
    plt.xlabel("zmienna losowa")
    plt.show()

def mediana(tab):
  tab.sort()
  if (len(tab) % 2 == 0):
    return (tab[int(len(tab)/2-1)] + tab[int(len(tab)/2)])/2
  else:
    return tab[int(len(tab)/2)]

def roz_wykladniczy(x,lambd):
    return -1/lambd*np.log(1-x) 

def roz_logistyczny(x,mi,s):
    return np.log(x/(1-x))*s+mi

def roz_cauchy(x,x0,gamma):
    return gamma*np.tan(np.pi*(x-0.5))+x0

def roz_normalny(x,mi,delta):
    return delta*np.sqrt(2)*special.erfinv(2*x-1)

def roz_trojkatny(x):
    tmp=random.uniform(0,1)
    if (tmp >= 0.5):
       return 1-np.sqrt(1-x)
    else:
        return -1+np.sqrt(1-x)

def roz_laplace(triantab,fixedtab,mi,b):   
    for x in triantab:
        if (x>0.5):
          y=mi-b*np.log(2-2*x)
        else:
          y=mi+b*np.log(2*x)
        fixedtab.append(y)

def roz_trojkatny2(x):
    if (x < 0.5):
       return np.sqrt(2*x)-1
    else:
        return 1-np.sqrt(2-2*x)    