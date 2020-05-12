# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 14:34:01 2020

@author: Hubert
"""
from generatory import *
import random
import numpy as np
import array
import matplotlib.pyplot as plt
import math
from scipy import signal
from scipy import special
import scipy.stats as st
from statsmodels.distributions.empirical_distribution import ECDF

def roz_normalny(x,mi,delta):
    return delta*np.sqrt(2)*special.erfinv(2*x-1)+mi


#estymator wartosci oczekiwaniej (nieobciazony)
def est_wo(tab):
    suma=0
    for x in tab:
        suma+=x
    return (suma/len(tab))

#estymator wariancji (obciazony)
def est_varob(tab):
    suma=0
    war_oczekiwana=est_wo(tab)
    for x in tab:
        suma+=(x**2-war_oczekiwana**2)
    return (suma/len(tab))

#estymator wariancji (nieobciazony)
def est_varnob(tab):
    suma=0
    war_oczekiwana=est_wo(tab)
    for x in tab:
        suma+=(x-war_oczekiwana)**2
    return (suma/(len(tab)-1))

#estymator kowariancji
def est_cov(tab1,tab2):
    cov_tab=[]
    for i in range(0,len(tab1)):
        cov_tab.append(tab1[i]*tab2[i])
    ex=est_wo(tab1)
    ey=est_wo(tab2)
    exy=est_wo(cov_tab)
    return (exy-ex*ey)
        
#estymator korelacji (wspolczynnik korelacji liniowej Pearsona)
def est_corr(tab1,tab2):
    cov=est_cov(tab1,tab2)
    return cov/(np.sqrt(est_varnob(tab1)*est_varob(tab2)))

#obciazenie estymatora (bias)
def bias(tab, mi):
    return (est_wo(tab)-mi)

#blad sredniokwadratowy estymatora (empiryczny)
def multitab_generator(N,L,tabtab):
    resulttab=[]
    testtab= []
    for i in range(0,L):
        gen_sawtooth(testtab, N, 88000)
        for x in testtab:
            y=roz_normalny(x, 5, 1)
            resulttab.append(y)
        tabtab.append(resulttab.copy())
        testtab.clear()
        resulttab.clear()
        
def MSE_err(N,L):
    mse_errtab=[]
    uni_tab=[] #tablica z roz. jednostajnym do tworzenia kolejnych rozkladow
    testtab=[] #tablica z rozkładem do estymacji 
    for n in N:
        mse_sum=0
        for i in range(0,L):
            gen_sawtooth(uni_tab,n,88000)
            for x in uni_tab:
                y=roz_normalny(x,5,1)
                testtab.append(y)
            mse_sum+=(est_wo(testtab)-5)**2
            uni_tab.clear()
            testtab.clear()
        err=mse_sum/L
        mse_errtab.append(err)
    return mse_errtab
        
    

if __name__=='__main__':
    """testtab=[]
    resulttab=[]
    resulttab2=[]
    gen_sawtooth(testtab, 10000, 88000)
    for x in testtab:
        y=roz_normalny(x, 5, 1)
        z=roz_wykladniczy(x, 0.5)
        resulttab.append(y)
        resulttab2.append(z)
    
    plt.hist(resulttab, bins=500)
    plt.title('Histogram rozkładu')
    plt.ylabel('liczba wystąpień próbki')
    plt.xlabel("zmienna losowa")
    plt.show()
    
    print("wariancja z estymatora obciążonego: %.10f" % est_varob(resulttab))  
    print("wariancja z estymatora nieobciążonego: %.10f" % est_varnob(resulttab))      
    print("wartosc oczekiwana z estymatora: %.10f" % est_wo(resulttab))
    print("dolne ogranicznenie wariancji: %.10f" % (est_varnob(resulttab)/len(resulttab)))
    print("kowariancja z dwoch rozkladow: %.10f" % est_cov(resulttab,resulttab2))
    print("współczynnik korelacji pearsona z dwoch rozkladow: %.10f" % est_corr(resulttab,resulttab2))
    testtab.clear()
    resulttab.clear()
    resulttab2.clear()
    multitab= []
    
    mi=5
    
    multitab_generator(100,10,multitab)
    mse_suma=0
    for i in range(0,10):
        wo=est_wo(multitab[i])
        y=(wo-mi)**2
        mse_suma+=y
       """
    #numbers_list = [10, 50, 100, 250, 500, 1000, 2000]
    N=[]
    for i in range(10,1000,25):
        N.append(i)
    L=10
    resulttab=[]
    resulttab=MSE_err(N,L)
    plt.plot(N,resulttab)
    plt.title("Błąd empiryczny MSE przy L=10")
    plt.xlabel("liczba próbek n")
    plt.ylabel("wartość błędu")
    plt.show()
    #print("MSE empiryczny dla normalnego L=10 N=100 mi=5: %.10f" %mse_suma)


