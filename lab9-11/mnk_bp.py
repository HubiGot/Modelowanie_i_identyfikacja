# -*- coding: utf-8 -*-
"""
Created on Fri May 29 13:23:19 2020

@author: Hubert
"""

import numpy as np
import matplotlib.pyplot as plt
#from labgotowe import random


D = 10
N = 200
L = 10
var = 1
tab_n = np.linspace(16, 100, 42).astype(int)
b = 0.5


if __name__ == "__main__":
    # generate a
    a = np.zeros(D)
    for i in range(D):
        a[i] = 1
    # generate Xn
    Xn = np.zeros((N, 10), dtype=float)
    for i in range(100):
        for j in range(D):
            Xn[i][j] = np.random.normal(0, 1)
    # generate Z
    epsilon = np.zeros(N + 1)
    for i in range(len(epsilon)):
        epsilon[i] = np.random.normal(0, var)
    Zn = np.zeros(N)
    for i in range(N):
        Zn[i] = epsilon[i + 1] - b * epsilon[i]
    Yn = np.add(Xn.dot(a), Zn)
    a_est = np.dot(np.dot(np.linalg.inv(np.dot(Xn.transpose(), Xn)), Xn.transpose()), Yn)
    print(a_est)
    #covZ
    R = np.zeros((N,N),dtype=float)
    for i in range(N):
        for j in range(N):
            if i==j:
                R[i][j] = (1+pow(b,2))*var
            elif np.absolute(i-j)==1:
                R[i][j] = b*var
    #cova
    cova = np.dot(np.dot(np.dot(np.dot(np.linalg.inv(np.dot(Xn.transpose(),Xn)),Xn.transpose()), R),Xn),np.linalg.inv(np.dot(Xn.transpose(),Xn)))
    error_tab = []
    Zn = np.zeros((L,N),dtype=float)
    tab_epsilon = np.zeros((L,N+1),dtype=float)
    for l in range(L):
        for i in range(N):
            tab_epsilon[l][i] = np.random.normal(0, var)
    for l in range(L):
        for i in range(N):
            Zn[l][i] = tab_epsilon[l, i + 1] - b * tab_epsilon[l, i]
    for n in tab_n:
        x_tmp = Xn[:n,:]
        error = 0
        for l in range(L):
            z_tmp = Zn[l][:n]
            y_tmp = np.add(x_tmp.dot(a), z_tmp)
            a_estimated_tmp = np.dot(np.dot(np.linalg.inv(np.dot(x_tmp.transpose(),x_tmp)),x_tmp.transpose()),y_tmp)
            error += 1/L * pow((np.linalg.norm(a_estimated_tmp-a)),2)
        error_tab.append(error)
    plt.figure(2)
    plt.plot(tab_n, error_tab)
    plt.show()