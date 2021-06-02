# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 10:20:48 2020

@author: jwilliams
"""

import math
import scipy.stats as sci
import pandas as pd 
import matplotlib.pyplot as plt
import pywt
import numpy as np 
from numpy import linalg as LA
import statsmodels.api as sm
from scipy import stats
import bisect 
import WANOVACent as wc
import csv 
import random
import numpy as np

def wtest(s, p, sig, e, runs):
    test, signal = wtc(s, p, sig, e, runs)
    MSEV = []; MSEJ = []; MSEM = []
    for i in range(np.shape(test)[0]):
        MSEV.append(np.sum((vShrink(test[i])-signal)**2))    
        #MSEJ.append(np.sum((HMMWTest2(pywt.waverec(BT1D(test[i])[0],'db6'))-signal)**2))
        MSEM.append(np.sum((jwilliThresh(test[i])-signal)**2))
    #MSEJ = np.asarray(MSEJ); MSEV = np.asarray(MSEV)
    MSEM = np.asarray(MSEM);
    #print(np.mean(MSEJ/MSEV)); print(np.std(MSEJ/MSEV));
    print(np.mean(MSEM/MSEV)); print(np.std(MSEM/MSEV));
    return test[i], signal
   
def wtc(s, p, sig, e, runs):            
    powerCase = {
        1: 3,
        2: 5
        }
    
    signalCase = {
        1: Blocks(s),
        2: Bumps(s),
        3: Doppler(s),
        4: HeaviSine(s),
        5: Wave(s)
            }
    
    
    Test = []
    signal = signalCase[sig];
    sp = (np.sum(np.abs(signal))*(np.sum(np.abs(signal))/np.size(signal))) 
    power = powerCase[p];
    for r in range(runs):
        #error needs to be defined each run
        errorCase = {
        1: np.random.standard_t(3,s),
        2: np.random.standard_cauchy(s),
        3: np.random.lognormal(0,1,s),
        4: np.random.normal(0,1,s)
        }
        error = errorCase[e]
        npow = np.sum(np.abs(error))*(np.sum(np.abs(error))/np.size(error))
        k = (sp/npow)*(10**(-power/10))
        Test.append(signal+k*error)
    return Test, signal







                
                    
def Blockstest(size):
    #l = 1024
    #s = np.zeros((l,100))
    b = Blocks(size)
    Spower = np.sum(np.abs(b))*(np.sum(np.abs(b))/np.size(b))
    #n = (np.random.standard_cauchy(2048))
    #n = np.random.lognormal(0,1,size)
    n = np.random.normal(0,1,size)
    Npower = np.sum(np.abs(n))*(np.sum(np.abs(n))/np.size(n))
    k = (Spower/Npower)*(10**(-5/10))
    return b+n*k
    '''
   
    for i in range(1): 
        B = Blocks(size)
        e = np.random.standard_t(1,size)
        #e = np.random.standard_cauchy(2048)
        #e = np.random.lognormal(size=l)
        #e = np.random.normal(0,1,2048)
        #e[(e>50)]=0
        #e[e<-50]=0
        s = B+e*k
    #s[:,50] = BlocksO()+np.random.normal(0,1,2048)*.5
    #e = wanovaBoxPlot(s.T)
    #e = 0
    return s
    '''

def Bumpstest(size):
    b = Bumps(size)
    Spower = np.sum(np.abs(b))*(np.sum(np.abs(b))/np.size(b))
    #n = (np.random.standard_cauchy(2048))
    #n = np.random.lognormal(0,1,size)
    n = np.random.normal(0,1,1024)
    Npower = np.sum(np.abs(n))*(np.sum(np.abs(n))/np.size(n))
    k = (Spower/Npower)*(10**(-5/10))
    return b+k*n
'''    
for i in range(100): 
        B = Bumps()
        e = np.random.normal(0,1,2048)
        #e = np.random.standard_cauchy(2048)
        #e[(e>50)]=0
        #e[e<-50]=0
        #e = np.random.lognormal(0,1,2048)
        s[:,i] = B+e*k
    e = 0
    #e = wanovaBoxPlot(s.T)
    return s,e
'''
def HeaviSinetest():
    s = np.zeros((2048,100))
    b = HeaviSine()
    Spower = sum(abs(b))*(sum(abs(b))/np.size(b))
    #n = (np.random.standard_cauchy(2048))
    #n = np.random.lognormal(0,1,2048)
    n = np.random.normal(0,1,2048)
    Npower = sum(abs(n))*(sum(abs(n))/np.size(n))
    k = (Spower/Npower)*(10**(-5/10))
    s = np.zeros((2048,100))
    for i in range(100): 
        B = HeaviSine()
        e = np.random.normal(0,1,2048)
        #e = np.random.standard_cauchy(2048)
        #e[(e>50)]=0
        #e[e<-50]=0
        #e = np.random.lognormal(0,1,2048)
        s[:,i] = B+e*k
    #e = compWanova2(s.T)
    e = 0
    return s,e

def Dopplertest():
    s = np.zeros((2048,100))
    b = Doppler()
    Spower = sum(abs(b))*(sum(abs(b))/np.size(b))
    #n = (np.random.standard_cauchy(2048))
    #n = np.random.lognormal(0,1,2048)
    n = np.random.normal(0,1,2048)
    Npower = sum(abs(n))*(sum(abs(n))/np.size(n))
    k = (Spower/Npower)*(10**(-5/10))
    s = np.zeros((2048,100))
    for i in range(100): 
        B = Doppler()
        e = np.random.normal(0,1,2048)
        #e = np.random.standard_cauchy(2048)
        #e[(e>50)]=0
        #e[e<-50]=0
        e = np.random.lognormal(0,1,2048)
        s[:,i] = B+e*k
    #e = compWanova2(s.T)
    e=0
    return s,e

def Step():
     t= np.linspace(0, 1, 2048)
     I = np.ones(2048)
     for i in range(np.size(t)):
         if t[i] <= (1/3):
             I[i]=0
         if t[i] >= 0.75:
             I[i] = 0
     f = .2+.6*I*t
     return f;
    
def Wave(size):
    t= np.linspace(0, 1, size)
    f = .5+2*np.cos(4*np.pi*t)+.1*np.cos(24*np.pi*t)
    return f;

def Blocks(size):
    t_j = [0.1, 0.13, 0.15, 0.23, 0.25, 0.4,0.44,0.65,0.76,0.78, 0.81]
    h_j = [4,-5,3,-4,5,-4.2,2.1,4.3,-3.1,2.1,-4.2]
    t= np.linspace(0, 1, size)
    temp = np.zeros(np.size(t_j))
    f = np.zeros((np.size(t_j),np.size(t)))
    for j in range(np.size(t_j)):
        f[j,:]  = h_j[j]*((1+np.sign(t-t_j[j]))/2)
    return sum(f)             
    
def BlocksO():
    t_j = [0.1, 0.13, 0.15, 0.23, 0.25, 0.4,0.44,0.65,0.76,0.78, 0.81]
    h_j = [10,-5,3,-4,5,-4.2,2.1,4.3,-3.1,2.1,-4.2]
    t= np.linspace(0, 1, 1024)
    temp = np.zeros(np.size(t_j))
    f = np.zeros((np.size(t_j),np.size(t)))
    for j in range(np.size(t_j)):
        f[j,:]  = h_j[j]*((1+np.sign(t-t_j[j]))/2)
    return sum(f)  

def Bumps(size):
    t= np.linspace(0, 1, size)
    t_j = [0.1, 0.13, 0.15, 0.23, 0.25, 0.4,0.44,0.65,0.76,0.78, 0.81]
    h_j = [4,5,3,4,5, 4.2, 2.1,4.3,3.1,5.1,4.2]
    w_j = [.005,.005,.006,.01,.01,.03,.01,.01,.005,.008,.005]
    f = np.zeros((np.size(t_j),np.size(t)))
    for j in range(np.size(t_j)):
        f[j,:]  = h_j[j]*((1+abs((t-t_j[j])/w_j[j]))**-4)
    return sum(f) 

def BumpsO():
    t= np.linspace(0, 1, 2048)
    t_j = [0.125, 0.15, 0.175, 0.23, 0.25, 0.4,0.44,0.65,0.76,0.78, 0.81]
    h_j = [4000000000000,5,3,4,5, 4.2, 2.1,4.3,3.1,5.1,4.2]
    w_j = [.005,.005,.006,.01,.01,.03,.01,.01,.005,.008,.005]
    f = np.zeros((np.size(t_j),np.size(t)))
    for j in range(np.size(t_j)):
        f[j,:]  = h_j[j]*((1+abs((t-t_j[j])/w_j[j]))**-4)
    return sum(f) 

def HeaviSine(size):
     t= np.linspace(0, 1, size)
     f = np.zeros(np.size(t))
     f = 9*np.sin(4*np.pi*t)-np.sign(t-.3)-np.sign(.72-t)
     return f
     
    
def Doppler(size):
    t = np.linspace(0,1, size)
    f = (t*(1-t))**(.5)*np.sin((2*np.pi)*((1+.05)/(t+.05)))
    return f