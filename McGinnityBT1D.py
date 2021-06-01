# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 10:51:05 2020
@author: jwilliams
Brings in a signal and thresholds it 
"""
import pywt
import numpy as np 


def BT1D(Signal):
    E, O = Split1D(Signal)
    w = pywt.Wavelet('db6')
    mL = pywt.dwt_max_level(np.size(E), filter_len=w.dec_len)
    coeffsF = pywt.wavedec(Signal,'db6',level=mL)
    coeffsE = pywt.wavedec(E,'db6')
    coeffsO = pywt.wavedec(O,'db6')
    mini = []
    for j in range(1,int(np.shape(coeffsE)[0])):
        Obj = []
        Thresh = np.hstack([ThreshRange1D(coeffsE[j]), ThreshRange1D(coeffsO[j])])
        for t in enumerate(Thresh):
            score = sum((TRecon(coeffsE,t[1],j)-O)**2+(TRecon(coeffsO,t[1],j)-E)**2)+ sum((TRecon(coeffsE,t[1],j)-TRecon(coeffsO,t[1],j))**2)
            Obj.append(score)
        mini.append(np.argmin(Obj))
    coeffs, df = FRecon(coeffsF,  mini)
    #return pywt.waverec(coeffs,'db6'), df;
    return(coeffs)
def WFRecon(Signal, T):
    c = pywt.wavedec(Signal,'db6',level=2)
    for j in range(0,3):
        w =4
        #w = int(np.log(np.size(c[j])))
        for i in range(0,int(np.size(c[j])),w):
            SS = np.sum((c[j][i:i+w]**2))
            if SS <= T[j-1]:
                c[j][i:i+w] = 0
    return c 

def FRecon(c, t):
    df = 0
    for j in range(1,int(np.shape(c)[0])): 
        w=4
        #w = int(np.log(np.size(c[j])))
        thresh = t[j-1]*(1-((np.log(2))/np.log(np.size(c[j])/2**(j))))**-1
        for i in range(0,int(np.size(c[j])),w):
            SS = sum((c[j][i:i+w]**2))
            if SS <= thresh:
                c[j][i:i+w] = 0
                df = df+w
    return c, df;
            
def TRecon(coeffs, t,j):
     w = int(np.log(np.size(coeffs[j])))
     w = 4
     for i in range(1,int(np.size(coeffs[j])),w):
         SS = sum(coeffs[j][i:i+w]**2)
         if SS <= t:
             coeffs[j][i:i+w] = 0
     return pywt.waverec(coeffs, 'db6');
 
def ThreshRange1D(c): 
    w = int(np.log(np.size(c)))
    w=4
    maxi = []
    for i in range(0,int(np.size(c)),w):
            maxi.append(np.sum((c[i:i+w]**2)))
    return maxi;
    
def Split1D(Signal):
    E = Signal[::2]
    O = Signal[1::2]
    return E, O;



