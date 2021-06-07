# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 11:50:40 2020

@author: jwilliams
"""
def vShrink(signal):
    coeffs = pywt.wavedec(signal,'db6'); s = np.shape(coeffs)[0]; 
    stand = np.std(coeffs[s-1]); n = np.size(signal); t = stand*np.sqrt(2*np.log(n)) 
    for j in range(1,s):
        coeffs[j][np.abs(coeffs[j])<=t] = 0
        for i in range(np.size(coeffs[j])):
            if np.abs(coeffs[j][i]) > t:
                     if coeffs[j][i] > 0:
                         coeffs[j][i] = coeffs[j][i]-t
                     else:
                         coeffs[j][i] = coeffs[j][i]-t
    rv = pywt.waverec(coeffs,'db6');
    return rv

def visushrinkS(Signal):
    coeffs = pywt.wavedec(Signal,'haar')
    for res in range(np.shape(coeffs)[0]):
        sigma = np.std(coeffs[np.shape(coeffs)[0]-1])
        for i in range(1,np.size(coeffs[res])):
            if coeffs[res][i]< sigma*np.sqrt(2*np.log(np.size(Signal))):
                coeffs[res][i]=0
    recon = pywt.waverec(coeffs,'haar')
    return recon

def visushrinkC(coeffs,Signal):
    for res in range(np.shape(coeffs)[0]):
        sigma = np.std(coeffs[np.shape(coeffs)[0]-1])
        for i in range(1,np.size(coeffs[res])):
            if coeffs[res][i]< sigma*np.sqrt(2*np.log(np.size(Signal))):
                coeffs[res][i]=0
    return coeffs


def LevVshrink(coeffs):
        n = np.size(coeffs)
        sigma = np.std(coeffs)
        thresh = sigma*np.sqrt(2*np.log(n))
        coeffsR = copy.deepcopy(coeffs)
        coeffsR[np.abs(coeffsR)<abs(thresh)]=0
        return coeffsR, thresh