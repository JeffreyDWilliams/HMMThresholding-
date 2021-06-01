# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 17:55:16 2020

@author: jwilliams
"""

from hmmlearn import hmm 
                      
def TreeChopper(signal):
   recon = copy.deepcopy(signal); sel = np.zeros(np.size(signal));
   E,O = Split1D(D1); w = pywt.Wavelet('haar'); SST1 = []; SST2 = [];
   mL = pywt.dwt_max_level(np.size(E), filter_len=w.dec_len);
   cE = pywt.wavedec(E,'haar'); cO = pywt.wavedec(O,'haar')
   tE = copy.deepcopy(cE); tO = copy.deepcopy(cO)
   for res in range(0,mL+1):
       tE[res][:] = 0; tO[res][:] = 0
   BE = copy.deepcopy(tE); BO = copy.deepcopy(tO);
   for Tree in range(0,np.size(cE[mL])):
       WE = copy.deepcopy(tE); WO = copy.deepcopy(tO);
       for lev in range(1,mL):
           coeff = np.ceil(Tree/(2**(mL-(lev))))-1
           WE[lev][coeff] = cE[lev][coeff]
           WO[lev][coeff] = cO[lev][coeff]
       rO = pywt.waverec(WO,'haar'); rE = pywt.waverec(WE,'haar')
       e_hat = np.asarray([np.mean(rE[i:i+2]) for i in range(np.size(rE)-1)]);
       o_hat = np.asarray([np.mean(rO[i:i+2]) for i in range(np.size(rO)-1)]);
       SS = .5*(np.sum((O[:np.size(O)-1]-e_hat)**2) + np.sum((E[1:]-o_hat)**2))
       SST1.append([SS,SS]);
   
   SST1 = np.asarray(np.hstack(SST1));
   mS = np.mean(SST1)
   '''
                    SST1[0:10]=SST1[11]
   SST1[SST1>np.percentile(SST1,90)]=0; SST1[SST1>0]=.8;
   SST1 = SST1+.2
   r  = HMMBThresh(signal, SST1)
   '''
   mod=[]
   resAdd = np.zeros(np.size(SST1))
   mod = hmm.GaussianHMM(n_components=20,n_iter = 1000000)
   mod.startprob_ = [.1,.1,.1,.1,.1,.1,.1,.1,.1,.1]
   mod.fit([[SST1[n]] for n in range(np.size(SST1))])
   fm1 = np.asarray(mod.predict([[SST1[n]] for n in range(np.size(SST1))]))
   fm1[0:10] = np.mean(fm1[0:20])
   for i in range(20):
       SST1[fm1==i] = min(1,np.mean(SST1[fm1==i])/mS)
       if np.mean(SST1[fm1==i])>np.percentile(SST1,80):
           resAdd[fm1==i]=1
   r = HMMBThresh(signal, SST1, resAdd)
   return r
   
def HMMBThresh(signal, *r, **resAdd):
    r =1
    resAdd =0
    r = np.ones(np.size(signal)/2)
    resAdd = np.zeros(np.size(signal)/2);
    E,O = Split1D(signal)
    coeffsE = pywt.wavedec(E,'db8'); coeffsFE = copy.deepcopy(coeffsE);
    coeffsO =  pywt.wavedec(O,'db8'); coeffsFO = copy.deepcopy(coeffsO);
    Erec = pywt.waverec(coeffsE,'db8');Orec = pywt.waverec(coeffsO,'db8'); 
    ERec = [np.mean(Erec[i:i+2]) for i in range(np.size(Erec)-1)]; ORec = [np.mean(Orec[i:i+2]) for i in range(np.size(Orec)-1)]
    SQ=.5*(np.sum(r[:np.size(O)-1]*(ERec-O[:np.size(O)-1])**2) +  np.sum(r[1:]*(ORec-E[1:])**2))
    s = np.shape(coeffsE)[0]
    models = []; threshT=[]; compons = [];
    for j in range(1,s):
        comps = np.max([4,(s-j)+2]);
        #comps = 11
        compons.append(comps)
        mod = []; mod = hmm.GaussianHMM(n_components=comps,n_iter = 1000000);
        probs = [];
        temp = []
        #Generate probability of starting in each node 
        '''
        for i in range(comps):
            if i == 0:
                probs.append(.8+.2*((j)/s))
                val = (1-probs[0])/(comps)
            probs.append(val)
        '''
        #Assign start probabilities
        #mod.startprob_ = np.hstack(probs)
        #Determine the nodes where each wavelet coefficients fall 
        for i in range(10):
            mod.startprob_ = [np.random.triangular(0.0,0.1,.4) for n in range(comps)]
            mod.fit([[coeffsE[j][n]]for n in range(np.size(coeffsE[j]))]);
            mod.fit([[coeffsO[j][n]]for n in range(np.size(coeffsO[j]))]);
        models.append(mod)
        tempE = np.asarray(np.hstack(mod.predict([[coeffsE[j][n]]for n in range(np.size(coeffsE[j]))])));
        tempO = np.asarray(np.hstack(mod.predict([[coeffsO[j][n]]for n in range(np.size(coeffsO[j]))])));  
        thresh = []
        for i in range(comps):
            coeffsEt = (copy.deepcopy(coeffsE)); coeffsOt = (copy.deepcopy(coeffsO))
            coeffsEt[j][tempE==i]=0; coeffsOt[j][tempO==i]=0;
            Erec = pywt.waverec(coeffsEt,'db8'); Orec = pywt.waverec(coeffsOt,'db8')
            ERec = [np.mean(Erec[i:i+2]) for i in range(np.size(Erec)-1)]; ORec = [np.mean(Orec[i:i+2]) for i in range(np.size(Orec)-1)]
            SQt = .5*(np.sum(r[:np.size(O)-1]*(ERec-O[:np.size(O)-1])**2) +  np.sum(r[1:]*(ORec-E[1:])**2))
            if SQt<=SQ:
                thresh.append(i)
        threshT.append(thresh)
    c = HMMRecon(signal,threshT, models, resAdd)
    return pywt.waverec(c,'db8')

def HMMRecon(signal, threshT, models,resAdd):
    coeffs = pywt.wavedec(signal,'db8')
    s = np.shape(coeffs)[0]
    for j in range(1,s):
         if j == s-1:
             temp = np.asarray(np.hstack(models[j-2].predict([[coeffs[j][n]]for n in range(np.size(coeffs[j]))])));
             tcomp = threshT[j-2]
             #coeffs[j][:]=0
         else:
             temp = np.asarray(np.hstack(models[j-1].predict([[coeffs[j][n]]for n in range(np.size(coeffs[j]))])));
             tcomp = threshT[j-1]
         for i in range(np.size(tcomp)):
             coeffs[j][temp==tcomp[i]]=0
         coeffs[s-1][resAdd==0] = 0
    return coeffs
        