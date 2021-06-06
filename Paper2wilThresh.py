# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 13:39:59 2020

@author: jwilliams
"""

def jwilliThresh(signal):
    models = []
    E,O = Split1D(signal)
    coeffsE = pywt.wavedec(E,'db8'); coeffsO =  pywt.wavedec(O,'db8'); 
    Erec = pywt.waverec(coeffsE,'db8'); Orec = pywt.waverec(coeffsO,'db8'); 
    ERec = [np.mean(Erec[i:i+2]) for i in range(np.size(Erec)-1)]; ORec = [np.mean(Orec[i:i+2]) for i in range(np.size(Orec)-1)]
    SQ=(np.sum((ERec-O[:np.size(O)-1])**2) +  np.sum((ORec-E[1:])**2))
    s = np.shape(coeffsE)[0]
    models = []; threshT=[]; compons = [];
    for j in range(1,s):
        w = int(np.log(np.size(coeffsE[j])));
        comps = 2
        mod = []; mod = hmm.GaussianHMM(n_components=comps,n_iter = 10000);
        BlocksE = []; BlocksO = []
        for i in range(0,np.size(coeffsE[j])-w,w):
            BlocksE.append(coeffsE[j][i:i+w])
            BlocksO.append(coeffsO[j][i:i+w])
        mod.fit([BlocksE[m] for m in range(np.shape(BlocksE)[0])])
        mod.fit([BlocksO[m] for m in range(np.shape(BlocksO)[0])])
        models.append(mod)
        tempE = np.asarray(np.hstack(mod.predict([(BlocksE[m]) for m in range(np.shape(BlocksE)[0])])));
        tempO = np.asarray(np.hstack(mod.predict([(BlocksO[m]) for m in range(np.shape(BlocksO)[0])])));  
        thresh = []
        for i in range(comps):
            coeffsEt = copy.deepcopy(coeffsE); coeffsOt = copy.deepcopy(coeffsO); 
            for b in range(np.size(tempE)):
                if tempE[b]==i:
                    coeffsEt[j][(b*w):(b*w)+w]=0
                if tempO[b]==i:
                    coeffsOt[j][(b*w):(b*w)+w]=0
            Erec = pywt.waverec(coeffsEt,'db8'); Orec = pywt.waverec(coeffsOt,'db8')
            ERec = [np.mean(Erec[i:i+2]) for i in range(np.size(Erec)-1)]; ORec = [np.mean(Orec[i:i+2]) for i in range(np.size(Orec)-1)]
            SQt = (np.sum((ERec-O[:np.size(O)-1])**2) +  np.sum((ORec-E[1:])**2))
            if SQt<=SQ:
                thresh.append(i)
        threshT.append(thresh)
    c = HMMRecon(signal,threshT, models)
    return pywt.waverec(c,'db8')

def HMMRecon(signal, threshT, models):
    coeffs = pywt.wavedec(signal,'db8')
    s = np.shape(coeffs)[0]; Blocks = []
    for j in range(1,s):
        w = int(np.log(np.size(coeffs[j])))
        Blocks = []
        for i in range(0,np.size(coeffs[j])-w,w):
            Blocks.append(coeffs[j][i:i+w])
        if j == s-1:
            coeffs[j][:]=0
        else:
            temp = np.hstack(models[j-1].predict([Blocks[m] for m in range(np.shape(Blocks)[0])]));
            tcomp = threshT[j-1]
        for i in range(np.size(tcomp)):
             for b in range(np.size(temp)):
                 if temp[b]==tcomp[i]:
                     coeffs[j][(b*w):(b*w)+w]=0
    return coeffs