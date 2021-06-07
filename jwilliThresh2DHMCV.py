# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 10:40:48 2020

@author: jwilliams
"""
ngimg = graphr
coeffs = pywt.wavedec2(ngimg, 'db6', level = 4)
for j in range(1,5):
    for dim in range(np.shape(coeffs[j])[0]):      
        comps= 9
        mod = []; mod = hmm.GaussianHMM(n_components=comps,n_iter = 100);
        ch = filters.roberts(coeffs[j][dim])
        w = int(np.log(np.sqrt(np.size(coeffs[j][dim])))); 
        chs = []
        for x in range(0,np.shape(coeffs[j][dim])[0]-w,w):
            for y in range(0,np.shape(coeffs[j][dim])[1]-w,w):
                chs.append(np.hstack(ch[x:x+w,y:y+w]))
        mod.fit([chs[m] for m in range(np.shape(chs)[0])])
        states = mod.predict([chs[m] for m in range(np.shape(chs)[0])])
        s=0
        blocks = copy.deepcopy(ch)
        for x in range(0,np.shape(coeffs[j][dim])[0]-w,w):
           for y in range(0,np.shape(coeffs[j][dim])[1]-w,w):
               blocks[x:x+w,y:y+w][:]=states[s] 
               s=s+1
        if x < np.shape(coeffs[j][dim])[0]-1:
            val = np.shape(blocks)[0]-x-1
            blocks[x:x+val,:] = blocks[x-val:x,:]
        if y < np.shape(coeffs[j][dim])[1]-1:
            val = np.shape(blocks)[1]-y-1
            blocks[:,y:y+val] = blocks[:,y-val:y]
        states = np.hstack(blocks);
        ch = np.hstack(ch); c = np.hstack(coeffs[j][dim]);
        var = np.var([coeffs[2][0],coeffs[2][1],coeffs[2][2]]);
        for s in range(comps):
            sthresh =1-((var*np.size(ch[states==s])*np.sqrt(4.5))/np.linalg.norm(ch[states==s])**2)
            if sthresh<0:
                sthresh = 0
            for i in range(np.size(c)):
                if states[i] ==s:
                    c[i] = c[i]*sthresh
        coeffs[j][dim][:] = np.reshape(c,[np.shape(coeffs[j][dim])[0],np.shape(coeffs[j][dim])[1]])
        if np.count_nonzero(coeffs[j][dim])/np.size(coeffs[j][dim])== 1:
            coeffs[j][dim][:]=0
r = pywt.waverec2(coeffs,'db6')
print(np.mean((gimg-r)**2))