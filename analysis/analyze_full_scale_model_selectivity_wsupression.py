'''
Extract S_NMDA, S_GABA values centered by time in which the trial reaches the threshold.

* Notes:
- After the decision is made, the S_* values are taken from simulation trials with inhibitory suppression.
'''

import os
import warnings
import numpy as np
import matplotlib.pyplot as plt
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

'''=====Define parameters here====='''
T = np.arange(-60, 15, 5)
T2 = np.arange(15, 60, 5)
M = len(T); M2 = len(T2)
tnum = 24 # number of trials to include ## 228
gmavals = [(50,1,5,0.0052,2.5)] 
mother = os.getcwd()
sup, thre, delay = 0.005, 0.17, 0.15 # suppresion, threshold, delay
'''================================'''

def unpack(v):
    return v[:Ne], v[Ne:2*Ne], v[2*Ne:2*Ne+Ni], v[2*Ne+Ni:N], v[N:]

# main
for k in range(len(gmavals)):

    Ne, p, nsig, mu0, nsigg = gmavals[k]
    Ni = Ne; N = 2*(Ne + Ni)
    fname = 'Ne='+str(Ne)+', p='+str(p)+', nsig='+str(nsig)+', mu0='+str(mu0)+', nsigg='+str(nsigg)
    
    N_win, N_lose = {}, {} # NMDA
    G_win, G_lose = {}, {} # GABA
    for n in range(M+M2):
        N_win[n] = []; N_lose[n] = []
        G_win[n] = []; G_lose[n] = []
    
    rts = []
    for rp in range(tnum):

        # extract simulation parameters for this trial
        resy = np.load(os.path.join(fname,'run'+str(1)+'_tr'+str(rp)+'.npy'))
        s1s, s2s, sg1s, sg2s, noises = unpack(resy)
        mean_s1, mean_s2 = np.mean(s1s,axis=0), np.mean(s2s,axis=0)
        rt, win = hpr.accuracy(mean_s1, mean_s2, thre=thre)
        
        if rt != None: # if reaction time exists

            # before decision is made
            rts.append(rt)
            cond1 = (win == 1) # if winner is population 1
            for ne in range(Ne):
                s1, s2, sg1, sg2 = s1s[ne], s2s[ne], sg1s[ne], sg2s[ne]
                
                for n in range(M):
                    t = T[n]; idx = int(rt*100 + t + 50)
                    cond2 = (idx >= 0) and (idx < 300) # if time point considered is within existing time points
                    
                    if cond1 and cond2:
                        N_win[n].append(s1[idx])
                        G_win[n].append(sg1[idx])
                        N_lose[n].append(s2[idx])
                        G_lose[n].append(sg2[idx])
                    elif not cond1 and cond2:
                        N_win[n].append(s2[idx])
                        G_win[n].append(sg2[idx])
                        N_lose[n].append(s1[idx])
                        G_lose[n].append(sg1[idx])
                        
            # after decision is made
            supname = 'sup='+str(sup)+'_thre='+str(thre)+'_delay='+str(delay)+'_'+'Ne='+str(Ne)+', p='+str(p)+', nsig='+str(nsig)+', mu0='+str(mu0)+', nsigg='+str(nsigg)
            resy = np.load(os.path.join(mother, supname,'run'+str(1)+'_tr'+str(rp)+'.npy'))
            s1s, s2s, sg1s, sg2s, noises = unpack(resy)
            
            for ne in range(Ne):
                s1, s2, sg1, sg2 = s1s[ne], s2s[ne], sg1s[ne], sg2s[ne]
                for r in range(M2):
                    t = T2[r]; idx = int(t)
                    
                    if cond1:
                        N_win[M+r].append(s1[idx])
                        G_win[M+r].append(sg1[idx])
                        N_lose[M+r].append(s2[idx])
                        G_lose[M+r].append(sg2[idx])
                    else:
                        N_win[M+r].append(s2[idx])
                        G_win[M+r].append(sg2[idx])
                        N_lose[M+r].append(s1[idx])
                        G_lose[M+r].append(sg1[idx])
                    
    print(np.mean(rts))
        
    aucs = []; aucgs = []
    for n in range(M+M2):
        
        # for excitatory
        dist1 = np.array(N_win[n])
        dist2 = np.array(N_lose[n])
            
        # calculating auc
        x = list(dist1)+list(dist2)
        y = [1]*len(dist1)+[0]*len(dist2)
        x = np.array(x).reshape((len(x),1))
        clf = LogisticRegression(solver="liblinear", random_state=0).fit(x, y)
        auc = roc_auc_score(y, clf.predict_proba(x)[:, 1])
        aucs.append(auc)
        
        # for inhibitory
        dist1 = np.array(G_win[n])
        dist2 = np.array(G_lose[n])
            
        # calculating auc
        x = list(dist1)+list(dist2)
        y = [1]*len(dist1)+[0]*len(dist2)
        x = np.array(x).reshape((len(x),1))
        clf = LogisticRegression(solver="liblinear", random_state=0).fit(x, y)
        auc = roc_auc_score(y, clf.predict_proba(x)[:, 1])
        aucgs.append(auc)
        
        # for plotting roc
        coors = []
        for cthre in np.logspace(np.log10(0.1),np.log10(0.7), 10):
            TP = len([d for d in dist1 if d > cthre])/len(dist1) # true positive
            FP = len([d for d in dist2 if d > cthre])/len(dist1) # false positive
            coors.append((FP,TP))
        
    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(111)
    plt.plot(list(T)+list(T2), np.array(aucs)-0.5, 'b')
    plt.plot(list(T)+list(T2), np.array(aucgs)-0.5, 'r')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.xlabel('Time relative to choice'); plt.ylabel('|AUC-0.5|')
    plt.tight_layout()
    plt.savefig(os.path.join(mother,'thre='+str(thre)+'_'+fname+'-sup.png'), dpi=300); plt.clf()

        
    
                
        
        