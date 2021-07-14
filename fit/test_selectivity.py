import os
import warnings
import numpy as np
import matplotlib.pyplot as plt
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

T = np.arange(-60, 60, 5)
M = len(T)
thre = 0.17
gmavals = [(50,1,5,0.0052,2.5)] 
mother = os.getcwd()

def unpack(v):
    return v[:Ne], v[Ne:2*Ne], v[2*Ne:2*Ne+Ni], v[2*Ne+Ni:N], v[N:]

for k in range(len(gmavals)):
    Ne, p, nsig, mu0, nsigg = gmavals[k]
    Ni = Ne; N = 2*(Ne + Ni)
    fname = 'Ne='+str(Ne)+', p='+str(p)+', nsig='+str(nsig)+', mu0='+str(mu0)+', nsigg='+str(nsigg)
    
    N_win, N_lose = {}, {}
    G_win, G_lose = {}, {}
    for n in range(M):
        N_win[n] = []; N_lose[n] = []
        G_win[n] = []; G_lose[n] = []
    
    rts = []
    for rp in range(228):
        resy = np.load(os.path.join(fname,'run'+str(1)+'_tr'+str(rp)+'.npy'))
        s1s, s2s, sg1s, sg2s, noises = unpack(resy)
        mean_s1, mean_s2 = np.mean(s1s,axis=0), np.mean(s2s,axis=0)
        rt, win = hpr.accuracy(mean_s1, mean_s2, thre=thre)
        
        if rt != None:
            rts.append(rt)
            cond1 = (win == 1)
            for ne in range(Ne):
                s1, s2, sg1, sg2 = s1s[ne], s2s[ne], sg1s[ne], sg2s[ne]
                for n in range(M):
                    t = T[n]; idx = int(rt*100 + t + 50)
                    cond2 = (idx >= 0) and (idx < 300)
                    
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
    print(np.mean(rts))
        
    aucs = []; aucgs = []
    for n in range(M):
        
        #For excitatory
        dist1 = np.array(N_win[n])
        dist2 = np.array(N_lose[n])
            
        #Calculating auc
        x = list(dist1)+list(dist2)
        y = [1]*len(dist1)+[0]*len(dist2)
        x = np.array(x).reshape((len(x),1))
        clf = LogisticRegression(solver="liblinear", random_state=0).fit(x, y)
        auc = roc_auc_score(y, clf.predict_proba(x)[:, 1])
        aucs.append(auc)
        
        #For inhibitory
        dist1 = np.array(G_win[n])
        dist2 = np.array(G_lose[n])
            
        #Calculating auc
        x = list(dist1)+list(dist2)
        y = [1]*len(dist1)+[0]*len(dist2)
        x = np.array(x).reshape((len(x),1))
        clf = LogisticRegression(solver="liblinear", random_state=0).fit(x, y)
        auc = roc_auc_score(y, clf.predict_proba(x)[:, 1])
        aucgs.append(auc)
        
        #For plotting roc
        coors = []
        for cthre in np.logspace(np.log10(0.1),np.log10(0.7), 10):
            TP = len([d for d in dist1 if d > cthre])/len(dist1)
            FP = len([d for d in dist2 if d > cthre])/len(dist1)
            coors.append((FP,TP))
        
    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(111)
    plt.plot(T, np.array(aucs)-0.5, 'b')
    plt.plot(T, np.array(aucgs)-0.5, 'r')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.xlabel('Time relative to choice'); plt.ylabel('|AUC-0.5|')
    plt.tight_layout()
    plt.show(); plt.clf()
    # plt.savefig(os.path.join(mother,'thre='+str(thre)+'_'+fname+'.png'), dpi=300); plt.clf()

        
    
                
        
        