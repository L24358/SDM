import os
import imageio
import numpy as np
import matplotlib.pyplot as plt
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
import dynalysis.classes as clss

DT = int(0.3*100)
T = np.arange(-DT, DT, 5)
thre = 0.6
gmavals = [[0.5,0.6,0.7,1,1]]

b = clss.branch('gif_selectivity', os.getcwd()); b.mkdir()

for k in range(len(gmavals)):
    gma, gma2, gma3, gma4, alpha = gmavals[k]
    fname = 'SMA4e_gma='+str(gma)+', '+str(gma2)+', '+str(gma3)+', '+str(gma4)+', '+str(alpha)
    traces = vsn.inquire_trace_npy(1, fname, 9999)
    
    resps_1, resps_2 = [], []
    for trace in traces:
        s1, s2 = trace[:2]
        rt, win = hpr.accuracy(s1, s2, thre=thre)
        
        resps = []
        for t in T:
            idx = int(rt*100 + t + 50)
            if not np.isnan(rt) and (idx > 0) and (idx < 300):
                resp = s1[idx]
                resps.append(resp)
        
        if len(resps) == len(T):
            if win == 1: resps_1.append(resps)
            else: resps_2.append(resps)
        
    for t in range(len(T)):
        fig, ax = plt.subplots(figsize=(10,8))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.xlabel('Responses', fontsize=30); plt.ylabel('pdf', fontsize=30)
        dist1 = np.array(resps_1)[:,t]
        dist2 = np.array(resps_2)[:,t]
        plt.hist(dist1, color='g', alpha=0.7)
        plt.hist(dist2, color='r', alpha=0.5)
        plt.xlim((0,0.7))
        plt.ylim((0,600))
        plt.xticks(fontsize=24); plt.yticks(fontsize=24)
        ax.text(0.45, 570, r"$\bf{t = " + '{:.2f}'.format(T[t]) + "}$", color='k', fontsize=30)
        plt.tight_layout()
        plt.savefig(os.path.join(b.pathlink, 't='+str(t)+'.png'), dpi=50)
        plt.close('all')
                
    images = []       
    for t in range(len(T)):
        images.append(imageio.imread(os.path.join(b.pathlink, 't='+str(t)+'.png')))
    imageio.mimsave(os.path.join(b.pathlink,'selectivity.gif'), images, duration=0.2)
    
                
        
        