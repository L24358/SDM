'''
Plot psychometric function in the style of Wong and Wang's paper for comparison purposes.
'''

import os
import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

'''=====Define parameters here====='''
gmavals = [[0.5,0.6,0.7,1,1]] 
thres = [0.6]
'''================================'''

nsig = 4
Q = '_ecomp_nsig='+str(nsig)
rtname = 'rtc'
mother = os.getcwd()
coherence = [0,3.2,6.4,12.8,25.6,51.2]
xaxis = [0] + [np.log(c) for c in coherence[1:]]

for thre in thres:
    print(thre)
    folder = 'files_4e_thre='+str(thre)
    storages = []
    files = ['apall'+Q,'bpall'+Q]
    for f in range(len(files)):
        infofile = os.path.join(bcs.datapath(), folder, files[f]+'.p')
        storages.append(pickle.load(open(infofile, 'rb')))
    apall,bpall = storages
    
    for k in range(len(gmavals)):
        gma, gma2, gma3, gma4, alpha = gmavals[k] 
        fname = 'SMA4e'+'_nsig='+str(round(nsig,3))+'_gma='+str(gma)+', '+str(gma2)+', '+str(gma3)+', '+str(gma4)+', '+str(alpha)
        aps, bps = apall[tuple(gmavals[k])], bpall[tuple(gmavals[k])]
        print(np.mean(aps))
        print(np.mean(bps))
        
        traces = [vsn.inquire_trace_npy(run, fname, 9999) for run in range(1,7)]
        
        accs, rtcs, rtws = hpr.psych_func(traces, suppress=True, thre=thre)[:3]
        print(rtcs)
        print(rtws)
        print('')
        
        fig = plt.figure(figsize=(4,3))
        ax = fig.add_subplot(111)
        plt.plot(xaxis, accs, 'k', marker='o', label='correct', fillstyle='full')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.xaxis.set_major_locator(ticker.FixedLocator(xaxis))
        plt.xlabel('coherence'); plt.ylabel('acc')
        plt.tight_layout()
        plt.savefig(os.path.join(mother,'acc_nsig='+str(nsig)+'.png'), dpi=600); plt.clf()
        
        fig = plt.figure(figsize=(4,3))
        ax = fig.add_subplot(111)
        plt.plot(xaxis[1:], rtcs[1:], 'k', marker='o', label='correct', fillstyle='full')
        plt.plot(xaxis[1:], rtws[1:], 'k', marker='o', label='incorrect', fillstyle='none')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.xaxis.set_major_locator(ticker.FixedLocator(xaxis))
        plt.xlabel('coherence'); plt.ylabel('RT')
        plt.tight_layout()
        plt.savefig(os.path.join(mother,'rt_nsig='+str(nsig)+'.png'), dpi=600); plt.clf()
