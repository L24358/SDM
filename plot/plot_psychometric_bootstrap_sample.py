'''
Plot averaged psychometric function with single bootstrap sample.
'''

import os
import pickle
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution

'''=====Define parameters here====='''
gmavals = [[0.3,1,1,0,2]]
'''================================'''

def pick_and_mean(l):
    l = -np.array(l)+2
    k = np.random.choice(l, len(l))
    return np.mean(k)

def pick_and_mean2(l):
    k = np.random.choice(l, len(l))
    return np.mean(k)

def Quad(T, a, b, c):
    return [a*t**2+b*t+c for t in T]
    
def bounds():
    return [(-100,100), (-100,100), (-100,100)]

# define parameters
mother = os.getcwd()
coherence = [0,3.2,6.4,12.8,25.6,51.2]
x = np.linspace(0,np.log(51.2),num=20)[1:]
xexp = [0] + [np.exp(i) for i in x]
x = [0] + list(x)
xdot = [0]+[np.log(c) for c in coherence[1:]]
B = 1 # bootstrap sample number

# main
for k in range(len(gmavals)):
    bootaccs, bootrtcs, bootrtws = [], [], []
    gma, gma2, gma3, gma4, alpha = gmavals[k]
    fname = 'SMA4e_gma='+str(gma)+', '+str(gma2)+', '+str(gma3)+', '+str(gma4)+', '+str(alpha)
    for rnum in range(1,7):
        traces = vsn.inquire_trace_npy(rnum, fname, 9999)
        accs, rtcs = hpr.psych_func_complete([traces])
        acc, rtc = accs[0], rtcs[0]
        bootacc = [pick_and_mean(acc) for b in range(B)]
        bootrtc = [pick_and_mean2(rtc) for b in range(B)]
        bootaccs.append(bootacc)
        bootrtcs.append(bootrtc)
        
    bootaccs = np.transpose(np.array(bootaccs))
    for b in range(B):
        popt, pcov = curve_fit(hpr.psyched, coherence, bootaccs[b])
        y = [hpr.psyched(c,*popt) for c in xexp]
        
    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(111)
    plt.plot(x, y, 'k--', label='fitted') # averaged result
    plt.plot(xdot, bootaccs[0], 'b.', label='bootstrap sample') # single bootstrap sample
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    xaxis = [0]+[np.log(c) for c in [3.2,6.4,12.8,25.6,51.2]]
    ax.xaxis.set_major_locator(ticker.FixedLocator(xaxis))
    plt.xlabel('coherence'); plt.ylabel('accuracy')
    plt.tight_layout()
    plt.savefig(os.path.join(mother,'fit_accuracy.png'), dpi=300)
        
    bootrtcs = np.transpose(np.array(bootrtcs))
    for b in range(B):
        def SSE(pm):
            warnings.filterwarnings('ignore')
            val = Quad(coherence, *pm)
            return np.sum((np.array(bootrtcs[b])-np.array(val))**2)
        popt = differential_evolution(SSE, bounds=bounds(), seed=3)['x']
        y = Quad(xexp,*popt)
        
    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(111)
    plt.plot(x, y, 'k--', label='fitted')
    plt.plot(xdot, bootrtcs[0], 'b.', label='bootstrap sample')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    xaxis = [0]+[np.log(c) for c in [3.2,6.4,12.8,25.6,51.2]]
    ax.xaxis.set_major_locator(ticker.FixedLocator(xaxis))
    plt.xlabel('coherence'); plt.ylabel('RT, correct')
    plt.tight_layout()
    plt.savefig(os.path.join(mother,'fit_rtc.png'), dpi=300)
        
        
    
    