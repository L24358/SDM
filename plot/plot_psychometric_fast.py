import os
import pickle
import warnings
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution

'''=====Define parameters here====='''
files = ['SMA4e_fast=10_gma=0.3, 1, 1, 0, 2',\
         'SMA4e_fast=10_gma=0.5, 1, 1, 0, 2',\
         'SMA4e_fast=10_gma=0.83, 1, 1, 0, 2']
color = sns.color_palette("hls", 5)
color = [color[0],color[2],color[4]]
'''================================'''

#define parameters
mother = os.getcwd()
x = [0]+[np.log(c) for c in [3.2,6.4,12.8,25.6,51.2]]
az, cz, wz, asz, csz, wsz = [], [], [], [], [], []
for k in range(len(files)):
    fname = files[k]
    traces = []
    for rnum in range(1,7): traces.append(vsn.inquire_trace_npy(rnum, fname, 9999))
    accs, rtcs, rtws, accstd, rtcstd, rtwstd = hpr.psych_func(traces, thre=0.35)
    az.append(accs); asz.append(accstd)
    cz.append(rtcs); csz.append(rtcstd)
    wz.append(rtws); wsz.append(rtwstd)
    
#plot accuracy
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
for g in range(len(az)):
    plt.errorbar(x, az[g], yerr=asz[g], ls='-', color=color[g],\
                  marker='.', linewidth=1, markersize=1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_major_locator(ticker.FixedLocator(x))
plt.xlabel('coherence (%)'); plt.ylabel('accuracy')
plt.tight_layout()
plt.savefig(os.path.join(mother,'accuracy_UMA_fast.png'), dpi=600); plt.clf()

#plot reaction time
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
for g in range(len(az)):
    plt.errorbar(x, cz[g], yerr=csz[g], ls='-', color=color[g],\
                  marker='.', linewidth=1, markersize=1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_major_locator(ticker.FixedLocator(x))
plt.xlabel('coherence (%)'); plt.ylabel('reaction time (s)')
plt.tight_layout()
plt.savefig(os.path.join(mother,'correct_UMA_fast.png'), dpi=600); plt.clf()

#plot reaction time
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
for g in range(len(az)):
    plt.errorbar(x, wz[g], yerr=wsz[g], ls='-', color=color[g],\
                  marker='.', linewidth=1, markersize=1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_major_locator(ticker.FixedLocator(x))
plt.xlabel('coherence (%)'); plt.ylabel('reaction time (s)')
plt.tight_layout()
plt.savefig(os.path.join(mother,'wrong_UMA_fast.png'), dpi=600); plt.clf()

#plot legend
for g in range(len(files)):
    plt.plot([0],[0], label='\u03B3 = '+str(files[g]), color=color[g])
plt.legend()
plt.savefig(os.path.join(mother,'legend_UMA_fast.png'), dpi=600); plt.clf()
    