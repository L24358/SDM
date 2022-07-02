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
nsig = 3.5
thre = 0.6
gmavals = [[0.5,0.6,0.7,1,1]]
color = sns.color_palette("hls", 5)
'''================================'''

#define parameters
mother = os.getcwd()
x = [np.log10(c) for c in [3.2,6.4,12.8,25.6,51.2]]
az, cz, wz, asz, csz, wsz = [], [], [], [], [], []
for k in range(len(gmavals)):
    gma, gma2, gma3, gma4, alpha = gmavals[k]
    fname = 'SMA4e'+'_nsig='+str(round(nsig,3))+'_gma='+str(gma)+', '+str(gma2)+', '+str(gma3)+', '+str(gma4)+', '+str(alpha)
    traces = []
    for rnum in range(2,7): traces.append(vsn.inquire_trace_npy(rnum, fname, 9999))
    accs, rtcs, rtws, accstd, rtcstd, rtwstd = hpr.psych_func(traces, thre=thre)
    az.append(accs); asz.append(accstd)
    cz.append(rtcs); csz.append(rtcstd)
    wz.append(rtws); wsz.append(rtwstd)
    
#plot accuracy
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
for g in range(len(gmavals)):
    plt.errorbar(x, az[g], yerr=asz[g], ls='-', color=color[g],\
                  marker='.', linewidth=1, markersize=1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_major_locator(ticker.FixedLocator(x))
plt.xlabel('coherence (%)'); plt.ylabel('accuracy')
plt.tight_layout()
plt.show(); plt.clf()
# plt.savefig(os.path.join(mother,'accuracy_UMA.png'), dpi=600); plt.clf()

#plot reaction time
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
for g in range(len(az)):
    plt.errorbar(x, cz[g], yerr=csz[g], ls='-', color=color[g],\
                  marker='.', linewidth=1, markersize=1)
    plt.errorbar(x, wz[g], yerr=wsz[g], ls='-', color=color[g],\
                  marker='.', linewidth=1, markersize=1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_major_locator(ticker.FixedLocator(x))
plt.xlabel('coherence (%)'); plt.ylabel('reaction time (s)')
plt.tight_layout()
plt.show(); plt.clf()
# plt.savefig(os.path.join(mother,'correct_UMA.png'), dpi=600); plt.clf()
