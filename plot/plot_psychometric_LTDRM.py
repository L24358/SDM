'''
Plot the psychometric function for LTDRM.
'''

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
agvals = list(np.load(os.path.join(bcs.datapath(),'files','agvals.npy')))
agvals = [agvals[x] for x in np.arange(2,38,7)]
color = sns.color_palette("hls", len(agvals)+1)
'''================================'''

# define parameters
mother = os.getcwd()
x = [0]+[np.log(c) for c in [3.2,6.4,12.8,25.6,51.2]]
az, cz, wz, asz, csz, wsz = [], [], [], [], [], []
for k in range(len(agvals)+1):
    if k == len(agvals): fname = 'SMAabv3Q'+'_gma='+str(0.5)+', '+str(0.6)+', '+str(0.7)+', '+str(1)+', '+str(1)
    else: fname = 'agval='+str(round(agvals[k],5))
    traces = []
    for rnum in range(1,7): traces.append(vsn.inquire_trace_npy(rnum, fname, 9999))
    accs, rtcs, rtws, accstd, rtcstd, rtwstd = hpr.psych_func(traces, thre=0.35)
    az.append(accs); asz.append(accstd)
    cz.append(rtcs); csz.append(rtcstd)
    wz.append(rtws); wsz.append(rtwstd)
    
# plot accuracy
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
for g in range(len(agvals)):
    plt.errorbar(x, az[g], yerr=asz[g], ls='-', color=color[g],\
                  marker='.', linewidth=1, markersize=1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_major_locator(ticker.FixedLocator(x))
plt.xlabel('coherence (%)'); plt.ylabel('accuracy')
plt.tight_layout()
plt.savefig(os.path.join(mother,'accuracy_UMA.png'), dpi=600); plt.clf()

# plot legend
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
for g in range(len(agvals)):
    plt.plot([0],[0], label='ag = '+str(round(agvals[g],5)), color=color[g])
plt.legend()
plt.savefig(os.path.join(mother,'legend_UMA.png'), dpi=600); plt.clf()
    