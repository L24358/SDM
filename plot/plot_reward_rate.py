'''
Plot reward rate as a function of selectivity.
'''

import os
import seaborn as sns
import aadm.helper as hpr
import aadm.vision as vsn
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

iii = 0
penalty = 0
thre = 0.6
mother = 'C:\\Users\\belle\\Desktop'

gma, gma3, gma4, alpha = 0.5, 0.7, 1, 1
gma2s = [1, 0.9, 0.8, 0.7, 0.6, 0.55, 0.5]
mrates = []
for gma2 in gma2s:
    fname = 'SMA4e_gma='+str(gma)+', '+str(gma2)+', '+str(gma3)+', '+str(gma4)+', '+str(alpha)
    
    rates = []
    for rnum in range(1,7):
        traces = vsn.inquire_trace_npy(rnum, fname, 9999)
        accs, rtcs, rtws = hpr.psych_func_complete2([traces], thre=thre)
        
        N = len(accs[0])
        N_w = len(rtws[0])
        reward = accs[0].count(1)
        rtc_total = sum(rtcs[0])
        rtw_total = sum(rtws[0])
        if rnum!=1: add = 0
        else: add = N_w*penalty
        rate = reward/(rtc_total + rtw_total + iii*(N-1) + add)
        rates.append(rate)
        
    mrate = np.mean(rates)
    mrates.append(mrate)
    
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.plot(gma2s, mrates, 'ko', linestyle='-')
ax.set_xlim(1.05, 0.45)
plt.xlabel('selectivity (1-gma2)'); plt.ylabel('reward rate')
plt.tight_layout()
plt.savefig(os.path.join(mother,'reward.png'), dpi=600); plt.clf()
        