'''
Plot deviation (lag) as a function of time.
'''

import os
import seaborn as sns
import aadm.helper as hpr
import aadm.vision as vsn
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

mother =os.getcwd()
gmalist = [0.3,0.75,0.83]
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = sns.color_palette("hls", 5)
colind = [0,2,4]

fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
for g in range(len(gmalist)):
    params = hpr.get_params_u1(gmalist[g])
    fname = 'SMA4e_gma='+str(gmalist[g])+', 1, 1, 0, 2'
    traces = vsn.inquire_trace_npy(1, fname, 9999)
    lags = hpr.slowpoke(traces, params)[-1]
    mean = np.mean(lags, axis=0); std = np.std(lags, axis=0)
    plt.plot(np.arange(0,3,0.01), mean, color=colors[colind[g]],\
                     linewidth=1)
    plt.fill_between(np.arange(0,3,0.01), mean-std, mean+std, color=colors[colind[g]],\
                     alpha=0.2, lw=0)
    
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('t (s)'); plt.ylabel('lag')
plt.tight_layout()
plt.savefig(os.path.join(mother,'lag.png'), dpi=600); plt.clf()

#plot legend
for g in range(len(gmalist)):
    plt.plot([0],[0], label='\u03B3 = '+str(gmalist[g]), color=colors[colind[g]])
plt.legend()
plt.savefig(os.path.join(mother,'legend_lag.png'), dpi=600); plt.clf()