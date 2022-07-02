'''
Plot velocity magnitude as a function of time.
'''

import os
import seaborn as sns
import aadm.helper as hpr
import aadm.vision as vsn
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

gmavals = [[0.3,1,1,0,2],[0.75,1,1,0,2],[0.83,1,1,0,2]]

mother =os.getcwd()
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = sns.color_palette("hls", 5)
colind = [0,2,4]

fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
for g in range(len(gmavals)):
    gma, gma2, gma3, gma4, alpha = gmavals[g]
    fname = 'SMA4e_gma='+str(gma)+', '+str(gma2)+', '+str(gma3)+', '+str(gma4)+', '+str(alpha)
    traces = vsn.inquire_trace_npy(1, fname, 9999)
    vels = hpr.get_flow_bunch(traces)
    mean = np.mean(vels, axis=0); std = np.std(vels, axis=0)
    plt.plot(np.arange(0,3,0.01)[50:150], mean[50:150], color=colors[colind[g]],\
                     linewidth=1)
    plt.fill_between(np.arange(0,3,0.01)[50:150], mean[50:150]-std[50:150], mean[50:150]+std[50:150], color=colors[colind[g]],\
                     alpha=0.2, lw=0)
    
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('t (s)'); plt.ylabel('lag')
plt.tight_layout()
plt.savefig(os.path.join(mother,'vel.png'), dpi=600); plt.clf()
