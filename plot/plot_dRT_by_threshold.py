'''
Plot Delta_RT as a function of threshold.
'''

import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
from scipy.stats import beta, gamma

rnum = 1
mother = os.getcwd()
agvals = [[0.5,0.6,0.7,1,1], [0.5,1,1,0,2], [0.75,1,1,0,2], [0.83,1,1,0,2]]
color = ['grey']+list(sns.color_palette("hls", 3))

fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
for k in range(len(agvals)):
    gma, gma2, gma3, gma4, alpha = agvals[k]
    fname = 'SMA4e_gma='+str(gma)+', '+str(gma2)+', '+str(gma3)+', '+str(gma4)+', '+str(alpha)
    print(fname)
    traces = vsn.inquire_trace_npy(rnum, fname, 9999)
    
    rtcss = []
    for thre in np.arange(0.1,0.7,0.1):
        accs, rtcs = hpr.psych_func_complete([traces], thre=thre)
        rtcss.append(np.mean(rtcs))
    
    rtcss = np.array(rtcss)
    plt.plot(np.arange(0.1,0.7,0.1)[1:], rtcss[1:]-rtcss[:-1], color=color[k], marker='.')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('till this thre'); plt.ylabel('RT')
plt.tight_layout()
plt.savefig(os.path.join(mother,'time change.png'), dpi=600); plt.clf()