'''
Fit reaction time to gamma distribution.
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
thre = 0.6
mother = os.getcwd()
agvals = [[0.5,0.6,0.7,1,1], [0.5,1,1,0,2], [0.75,1,1,0,2], [0.83,1,1,0,2]]
color = ['grey']+list(sns.color_palette("hls", 3))

# main
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)

for k in range(len(agvals)):
    gma, gma2, gma3, gma4, alpha = agvals[k]
    fname = 'SMA4e_gma='+str(gma)+', '+str(gma2)+', '+str(gma3)+', '+str(gma4)+', '+str(alpha)
    print(fname)
    traces = vsn.inquire_trace_npy(rnum, fname, 9999)
    accs, rtcs = hpr.psych_func_complete([traces], thre=thre)
        
    a, _, loc = gamma.fit(rtcs) # fit reaction time to gamma distribution
    x = np.linspace(gamma.ppf(0.01, a, scale=loc),gamma.ppf(0.99, a, scale=loc), 100)
    pdf = gamma.cdf(x, a, scale=loc)
    plt.plot(x, pdf, label=str(k), color=color[k])

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('RT'); plt.ylabel('scaled pdf')
plt.tight_layout()
plt.savefig(os.path.join(mother,'distribution.png'), dpi=600); plt.clf()

# plot legend
lbls = ['SM', 'UM, \u03B3=0.5',  'UM, \u03B3=0.75', 'UM, \u03B3=0.83']
for k in range(len(agvals)):
    plt.plot([0], [0], label=lbls[k], color=color[k])
plt.legend()
plt.savefig(os.path.join(mother,'legend_distribution.png'), dpi=300); plt.clf()