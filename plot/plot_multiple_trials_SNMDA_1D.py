import os
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

mother =os.getcwd()
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)

traces = vsn.inquire_trace_npy(1, 'SMA4e_gma=0.5, 0.6, 0.7, 1, 1', 10, mother=bcs.datapath())
for trace in traces:
    s1, s2, _, _, _, _ = trace
    if np.mean(s1)>np.mean(s2): x, y = s1, s2
    else: x, y = s2, s1
    plt.plot(np.arange(0,3,0.01),x, 'k', linewidth=0.7)
    plt.plot(np.arange(0,3,0.01),y, 'k--', linewidth=0.7)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('t (s)'); plt.ylabel('S_NMDA')
plt.xlim((-0.1,3)); plt.ylim((-0.1,0.7))
plt.tight_layout()
plt.savefig(os.path.join(mother,'single_trace_SNMDA.png'), dpi=600); plt.clf()

fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
plt.plot([0], [0], color='k', label='winner')
plt.plot([0], [0], color='grey', label='loser')
plt.legend()
plt.savefig(os.path.join(mother,'legend_SNMDA_traces.png'), dpi=300); plt.clf()