import os
import aadm.helper as hpr
import aadm.vision as vsn
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

mother =os.getcwd()
trace = vsn.inquire_trace_npy(1, 'SMA4e_gma=0.3, 1, 1, 0, 2', 1)[0]
s1, s2, _, _, _, _ = trace

fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
plt.plot(np.arange(0,3,0.01),s1, 'r',linewidth=1.5)
plt.plot(np.arange(0,3,0.01),s2, 'g',linewidth=1.5)
plt.plot([0.5,0.5],[-0.1,0.2],color='grey',linestyle='--',linewidth=1.5)
plt.plot([-0.1,2.5],[0.35,0.35],color='grey',linestyle='--',linewidth=1.5)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('t (s)'); plt.ylabel('S_NMDA')
plt.xlim((-0.1,3)); plt.ylim((-0.1,0.7))
plt.tight_layout()
plt.savefig(os.path.join(mother,'single_trace.png'), dpi=600); plt.clf()