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
    plt.plot(s1, s2, 'g',linewidth=0.8)
    
traces = vsn.inquire_trace_npy(1, 'SMA4e_gma=0.5, 1, 1, 0, 2', 10, mother=bcs.datapath())
for trace in traces:
    s1, s2, _, _, _, _ = trace
    plt.plot(s1, s2, 'r',linewidth=0.8)
    
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('t (s)'); plt.ylabel('S_NMDA')
plt.xlim((0,0.7)); plt.ylim((0,0.7))
plt.tight_layout()
plt.savefig(os.path.join(mother,'single_trace_SNMDA_2d.png'), dpi=600); plt.clf()