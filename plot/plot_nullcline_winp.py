'''
Plot nullcline with input.
'''

import os
import aadm.basics as bcs
import aadm.helper as hpr
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

mother = os.getcwd()
x0 = bcs.init_x0(2, d=0.1)
params = hpr.get_params_u1(0.3)
res = hpr.get_nullcline_u2(0, 0.6, x0, params, 1)
s1 = res[0]; s2 = res[1][0]
s1, s2 = hpr.sort_nullcline(s1, s2)

# traces = hpr.inquire_trace(1, 'UMA_gma='+str(0.3), 10, 0)
# s1r, s2r, _, _, _ = traces[0]
# s1w, s2w, _, _, _ = traces[6]

fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
plt.plot(s1,s2,'k',linewidth=1.5)
plt.plot(s2,s1,'k',linewidth=1.5)
# plt.plot(s1r, s2r, 'g.', markersize=1)
# plt.plot(s1w, s2w, 'r.', markersize=1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('S\u2081'); plt.ylabel('S\u2082')
plt.tight_layout()
plt.savefig(os.path.join(mother,'nullcline_inp.png'), dpi=600); plt.clf()