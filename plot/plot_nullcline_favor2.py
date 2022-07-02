'''
Plot nullcline in which population 2 is favored.

* Notes:
- Make smoother: https://stackoverflow.com/questions/5283649/plot-smooth-line-with-pyplot
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

fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)

res = hpr.get_nullcline_u2(12.8, 0.6, x0, params, 1)
s1 = res[0]; s2 = res[1][0]
s1, s2 = hpr.sort_nullcline(s1, s2)
plt.plot(s1,s2,'k',linewidth=1.5)

res = hpr.get_nullcline_u2(12.8, 0.6, x0, params, 2)
s1 = res[0]; s2 = res[1][0]
s1, s2 = hpr.sort_nullcline(s1, s2)
plt.plot(s2,s1,'k',linewidth=1.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('S\u2081'); plt.ylabel('S\u2082')
plt.tight_layout()
plt.savefig(os.path.join(mother,'nullcline_favor2.png'), dpi=600); plt.clf()