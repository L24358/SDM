import os
import aadm.basics as bcs
import aadm.helper as hpr
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

mother = os.getcwd()
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)

x0 = bcs.init_x0(1, d=0.01)
res = hpr.get_nullcline_two(0, 0, x0, [0.2609,0.0497,0.3255], 1)
s1 = res[0]; s2 = res[1][0]
s1, s2 = hpr.sort_nullcline(s1, s2)
plt.plot(s1,s2,'k',linewidth=1.5)
plt.plot(s2,s1,'k',linewidth=1.5)

res = hpr.get_nullcline_two(0, 0, x0, [0.262, 0.02 , 0.324], 1, minn=-0.1)
s1 = res[0]; s2 = res[1][0]
s1, s2 = hpr.sort_nullcline(s1, s2)
plt.plot(s1,s2,color='gray',linewidth=1.2,linestyle='--')
plt.plot(s2,s1,color='gray',linewidth=1.2,linestyle='--')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('s1'); plt.ylabel('s2')
plt.tight_layout()
plt.savefig(os.path.join(mother,'nullcline_noinp.png'), dpi=600); plt.clf()