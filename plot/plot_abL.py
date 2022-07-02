'''
Fits aL(t), bL(t) values.
'''
import os
import gc
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
import dynalysis.classes as clss
from mpl_toolkits import mplot3d
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
from matplotlib import cm

gmaval = [0.5,0.6,0.7,1,1]
gma, gma2, gma3, gma4, alpha = gmaval
params = hpr.get_params_s2(gma, gma2, gma3, gma4, alpha)
color = sns.color_palette("hls", 3)

ass, bss = [], []
for ss in np.arange(0.1,0.35,0.1): # different starting points

    s1, s2, sg1, sg2 = vsn.sim_SM_nonoise(ss, params, d = 0.05)
    sminus = (s1-s2)/2; sgminus = (sg1-sg2)/2
    period = []; period2 = []

    for t0 in np.arange(0,9500,250): # different time periods
        popt, pcov = curve_fit(hpr.affine, sminus[t0:t0+250], sgminus[t0:t0+250])
        period.append(popt[1]); period2.append(popt[0])
    ass.append(period2); bss.append(period)
    
# plot aL
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
for p in range(len(ass)):
    plt.plot(np.arange(0,9500,250)/10000, ass[p], color=color[p])
plt.xlabel('t'); plt.ylabel('ag')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig('ag_evolution.png', dpi=600)

# plot bL
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
for p in range(len(bss)):
    plt.plot(np.arange(0,9500,250)/10000, bss[p], color=color[p])
plt.plot([0,9250/10000],[0,0],color='grey',linewidth=0.8, linestyle='--')
plt.xlabel('t'); plt.ylabel('bg')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig('bg_evolution.png', dpi=600)

# plot legend
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
splus = [0.1,0.2,0.3]
for p in range(len(bss)):
    plt.plot([0,0], [0,0], color=color[p], label=r'$S_{+}$'+' = '+str(splus[p]))
plt.plot([0,0],[0,0],color='grey',linewidth=0.8, linestyle='--', label='y = 0')
plt.xlabel('t'); plt.ylabel('bg')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.legend()
plt.tight_layout()
plt.savefig('legend_bg_evolution.png', dpi=300)