'''
Plot accuracy and reaction time parameters in the phase space of gamma2 and gamma3.
'''

import os
import pickle
import warnings
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
import dynalysis.classes as clss
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution

'''=====Define parameters here====='''
gma2s =  gma3s = np.arange(0.3,1,0.1) #1.1
thre = 0.6
'''================================'''

mother = os.getcwd()
storages = []
folder = 'files_4e_thre='+str(thre)
rtname, Q = 'fst', ''
files = [rtname+'all'+Q,rtname+'std'+Q,'ap'+Q,'apstd'+Q]
for f in range(len(files)):
    infofile = os.path.join(bcs.datapath(), folder, files[f]+'.p')
    storages.append(pickle.load(open(infofile, 'rb')))
rtc,rtcstd,ap,apstd = storages

def get_distance(dx, dy, scale):
    res = (dx*scale)**2+dy**2
    return np.log10(res**0.5)

#Scale by average
aps, rtcs = [], []
for gma2 in gma2s:
    temp = []
    if round(gma2,1) != 1.0: gma2_rev = round(gma2,1)
    else: gma2_rev = 1
    for gma3 in gma3s:
        if round(gma3,1) != 1.0: gma3_rev = round(gma3,1)
        else: gma3_rev = 1
        if (0.5,gma2_rev,gma3_rev,1,1) in rtc.keys():
            aps.append(ap[(0.5,gma2_rev,gma3_rev,1,1)])
            rtcs.append(rtc[(0.5,gma2_rev,gma3_rev,1,1)])
            
scale = np.mean(aps)/np.mean(rtcs)
dx_um = rtc[(0.5,1,1,1,1)]
dy_um = ap[(0.5,1,1,1,1)]
scale = dy_um/dx_um

#Get matrix
matrix, matrix2 = [], []
for gma2 in gma2s:
    temp, temp2 = [], []
    if round(gma2,1) != 1.0: gma2_rev = round(gma2,1)
    else: gma2_rev = 1
    
    for gma3 in gma3s:
        if round(gma3,1) != 1.0: gma3_rev = round(gma3,1)
        else: gma3_rev = 1
            
        if (0.5,gma2_rev,gma3_rev,1,1) not in rtc.keys(): temp.append(.7); temp2.append(17)
        else:
            dx = rtc[(0.5,gma2_rev,gma3_rev,1,1)]
            dy = ap[(0.5,gma2_rev,gma3_rev,1,1)]
            temp.append(dx); temp2.append(dy)

    matrix.append(temp); matrix2.append(temp2)
    
ticks = [1,3,5]
fig = plt.figure(figsize=(4,3))
ax = sns.heatmap(matrix, cmap="YlGnBu")
ax.set_xticks(ticks)
ax.set_xticklabels([round(gma3s[i],2) for i in ticks])
ax.set_yticks(ticks)
ax.set_yticklabels([round(gma2s[i],2) for i in ticks])
ax.invert_yaxis()
plt.xlabel('gma3'); plt.ylabel('gma2')
plt.tight_layout()
plt.savefig(os.path.join(mother,'phase_selective_rtc.png'), dpi=600); plt.clf()

fig = plt.figure(figsize=(4,3))
ax = sns.heatmap(matrix2, cmap="YlGnBu")
ax.set_xticks(ticks)
ax.set_xticklabels([round(gma3s[i],2) for i in ticks])
ax.set_yticks(ticks)
ax.set_yticklabels([round(gma2s[i],2) for i in ticks])
ax.invert_yaxis()
plt.xlabel('gma3'); plt.ylabel('gma2')
plt.tight_layout()
plt.savefig(os.path.join(mother,'phase_selective_ap.png'), dpi=600); plt.clf()