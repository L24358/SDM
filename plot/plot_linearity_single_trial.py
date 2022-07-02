'''
Plot linear relation for a single trial of UM, SM plus; and seemingly linear relation of a single trial of SM minus.
'''

import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn

mother = os.getcwd()
color = sns.color_palette("hls", 5)
galists = [[0.1,1,1,0,2],[0.3,1,1,0,2],[0.5,1,1,0,2],\
            [0.75,1,1,0,2], [0.83,1,1,0,2]]
cnt = 0; endt = 100; sd = 1
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)

# plot UM
for ga in galists:
    print(ga)
    gma, gma2, gma3, gma4, alpha = ga
    fname = 'SMA4e_gma='+str(bcs.rn(gma))+', '+str(bcs.rn(gma2))+', '+str(bcs.rn(gma3))+', '+str(bcs.rn(gma4))+', '+str(bcs.rn(alpha))
    trace = vsn.inquire_trace_npy(1, fname, 10)[sd]
    s1, s2, sg1, sg2, _, _ = trace
    splus = (s1+s2)/2; sgplus = (sg1+sg2)/2
    sminus = (s1-s2)/2; sgminus = (sg1-sg2)/2
    ax.plot(splus[50:endt], sgplus[50:endt], color=color[cnt]); cnt += 1
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('sg'); plt.ylabel('splus')
plt.tight_layout()
plt.savefig(os.path.join(mother,'linearity_UM.png'), dpi=600); plt.clf()
        
# plot SM
galists = [[0.3,0.6,0.5,1,1],[0.3,0.6,0.6,1,1],[0.3,0.6,0.7,1,1],\
            [0.3,0.6,0.8,1,1]]
cnt = 0
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
for ga in galists:
    print(ga)
    gma, gma2, gma3, gma4, alpha = ga
    fname = 'SMA4e_gma='+str(bcs.rn(gma))+', '+str(bcs.rn(gma2))+', '+str(bcs.rn(gma3))+', '+str(bcs.rn(gma4))+', '+str(bcs.rn(alpha))
    trace = vsn.inquire_trace_npy(1, fname, 10)[sd]
    s1, s2, sg1, sg2, _, _ = trace
    splus = (s1+s2)/2; sgplus = (sg1+sg2)/2
    sminus = (s1-s2)/2; sgminus = (sg1-sg2)/2
    ax.plot(splus[50:endt], sgplus[50:endt], color=color[cnt]); cnt += 1
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('sg'); plt.ylabel('splus')
plt.tight_layout()
plt.savefig(os.path.join(mother,'linearity_SM.png'), dpi=600); plt.clf()
 
# plot SM, minus direction
galists = [[0.3,0.6,0.5,1,1],[0.3,0.6,0.6,1,1],[0.3,0.6,0.7,1,1],\
            [0.3,0.6,0.8,1,1]]
cnt = 0
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
for ga in galists:
    print(ga)
    gma, gma2, gma3, gma4, alpha = ga
    fname = 'SMA4e_gma='+str(bcs.rn(gma))+', '+str(bcs.rn(gma2))+', '+str(bcs.rn(gma3))+', '+str(bcs.rn(gma4))+', '+str(bcs.rn(alpha))
    trace = vsn.inquire_trace_npy(1, fname, 10)[sd]
    s1, s2, sg1, sg2, _, _ = trace
    splus = (s1+s2)/2; sgplus = (sg1+sg2)/2
    sminus = (s1-s2)/2; sgminus = (sg1-sg2)/2
    ax.plot(sminus[50:endt], sgminus[50:endt], color=color[cnt]); cnt += 1
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('sg'); plt.ylabel('splus')
plt.tight_layout()
plt.savefig(os.path.join(mother,'linearity_SM_minus.png'), dpi=600); plt.clf()
 