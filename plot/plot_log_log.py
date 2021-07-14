import os
import gc
import time
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
import dynalysis.basics as dbcs
import dynalysis.classes as clss
from math import sqrt
from scipy.optimize import curve_fit

'''=====Define parameters here====='''
gmavals = [[0.3,0.6,0.5,1,1], [0.3,0.6,0.7,1,1], [0.3,0.6,0.7,1,1],\
           [0.3,0.6,0.8,1,1]] 
'''================================'''

mother = os.getcwd()
color = sns.color_palette('hls', 5)

fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
for k in range(len(gmavals)):
    
    gma, gma2, gma3, gma4, alpha = gmavals[k]
    fname = 'SMA4e_gma='+str(gma)+', '+str(gma2)+', '+str(gma3)+', '+str(gma4)+', '+str(alpha)

    print('Doing '+fname)
    params = hpr.get_params_s2(gma, gma2, gma3, gma4, alpha)
    gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params
    a, b = vsn.get_ab_SM(params)
    ag, bg = vsn.get_abg_SM(params)
    
    s1, s2, sg1, sg2 = vsn.sim_SM_nonoise(0.25, params, d=0.01)
    sminus = (s1-s2)/2; sgminus = (sg1-sg2)/2
    sminus = np.sign(np.mean(sminus))*sminus
    sgminus = np.sign(np.mean(sgminus))*sgminus
    poptl, pcovl = curve_fit(hpr.Quad, sminus, sgminus)
    linear = hpr.affine(sminus, poptl[1], poptl[2])
    residual = sgminus - linear
    popt, pcov = curve_fit(hpr.Quad, sminus, residual)
    
    plt.plot(np.log(sminus[5000:]), np.log(residual[5000:]), color=color[k])
    poptq, pcovq = curve_fit(hpr.affine, np.log(sminus[5000:]), np.log(residual[5000:]))

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('log(sgminus)'); plt.ylabel('log(sminus)')
plt.tight_layout()
plt.savefig(os.path.join(mother,'log_log.png'), dpi=600); plt.clf()
 