'''
Plot aQ as a function of bQ for different SMs.
'''

import os
import gc
import time
import pickle
import numpy as np
import matplotlib.pyplot as plt
import concurrent.futures
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
import dynalysis.basics as dbcs
import dynalysis.classes as clss
from math import sqrt
from scipy.optimize import curve_fit

'''=====Define parameters here====='''
gmavals = [[0.3,0.6,0.5,1,1], [0.3,0.6,0.7,1,1], [0.3,0.6,0.7,1,1],\
           [0.3,0.6,0.8,1,1], [0.3,0.8,1.25,1,1]] 
'''================================'''

mother = os.getcwd()
nn = []; nn2 = [] # nn: aQ, nn2: bQ
asim, bsim = [], []
astd, bstd = [], []
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
    try:
        popt, pcov = curve_fit(hpr.Quad, sminus, sgminus)
        nn.append(popt[0]); nn2.append(popt[1])
    except:
        nn.append(0); nn2.append(popt[1])
        
    traces = vsn.inquire_trace_npy(1, fname, 9999)
    alist, blist = [], []
    for trace in traces:
        s1, s2, sg1, sg2, _, _ = trace
        splus = (s1-s2)/2; sgplus = (sg1-sg2)/2
        y1 = splus[50:150]*np.sign(np.mean(splus[50:150]))
        y2 = sgplus[50:150]*np.sign(np.mean(splus[50:150]))
        ae, be, _ = curve_fit(hpr.Quad, y1, y2)[0]
        alist.append(ae); blist.append(be)
    asim.append(np.mean(alist)); bsim.append(np.mean(blist))
    astd.append(np.std(alist)); bstd.append(np.std(blist))

fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
ax.plot(np.array(nn), nn2, marker='.', linestyle='', color='r', markersize=2.5, linewidth=1) #quad
# ax.plot(nn2, marker='.', linestyle='', color='hotpink', markersize=2.5, linewidth=1)
ax.errorbar(asim, bsim, xerr=astd, yerr=bstd, color='b', marker='.', markersize=2.5, linestyle='', linewidth=1, label='simulation')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('sg'); plt.ylabel('splus')
plt.tight_layout()
plt.savefig(os.path.join(mother,'quad_coefficients.png'), dpi=600); plt.clf()
 