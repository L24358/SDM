import os
import numpy as np
import matplotlib.pyplot as plt
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
from mpl_toolkits import mplot3d
from scipy.optimize import curve_fit

mother = os.getcwd()
gmavals = [[0.3,0.6,0.6,1,1], [0.3,0.6,0.7,1,1], [0.3,0.6,0.8,1,1], [0.3,0.8,1.25,1,1],\
           [0.5,0.6,0.7,1,1], [0.5,0.6,0.8,1,1], [0.3,0.6,0.5,1,1], [0.5,1,1,0,2]]

def sim_sminus(params, ss=0.25, d=0.01):
    s1, s2, sg1, sg2 = vsn.downsample(vsn.sim_SM_nonoise(ss, params, d = d))
    splus = (s1+s2)/2; sgplus = (sg1+sg2)/2
    sminus = (s1-s2)/2; sgminus = (sg1-sg2)/2
    return splus, sgplus, sminus, sgminus

x, y, c = [], [], []
for k in range(len(gmavals)):
    gma, gma2, gma3, gma4, alpha = gmavals[k]
    params = hpr.get_params_s2(gma, gma2, gma3, gma4, alpha)
    fname = str(gma)+', '+str(gma2)+', '+str(gma3)+', '+str(gma4)+', '+str(alpha)
    
    gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params
    sp_inp = 307.5*geis*(1+gma2)
    sp_rel = 200+307.5*giis*(1+gma4)
    sm_inp = 307.5*geis*(1-gma2)
    sm_rel = 200+307.5*giis*(1-gma4)
    
    splus, sgplus, sminus, sgminus = sim_sminus(params)
    ssplus = (sp_inp*splus-88.5+307.5)/sp_rel
    ssminus = (sm_inp*sminus)/sm_rel
    dp = sgplus - ssplus
    dm = sgminus - ssminus
    
    fig = plt.figure(figsize=(4,3))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax1.plot(np.linspace(0,3,len(dp)), splus, 'b')
    ax1.plot(np.linspace(0,3,len(dm)), sminus, 'r')
    ax2.plot(np.linspace(0,3,len(dp)), dp*100, 'b')
    ax2.plot(np.linspace(0,3,len(dm)), dm*100, 'r')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax1.set_xlabel('t (s)'); ax1.set_ylabel('s+, s-')
    ax2.set_xlabel('t (s)'); ax2.set_ylabel('lag')
    plt.tight_layout()
    plt.savefig(os.path.join(mother,'lag_'+fname+'.png'), dpi=600); plt.clf()
    
    popt, pcov = curve_fit(hpr.Quad, sminus, sgminus)
    
    x.append(min(dm)*100)
    y.append(popt[0])
    c.append(sm_inp)
    
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
plt.plot(c, x, 'k.')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('sminp'); plt.ylabel('lag')
plt.tight_layout()
plt.savefig(os.path.join(mother,'sminp_lag.png'), dpi=600); plt.clf()
    
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
plt.plot(x, y, 'k.')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('lag'); plt.ylabel('aQ')
plt.tight_layout()
plt.savefig(os.path.join(mother,'lag_quad.png'), dpi=600); plt.clf()