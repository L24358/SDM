'''
Compares theoretical (noiseless) vs averaged (noisy) a, b and aQ, bQ values.
'''
import os
import numpy as np
import matplotlib.pyplot as plt
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
import matplotlib.ticker as ticker
from math import sqrt
from scipy.optimize import fsolve, curve_fit
from scipy.integrate import solve_ivp

mother = os.getcwd()

'''=====UM====='''
if False:
    atheory, btheory = [], []
    asim, bsim = [], []
    astd, bstd = [], []
    gmavals = [[0.1,1,1,0,2], [0.3,1,1,0,2], [0.5,1,1,0,2],\
            [0.75,1,1,0,2], [0.83,1,1,0,2]]

    for k in range(len(gmavals)):
        gma, gma2, gma3, gma4, alpha = gmavals[k]
        fname = 'SMA4e_gma='+str(gma)+', '+str(gma2)+', '+str(gma3)+', '+str(gma4)+', '+str(alpha)
        traces = vsn.inquire_trace_npy(1, fname, 9999)
        params = hpr.get_params_s2(gma, gma2, gma3, gma4, alpha)
        a, b = vsn.get_ab_SM(params); b = b/2 # obtain a, b values for plus direction for noiseless simulation
        atheory.append(a); btheory.append(b)
        
        # obtain averaged a, b values for noisy simulations
        alist, blist = [], []
        for trace in traces:
            s1, s2, sg1, sg2, _, _ = trace; splus = (s1+s2)/2
            splus = (s1+s2)/2; sgplus = (sg1+sg2)/2
            y1 = splus[50:100]
            y2 = sgplus[50:100]
            ae, be = curve_fit(hpr.affine, y1, y2)[0]
            alist.append(ae); blist.append(be)
        asim.append(np.mean(alist)); bsim.append(np.mean(blist))
        astd.append(np.std(alist)); bstd.append(np.std(blist))

    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.errorbar(asim, bsim, xerr=astd, yerr=bstd, color='b', marker='.', markersize=2.5, linestyle='', linewidth=1)
    plt.plot(atheory, btheory, color='r', marker='^', markersize=2.5, linestyle='')
    plt.xlabel('gma'); plt.ylabel('intercept')
    plt.tight_layout()
    plt.savefig(os.path.join(mother,'ab_UM.png'), dpi=600)
'''============'''

'''=====SM, plus====='''
if False:
    atheory, btheory = [], []
    asim, bsim = [], []
    astd, bstd = [], []
    gmavals = [[0.3,0.6,0.5,1,1], [0.3,0.6,0.6,1,1], [0.3,0.6,0.7,1,1],\
            [0.3,0.6,0.8,1,1], [0.3,0.8,1.25,1,1]]
        
    for k in range(len(gmavals)):
        gma, gma2, gma3, gma4, alpha = gmavals[k]
        fname = 'SMA4e_gma='+str(gma)+', '+str(gma2)+', '+str(gma3)+', '+str(gma4)+', '+str(alpha)
        traces = vsn.inquire_trace_npy(1, fname, 9999)
        params = hpr.get_params_s2(gma, gma2, gma3, gma4, alpha)
        a, b = vsn.get_ab_SM(params); b = b/2 # obtain a, b values for plus direction for noiseless simulation
        atheory.append(a); btheory.append(b)
        
        # obtain averaged a, b values for noisy simulations
        alist, blist = [], []
        for trace in traces:
            s1, s2, sg1, sg2, _, _ = trace
            splus = (s1+s2)/2; sgplus = (sg1+sg2)/2
            y1 = splus[50:100]
            y2 = sgplus[50:100]
            ae, be = curve_fit(hpr.affine, y1, y2)[0]
            alist.append(ae); blist.append(be)
        asim.append(np.mean(alist)); bsim.append(np.mean(blist))
        astd.append(np.std(alist)); bstd.append(np.std(blist))
        
    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.errorbar(asim, bsim, xerr=astd, yerr=bstd, color='b', marker='.', markersize=2.5, linestyle='', linewidth=1, label='simulation')
    plt.plot(atheory, btheory, color='r', marker='^', markersize=2.5, linestyle='', label='no noise')
    plt.xlabel('gma'); plt.ylabel('slope')
    plt.tight_layout()
    plt.savefig(os.path.join(mother,'ab_SMplus.png'), dpi=600)
'''=================='''

'''=====SM, minus====='''
if True:
    atheory, btheory = [], []
    asim, bsim = [], []
    astd, bstd = [], []
    gmavals = [[0.3,0.6,0.5,1,1], [0.3,0.6,0.6,1,1], [0.3,0.6,0.7,1,1],\
            [0.3,0.6,0.8,1,1], [0.3,0.8,1.25,1,1]]
        
    for k in range(len(gmavals)):
        gma, gma2, gma3, gma4, alpha = gmavals[k]
        fname = 'SMA4e_gma='+str(gma)+', '+str(gma2)+', '+str(gma3)+', '+str(gma4)+', '+str(alpha)
        traces = vsn.inquire_trace_npy(1, fname, 9999)
        params = hpr.get_params_s2(gma, gma2, gma3, gma4, alpha)
        a, b = vsn.get_abg_SM(params) # obtain aQ, bQ values for plus direction for noiseless simulation
        atheory.append(a); btheory.append(b)
        
        # obtain averaged aQ, bQ values for noisy simulations
        alist, blist = [], []
        for trace in traces:
            s1, s2, sg1, sg2, _, _ = trace
            sminus = (s1-s2)/2; sgminus = (sg1-sg2)/2
            y1 = sminus[50:100]
            y2 = sgminus[50:100]
            ae, be = curve_fit(hpr.affine, y1, y2)[0]
            alist.append(ae); blist.append(be)
        asim.append(np.mean(alist)); bsim.append(np.mean(blist))
        astd.append(np.std(alist)); bstd.append(np.std(blist))
        
    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.errorbar(asim, bsim, xerr=astd, yerr=bstd, color='b', marker='.', markersize=2.5, linestyle='', linewidth=1, label='simulation')
    plt.plot(atheory, btheory, color='r', marker='^', markersize=2.5, linestyle='', label='no noise')
    plt.xlabel('gma'); plt.ylabel('slope')
    plt.tight_layout()
    plt.savefig(os.path.join(mother,'ab_SMminus.png'), dpi=600)
'''==================='''