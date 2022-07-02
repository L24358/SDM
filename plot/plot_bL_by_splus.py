'''
Calculate crossing time as a function of splus, conditioned on the benchmark value (=0.05).
'''
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
from scipy.optimize import curve_fit

color = ['grey']
gmavals = [[0.5,0.6,0.7,1,1]]
rnum = 1
mother = os.getcwd()

def inn(x, rge):
    '''Determine whether x is in the range given by rge.'''
    if (x >= rge[0]) and (x < rge[1]): return True
    return False

def crossing(seq1, seq2):
    '''Returns crossing time.'''
    ind = None
    for t in range(len(seq1)//5-1):
        popt, pcov = curve_fit(hpr.affine, seq1[t*5:(t+1)*5], seq2[t*5:(t+1)*5])
        if popt[1] <= 0: ind = t; break
    return ind

def sp_by_bm(traces, bmark):
    '''Finds splus value and crossing time for each trial when it hits the benchmark value.'''
    dic = {(0.1,0.125):[], (0.125,0.15):[], (0.15,0.175):[], (0.175,0.5):[]}
    
    for tr in range(len(traces)):
        s1, s2, sg1, sg2 = traces[tr][:4]
        splus = (s1+s2)/2
        sminus = (s1-s2)/2; sgminus = (sg1-sg2)/2
        flag = True
        for t in range(len(sminus)):
            if abs(sminus[t]*2) > bmark: 
                ind = t
                sign = np.sign(np.mean(sminus))
                flag = False; break
        if flag: print('flagged!')
        
        cind = crossing(sminus[ind:]*sign, sgminus[ind:]*sign)
        for key in dic.keys():
            if inn(splus[ind], key) and (cind!=None): dic[key].append(cind*0.05+0.025); break
    
    return dic

# main
xaxis = np.arange(0.1,0.2,0.025)
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
for k in range(len(gmavals)):
    gma, gma2, gma3, gma4, alpha = gmavals[k]
    params = hpr.get_params_s2(gma, gma2, gma3, gma4, alpha)
    fname = 'SMA4e_gma='+str(gma)+', '+str(gma2)+', '+str(gma3)+', '+str(gma4)+', '+str(alpha)
    traces = vsn.inquire_trace_npy(rnum, fname, 1500)
    dic = sp_by_bm(traces, 0.05) # key: splus interval, value: crossing time
    
    ass, astd = [], []
    for key in dic.keys():
        if list(dic[key]) != []:
            ass.append(np.mean(dic[key]))
            astd.append(np.std(dic[key])/np.sqrt(len(dic[key])-1))
        else:
            ass.append(0); astd.append(0)
    plt.errorbar(xaxis, ass, yerr = astd, marker='.', color=color[k], linewidth=0.8)
        
plt.xlabel('splus'); plt.ylabel('bL crossing (s)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(mother,'bLcrossing_sim.png'), dpi=600)
    