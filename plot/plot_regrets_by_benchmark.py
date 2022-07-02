'''
Plot the percentage of regret trials as a function of benchmark.
'''

import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn

def regret(traces, bmark, thre):
    nocur, nowin = 0, 0
    no_regrets, corrected, wronged = [], [], []
    for tr in range(len(traces)):
        s1, s2 = traces[tr][0], traces[tr][1]
        sminus = s1 - s2
        flag = True
        for t in range(len(sminus)):
            if abs(sminus[t]) > bmark:
                current = np.sign(sminus[t])
                flag = False; break
        if flag: current = None; nocur += 1
        
        rt, win = hpr.accuracy(s1, s2, thre=0.35)
        if win != None and current != None:
            winner = -2*win + 3
            if winner == current: no_regrets.append(winner)
            elif winner > current: corrected.append(tr)
            else: wronged.append(tr)
        else: nowin += 1
    return no_regrets, corrected, wronged, nocur, nowin

rnum = 1
mother = os.getcwd()
# agvals = list(np.load(os.path.join(bcs.datapath(),'files','agvals.npy'))) # plot regret for LTDRM instead
# agvals = [agvals[x] for x in np.arange(2,38,7)]
agvals = [[0.5,0.6,0.7,1,1], [0.5,1,1,0,2], [0.83,1,1,0,2]]
color = ['grey']+[sns.color_palette("hls", 3)[0], sns.color_palette("hls", 3)[2]]

# main
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
bmarks = list(np.arange(0.0,0.1,0.01)[1:])#+list(np.arange(0.1,0.2,0.05)) # benchmarks
for k in range(len(agvals)):
    # fname = 'agval='+str(round(agvals[k],5)) # plot regret for LTDRM instead
    gma, gma2, gma3, gma4, alpha = agvals[k]
    fname = 'SMA4e_gma='+str(gma)+', '+str(gma2)+', '+str(gma3)+', '+str(gma4)+', '+str(alpha)
    traces = vsn.inquire_trace_npy(rnum, fname, 1500)
    
    cors, wrs, ones = [], [], []
    for bmark in bmarks:
        no_regrets, corrected, wronged, nocur, nowin = regret(traces, bmark, 0.6)
        cors.append(len(corrected))
        wrs.append(len(wronged))
        ones.append(len(no_regrets)/len(traces))

    plt.plot(bmarks, ones, marker='.', color=color[k])
    
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('bmarks'); plt.ylabel('no regrets')
plt.tight_layout()
plt.savefig(os.path.join(mother,'regrets.png'), dpi=600); plt.clf()

lbls = ['SM', 'UM, \u03B3=0.5', 'UM, \u03B3=0.83']
for k in range(len(agvals)):
    plt.plot([0], [0], label=lbls[k], color=color[k])
plt.legend()
plt.savefig(os.path.join(mother,'legend_regrets.png'), dpi=300); plt.clf()