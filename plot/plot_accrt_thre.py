'''
Plots speed-accuracy trade-off plane for different thresholds.
'''
import os
import pickle
import matplotlib.pyplot as plt
import aadm.basics as bcs
import aadm.helper as hpr
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

'''=====Define parameters here====='''
gmavals = [[0.1,1,1,0,2],[0.3,1,1,0,2],[0.5,1,1,0,2],[0.7,1,1,0,2],\
           [0.75,1,1,0,2],[0.8,1,1,0,2], [0.83,1,1,0,2],
           [0.3,0.6,0.5,1,1], [0.3,0.6,0.6,1,1],\
           [0.5,0.6,0.7,1,1],[0.5,0.6,0.8,1,1],\
           [0.5,0.6,0.7,1.05,1]] 
thres = [0.35,0.4,0.45,0.5,0.55,0.6]

Q = ''
rtname = 'rtc'
'''================================'''

ultimate = {} # key: threshold, value: [rtc, rtcstd, ap, apstd]
for thre in thres: ultimate[thre] = []

for thre in thres:
    print(thre)
    folder = 'files_4e_thre='+str(thre)
    storages = []
    files = [rtname+'all'+Q,rtname+'std'+Q,'ap'+Q,'apstd'+Q]

    for f in range(len(files)):
        infofile = os.path.join(bcs.datapath(), folder, files[f]+'.p')
        storages.append(pickle.load(open(infofile, 'rb')))
    rtcall,rtcstd,ap,apstd = storages
    
    for k in range(len(gmavals)):
        x, xstd, y, ystd = rtcall[tuple(gmavals[k])], rtcstd[tuple(gmavals[k])],\
            ap[tuple(gmavals[k])], apstd[tuple(gmavals[k])]
        ultimate[thre].append([x,xstd,y,ystd])

# depends on how gmavals is setup
N = 7 
M = 9
mks = ['o']*N + ['d','^'] + ['d','^','s']
color = ['mediumaquamarine']*N + ['hotpink']*(M-N) + ['orange']*(len(gmavals)-M)

# plot speed-accuracy trade-off plane for different thresholds
for thre in thres:
    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(111)
    x, xstd, y, ystd = zip(*ultimate[thre])
    for k in range(len(gmavals)):
        plt.errorbar(x[k], y[k], xerr=xstd[k], yerr=ystd[k], color=color[k], linestyle='',\
                     marker=mks[k], markersize=6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.xlabel('rt'); plt.ylabel('acc')
    plt.tight_layout()
    figname = 'accrt_thre='+str(thre)+'.png'
    #plt.savefig(os.path.join(bcs.graphpath(),figname),dpi=600); plt.clf()
    plt.show(); plt.clf()
        
# plot legend
color = ['mediumaquamarine'] + ['hotpink']*(M-N) + ['orange']*(len(gmavals)-M)
mk = ['o'] + ['d','^'] + ['d','^','s']
labels = ['UM']
for k in range(N, len(gmavals)): labels.append('SM, '+str(tuple(gmavals[k])))
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
for g in range(len(mk)):
    plt.plot([0],[0], label=labels[g], color=color[g], marker=mk[g])
plt.legend()
plt.savefig(os.path.join(bcs.graphpath(),'legend_thre.png'), dpi=600); plt.clf()