'''
Plot averaged 1D energy landscape.
'''

import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import aadm.basics as bcs
import aadm.vision as vsn

def isin(x, low, high):
    if x >= low and x < high: return True
    return False

'''=====Define parameters here====='''
Q = '_bg=0' # _bg=0 or space
mother = os.getcwd()
agvals = list(np.load(os.path.join(bcs.datapath(),'files','agvals.npy')))
agvals = [agvals[x] for x in np.arange(2,38,7)]
time = [np.arange(0,9500,250)[x]/10000 for x in np.arange(2,38,7)]
color = sns.color_palette("hls", 6)
bks = np.arange(-0.9,0.905,0.05)
'''================================'''

'''=====LRM for SM: time-independent landscape====='''
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
fname = 'SMAabv3Q_gma=0.5, 0.6, 0.7, 1, 1'
traces = vsn.inquire_trace_npy(1, fname, 9999)

pts = []
for trace in traces:
    s1, s2 = trace[0], trace[1]
    xs = s1-s2
    chosen = np.array([xs[int(tm*100+50)] for tm in time])
    pts.append(chosen*np.sign(np.mean(chosen)))
pts = np.mean(pts, axis=0)

for pt in pts:
    plt.plot([pt,pt],[0.12,-0.12], color='silver', linestyle='--', linewidth=0.8)
plt.plot([-0.9,0.9], [0,0], color='silver', linestyle='--', linewidth=0.8)


'''=====LRM for UM: time-independent landscape====='''
fname = 'SMAabv3_gma=0.5, 1, 1, 0, 2'
traces = vsn.inquire_trace_npy(1, fname, 9999)
dic = {}

for bk in range(len(bks)-1):
    dic[(bks[bk], bks[bk+1])] = []

for trace in traces:
    s1, s2 = trace[0], trace[1]
    xs = s1-s2
    vs = xs[1:]-xs[:-1]
    for t in range(len(xs)-1):
        for key in dic.keys():
            if isin(xs[t], key[0], key[1]):
                dic[key].append(vs[t]); break

avg_v = []
for key in dic.keys(): 
    if dic[key] != []: avg_v.append(np.mean(dic[key]))
    else: avg_v.append(0)
U = np.array([-sum(avg_v[:i]) for i in range(1,len(avg_v))])
U = [0]+list(U)
xax = np.array([np.mean(key) for key in dic.keys()])
    
start, end = 0, -1
for u in range(len(U)):
    if abs(U[u]) > 0: start = u; break
for u in range(len(U)//2, len(U)):
    if abs(U[u]) < 0.001: end = u; break
end = max([len(U)-start, end])
plt.plot(xax[start:end], U[start:end], color='grey')
'''================================================'''

'''=====LTDRM====='''
for k in range(len(agvals)):
    print(k)
    fname = 'agval='+str(round(agvals[k],5))+Q
    traces = vsn.inquire_trace_npy(1, fname, 9999)
    
    dic = {}
    for bk in range(len(bks)-1):
        dic[(bks[bk], bks[bk+1])] = []
    
    for trace in traces:
        s1, s2 = trace[0], trace[1]
        xs = s1-s2
        vs = xs[1:]-xs[:-1]
        
        for t in range(len(xs)-1):
            flag = True
            for key in dic.keys():
                if isin(xs[t], key[0], key[1]):
                    dic[key].append(vs[t])
                    flag = False
                    break
            if flag: print(xs[t])
       
    avg_v = []
    for key in dic.keys():
        if dic[key] != []: avg_v.append(np.mean(dic[key]))
        else: avg_v.append(0)
    
    U = np.array([-sum(avg_v[:i]) for i in range(1,len(avg_v))])
    U = [0]+list(U)
    xax = np.array([np.mean(key) for key in dic.keys()])
    
    start, end = 0, -1
    for u in range(len(U)):
        if abs(U[u]) > 0: start = u; break
    for u in range(len(U)//2, len(U)):
        if abs(U[u]) < 0.001: end = u; break
    if end != -1:
        end = max([len(U)-start, end])
    plt.plot(xax[start:end], U[start:end], color=color[k])
    
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('x'); plt.ylabel('U')
plt.tight_layout()
plt.savefig(os.path.join(mother,'landcape_1d'+Q+'.png'), dpi=600); plt.clf()

for k in range(len(time)):
    plt.plot([0], [0], label='t = '+str(time[k]), color=color[k])
plt.legend()
plt.savefig(os.path.join(mother,'legend_landcape_1d'+Q+'.png'), dpi=300); plt.clf()
'''==============='''