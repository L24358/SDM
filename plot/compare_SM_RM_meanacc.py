import os
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import aadm.basics as bcs

'''=====Define parameters here====='''
gmavals = [[0.1,1,1,0,2],[0.3,1,1,0,2],[0.5,1,1,0,2],\
           [0.75,1,1,0,2],[0.8,1,1,0,2], [0.83,1,1,0,2],
           [0.3,0.6,0.5,1,1], [0.3,0.6,0.6,1,1], [0.3,0.6,0.7,1,1],\
           [0.3,0.6,0.8,1,1], [0.3,0.8,1.25,1,1]]
'''================================'''

storage = []
files = ['ap_abv3', 'apstd_abv3']
for f in range(len(files)):
    info = pickle.load(open(os.path.join(bcs.datapath(),'files',files[f]+'.p'), 'rb'))
    storage.append(info)
Au, Austd = storage

run = 1
rtc1s, rtc2s = [], []
rtc1std, rtc2std = [], []
for k in range(len(gmavals)):
    rtc1s.append(Au[tuple(gmavals[k])])
    rtc2s.append(Au[tuple(gmavals[k])])
    rtc1std.append(Austd[tuple(gmavals[k])])
    rtc2std.append(Austd[tuple(gmavals[k])])
    
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
minn, maxx = min(rtc1s+rtc2s), max(rtc1s+rtc2s)
plt.plot([minn,maxx],[minn,maxx], linestyle='--', color='gray')
plt.errorbar(rtc1s[:6], rtc2s[:6], color='mediumaquamarine', linestyle='',\
             xerr=rtc1std[:6], yerr=rtc2std[:6])
plt.errorbar(rtc1s[6:], rtc2s[6:], color='hotpink', linestyle='',\
             xerr=rtc1std[6:], yerr=rtc2std[6:])
# plt.errorbar(rtc1s[8:], rtc2s[8:], color='orange', linestyle='',\
#              xerr=rtc1std[8:], yerr=rtc2std[8:])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('full'); plt.ylabel('abv2')
plt.tight_layout()
figname = 'sameseed_compare_meanacc.png'
plt.savefig(os.path.join(bcs.graphpath(),figname),\
            dpi=600); plt.clf()