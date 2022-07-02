import os
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import aadm.basics as bcs

'''=====Define parameters here====='''
gmavals = [[0.3,0.6,0.5,1,1], [0.3,0.6,0.6,1,1], [0.3,0.6,0.7,1,1],\
           [0.3,0.6,0.8,1,1]]
'''================================'''

foldername = 'files_4e_thre='+str(0.35)
storage = []
files = ['fullacc_4eQ', 'fullacc_abv3Q', 'fullrtc_4eQ', 'fullrtc_abv3Q']
for f in range(len(files)):
    info = pickle.load(open(os.path.join(bcs.datapath(),foldername,files[f]+'.p'), 'rb'))
    storage.append(info)
accdic1, accdic2, rtcdic1, rtcdic2 = storage

run = 1
rtc1s, rtc2s = [], []
rtc1std, rtc2std = [], []
for k in range(len(gmavals)):
    rtc1 = list(rtcdic1[tuple(gmavals[k])][run-1])
    rtc2 = list(rtcdic2[tuple(gmavals[k])][run-1])
    rtc1s.append(np.mean(rtc1)); rtc2s.append(np.mean(rtc2))
    rtc1std.append(np.std(rtc1)/np.sqrt(len(rtc1)-1)); rtc2std.append(np.std(rtc2)/np.sqrt(len(rtc2)-1))
    
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
minn, maxx = min(rtc1s+rtc2s), max(rtc1s+rtc2s)
plt.plot([minn,maxx],[minn,maxx], linestyle='--', color='gray')
plt.errorbar(rtc1s, rtc2s, color='hotpink', linestyle='',\
             xerr=rtc1std, yerr=rtc2std)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('full'); plt.ylabel('abv2')
plt.tight_layout()
figname = 'sameseed_compare_meanrtcQ.png'
plt.savefig(os.path.join(bcs.graphpath(),figname),\
            dpi=600); plt.clf()