'''
Plot the trial-by-trial comparison of the full model and reduced model.
'''

import os
import pickle
import matplotlib.pyplot as plt
import aadm.basics as bcs

'''=====Define parameters here====='''
gmavals = [[0.3,0.6,0.5,1,1], [0.3,0.6,0.7,1,1], [0.3,0.6,0.7,1,1],\
           [0.3,0.6,0.8,1,1], [0.3,0.8,1.25,1,1], [0.3,1,1,0,2]]
'''================================'''

Q = ''
storage = []
files = ['fullacc_4e'+Q, 'fullacc_abv3'+Q, 'fullrtc_4e'+Q, 'fullrtc_abv3'+Q]
for f in range(len(files)):
    info = pickle.load(open(os.path.join(bcs.datapath(),'files',files[f]+'.p'), 'rb'))
    storage.append(info)
accdic1, accdic2, rtcdic1, rtcdic2 = storage

run = 1
for k in range(len(gmavals)):
    rtc1 = list(rtcdic1[tuple(gmavals[k])][run-1])
    rtc2 = list(rtcdic2[tuple(gmavals[k])][run-1])
    
    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(111)
    plt.plot(rtc1, rtc2, color='k', marker='.', linestyle='',\
             markersize=3)
    minn, maxx = min(rtc1+rtc2), max(rtc1+rtc2)
    plt.plot([minn,maxx],[minn,maxx], linestyle='--', color='gray')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.xlabel('full'); plt.ylabel('abv2')
    plt.tight_layout()
    to_str = [str(gma) for gma in gmavals[k]]
    figname = 'sameseed_'+'_'.join(to_str)+Q+'.png'
    plt.savefig(os.path.join(bcs.graphpath(),figname),\
                dpi=600); plt.clf()