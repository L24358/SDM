'''
Plot traces of simulations in gif.
'''

import os
import imageio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
import dynalysis.classes as clss

mother =os.getcwd()
b = clss.branch('gif_trace', mother); b.mkdir()

fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)

traces1 = vsn.inquire_trace_npy(1, 'SMA4e_gma=0.5, 0.6, 0.7, 1, 1', 10, mother=bcs.datapath())
traces2 = vsn.inquire_trace_npy(1, 'SMA4e_gma=0.5, 1, 1, 0, 2', 10, mother=bcs.datapath())

T = 300
for t in range(0,T,5):
    fig, ax = plt.subplots(figsize=(8,8))
    for trace in traces1:
        s1, s2, _, _, _, _ = trace
        plt.plot(s1[:t], s2[:t], 'g', linewidth=1)
    for trace in traces2:
        s1, s2, _, _, _, _ = trace
        plt.plot(s1[:t], s2[:t], 'r', linewidth=1)
       
    if False:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.xlabel('$S_1$')
        plt.ylabel('$S_2$')
    else: #turn off all axes
        plt.axis('off')
    plt.xlim((0,0.7)); plt.ylim((0,0.7))
    plt.tight_layout()
    plt.savefig(os.path.join(b.pathlink,'t={:.2f}.png'.format(t)), dpi=200)
    plt.close("all")
    
images = []       
for t in range(0,T,5):
    images.append(imageio.imread(os.path.join(b.pathlink, 't={:.2f}.png'.format(t))))
imageio.mimsave(os.path.join(b.pathlink,'trace.gif'), images, duration=0.2)
    
                