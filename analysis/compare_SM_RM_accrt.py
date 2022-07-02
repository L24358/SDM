'''
Get accuracy and reaction time information for both full and reduced models for comparison purposes.

* Notes:
- abv3: matches the noise of 4e.
'''

import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn

'''=====Define parameters here====='''
gmavals = [[0.1,1,1,0,2],[0.3,1,1,0,2],[0.5,1,1,0,2],\
           [0.75,1,1,0,2],[0.8,1,1,0,2], [0.83,1,1,0,2],
           [0.3,0.6,0.5,1,1], [0.3,0.6,0.6,1,1], [0.3,0.6,0.7,1,1],\
           [0.3,0.6,0.8,1,1], [0.3,0.8,1.25,1,1], [0.1,0.6,0.5,1,1],
           [0.5,0.6,0.7,1,1], [0.5,0.6,0.7,1.05,1]]
enums = [1000,1000,1000,500,200,200] # number of trials to take into consideration for each coherence
Q = 'Q' # identifier for reduced model type
'''================================'''

# load pre-existing information
accdic1, rtcdic1 = {}, {}
accdic2, rtcdic2 = {}, {}
to_store = [accdic1, accdic2, rtcdic1, rtcdic2]
files = ['fullacc_4e'+Q, 'fullacc_abv3'+Q, 'fullrtc_4e'+Q, 'fullrtc_abv3'+Q]
for f in range(len(files)):
    info = pickle.load(open(os.path.join('files',files[f]+'.p'), 'rb'))
    to_store[f].update(info)

# main
for k in range(len(gmavals)):
    gma, gma2, gma3, gma4, alpha = gmavals[k]
    fname2 = 'SMAabv3'+Q+'_gma='+str(gma)+', '+str(gma2)+', '+str(gma3)+', '+str(gma4)+', '+str(alpha)
    fname1 = 'SMA4e_gma='+str(gma)+', '+str(gma2)+', '+str(gma3)+', '+str(gma4)+', '+str(alpha)
    
    # skip over processed files
    do1, do2 = True, True 
    if tuple(gmavals[k]) in accdic1.keys(): do1 = False
    if tuple(gmavals[k]) in accdic2.keys() or not\
        os.path.exists(os.path.join(bcs.datapath(), fname2)): do2 = False

    atraces1, atraces2 = [], []
    if do1 or do2: print('Extracting accuracy and reaction time information for '+fname1)

    # extract information for full model
    if do1:
        for rnum in range(1,7):
            traces1 = vsn.inquire_trace_npy(rnum, fname1, enums[rnum-1])
            atraces1.append(traces1)
        accs1, rtcs1 = hpr.psych_func_complete(atraces1, thre=0.35)
        accdic1[tuple(gmavals[k])] = accs1
        rtcdic1[tuple(gmavals[k])] = rtcs1

    # extract information for reduced model
    if do2:
        for rnum in range(1,7):
            traces2 = vsn.inquire_trace_npy(rnum, fname2, enums[rnum-1])
            atraces2.append(traces2)
        accs2, rtcs2 = hpr.psych_func_complete(atraces2, thre=0.35)
        if accs2[0]!=[]: accdic2[tuple(gmavals[k])] = accs2
        if rtcs2[0]!=[]: rtcdic2[tuple(gmavals[k])] = rtcs2

    # rearrange & store information
    for f in range(len(files)):
        afile = open(os.path.join('./files',files[f]+'.p'), 'wb') 
        pickle.dump(to_store[f], afile) 
        afile.close()
    
        
        