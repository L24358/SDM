'''
Parametrize accuracy and reaction time.
'''

import os
import pickle
import warnings
import numpy as np
import matplotlib.pyplot as plt
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
import dynalysis.classes as clss
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution

'''=====Define parameters here====='''
gmavals = [[0.1,1,1,0,2],[0.3,1,1,0,2],[0.5,1,1,0,2],[0.7,1,1,0,2],\
           [0.75,1,1,0,2],[0.8,1,1,0,2], [0.83,1,1,0,2],
           [0.1,0.6,0.4,1,1],[0.1,0.6,0.45,1,1], [0.1,0.6,0.5,1,1],[0.1,0.6,0.6,1,1],\
           [0.1,0.6,0.7,1,1],[0.1,0.6,0.8,1,1],\
           [0.3,0.6,0.5,1,1], [0.3,0.6,0.6,1,1], [0.3,0.6,0.7,1,1],\
           [0.3,0.6,0.8,1,1], [0.3,0.8,1.25,1,1],\
           [0.5,0.6,0.7,1,1],[0.5,0.6,0.8,1,1], [0.75,0.9,0.9,1,1],\
           [0.5,0.6,0.7,1.05,1],[0.5,0.6,0.7,0.5,1],[0.5,0.5,0.8,1,1],\
           [0.5,0.9,0.9,1.2,1], [0.5,0.6,0.7,0,2], [0.5,0.6,0.6,0,2],\
           [0.5,0.7,0.7,1,1],[0.5,0.8,0.7,1,1],[0.5,0.9,0.7,1,1],\
           [0.5,1,0.7,1,1],[0.9,1,1,0,2]]

B = 1000 # number of bootstrap samples
thre = 0.6 # threshold value

coherence = [0,3.2,6.4,12.8,25.6,51.2] # coherence values used in simulations
x = np.linspace(0,np.log(51.2),num=20) # coherence values used in plotting
xexp = [0] + [np.exp(i) for i in x[1:]]
'''================================'''

def pick_and_mean(l):
    '''Sample from list (for accurracy).'''
    l = -np.array(l)+2
    k = np.random.choice(l, len(l))
    return np.mean(k)

def pick_and_mean2(l):
    '''Sample from list (for reaction time).'''
    k = np.random.choice(l, len(l))
    return np.mean(k)

def quad(T, a, b, c):
    return [a*t**2+b*t+c for t in T]
    
def bounds():
    return [(-100,100), (-100,100), (-100,100)]

# define path and folder
mother = os.getcwd() # parent path
foldername = 'files_4e_thre='+str(thre)
bran = clss.branch(foldername, mother); bran.mkdir()

# load pre-existing information
rtcall, rtcstd, rtcmore = {}, {}, {} # rtc: defined below, all: mean values, std: std values, more: the whole sequence
Auall, Austd, Aumore = {}, {}, {} # Au: defined below
apall, apstd, apmore = {}, {}, {} # ap: defined below
to_store = [Auall,Austd,Aumore,rtcall,rtcstd,rtcmore,apall,apstd,apmore]
files = ['Auall','Austd','Aumore','rtcall','rtcstd','rtcmore', 'apall','apstd','apmore']
for f in range(len(files)):
    info = pickle.load(open(os.path.join(foldername,files[f]+'.p'), 'rb')) # read pre-existing files
    to_store[f].update(info) # update dictionaries

# main
for k in range(len(gmavals)):

    bootaccs, bootrtcs, bootrtws = [], [], []
    gma, gma2, gma3, gma4, alpha = gmavals[k]

    if tuple(gmavals[k]) in Auall.keys(): # skip over processed files
        print('Have info on '+str(gmavals[k])+'!')

    else:
        fname = 'SMA4e_gma='+str(gma)+', '+str(gma2)+', '+str(gma3)+', '+str(gma4)+', '+str(alpha)
        print('Fitting parameters for '+fname)

        # load data
        for rnum in range(1,7):
            traces = vsn.inquire_trace_npy(rnum, fname, 9999)
            accs, rtcs = hpr.psych_func_complete([traces], thre=thre)
            acc, rtc = accs[0], rtcs[0]
            bootacc = [pick_and_mean(acc) for b in range(B)]
            bootrtc = [pick_and_mean2(rtc) for b in range(B)]
            bootaccs.append(bootacc)
            bootrtcs.append(bootrtc)
            
        # fit parameters for accuracy
        bootaccs = np.transpose(np.array(bootaccs))
        Aus = []; aps = []
        for b in range(B):
            popt, pcov = curve_fit(hpr.psyched, coherence, bootaccs[b])
            y = [hpr.psyched(c,*popt) for c in xexp]
            Aus.append(sum(y)/len(xexp)) # DEFINE Au: the area under curve for accuracy v. coherernce
            aps.append(popt[0]) # DEFINE ap: the fitted parameter "a" for the accuracy v. coherence curve
            
        # fit parameters for reaction time
        bootrtcs = np.transpose(np.array(bootrtcs))
        brtcs = []
        for b in range(B):
            def SSE(pm):
                warnings.filterwarnings('ignore')
                val = quad(coherence, *pm)
                return np.sum((np.array(bootrtcs[b])-np.array(val))**2)
            popt = differential_evolution(SSE, bounds=bounds(), seed=3)['x']
            y = quad(xexp,*popt)
            brtcs.append(sum(y)/len(xexp)) # DEFINE rtc: the area under curve for reaction time v. coherence
            
        # append results to dictionary
        Auall[tuple(gmavals[k])] = np.mean(Aus)
        Austd[tuple(gmavals[k])] = np.std(Aus)
        Aumore[tuple(gmavals[k])] = Aus
        apall[tuple(gmavals[k])] = np.mean(aps)
        apstd[tuple(gmavals[k])] = np.std(aps)
        apmore[tuple(gmavals[k])] = aps
        rtcall[tuple(gmavals[k])] = np.mean(brtcs)
        rtcstd[tuple(gmavals[k])] = np.std(brtcs)
        rtcmore[tuple(gmavals[k])] = brtcs
    
        # rearrange & store information
        for f in range(len(files)):
            afile = open(os.path.join('./'+str(foldername),files[f]+'.p'), 'wb') 
            pickle.dump(to_store[f], afile) 
            afile.close()
    