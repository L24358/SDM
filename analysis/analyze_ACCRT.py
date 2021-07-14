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
           [0.5,1,0.7,1,1],[0.9,1,1,0,2]] #[0.5,0.9,0.9,1.5,1]
'''================================'''

def pick_and_mean(l):
    l = -np.array(l)+2
    k = np.random.choice(l, len(l))
    return np.mean(k)

def pick_and_mean2(l):
    k = np.random.choice(l, len(l))
    return np.mean(k)

def Quad(T, a, b, c):
    return [a*t**2+b*t+c for t in T]
    
def bounds():
    return [(-100,100), (-100,100), (-100,100)]

#define parameters
mother = os.getcwd()
coherence = [0,3.2,6.4,12.8,25.6,51.2]
x = np.linspace(0,np.log(51.2),num=20)[1:]
xexp = [0] + [np.exp(i) for i in x]
x = [0] + list(x)
B = 1000
rtcall, rtcstd, rtcmore = {}, {}, {}
Au, Austd, Auall = {}, {}, {}
ap, apstd, apall = {}, {}, {}

#load info
thre = 0.6
foldername = 'files_4e_thre='+str(thre)
bran = clss.branch(foldername, mother); bran.mkdir()
to_store = [Au,Austd,Auall,rtcall,rtcstd,rtcmore,ap,apstd,apall]
files = ['Au','Austd','Auall','rtcall','rtcstd','rtcmore', 'ap','apstd','apall']
for f in range(len(files)):
    info = pickle.load(open(os.path.join(foldername,files[f]+'.p'), 'rb'))
    to_store[f].update(info)

for k in range(len(gmavals)):
    bootaccs, bootrtcs, bootrtws = [], [], []
    gma, gma2, gma3, gma4, alpha = gmavals[k]
    if tuple(gmavals[k]) in Au.keys():
        print('Have info on '+str(gmavals[k])+'!')
    else:
        fname = 'SMA4e_gma='+str(gma)+', '+str(gma2)+', '+str(gma3)+', '+str(gma4)+', '+str(alpha)
        print('Doing '+fname)
        for rnum in range(1,7):
            traces = vsn.inquire_trace_npy(rnum, fname, 9999)
            accs, rtcs = hpr.psych_func_complete([traces], thre=thre)
            acc, rtc = accs[0], rtcs[0]
            bootacc = [pick_and_mean(acc) for b in range(B)]
            bootrtc = [pick_and_mean2(rtc) for b in range(B)]
            bootaccs.append(bootacc)
            bootrtcs.append(bootrtc)
            
        bootaccs = np.transpose(np.array(bootaccs))
        Aus, Ats = [], []; aps = []
        for b in range(B):
            popt, pcov = curve_fit(hpr.psyched, coherence, bootaccs[b])
            y = [hpr.psyched(c,*popt) for c in xexp]
            Aus.append(sum(y)/len(xexp))
            weight = [hpr.triang(xe) for xe in xexp]
            Ats.append(np.dot(weight, y)/sum(weight))
            aps.append(popt[0])
            
        bootrtcs = np.transpose(np.array(bootrtcs))
        brtcs = []
        for b in range(B):
            def SSE(pm):
                warnings.filterwarnings('ignore')
                val = Quad(coherence, *pm)
                return np.sum((np.array(bootrtcs[b])-np.array(val))**2)
            popt = differential_evolution(SSE, bounds=bounds(), seed=3)['x']
            y = Quad(xexp,*popt)
            brtcs.append(sum(y)/len(xexp))
            
        Au[tuple(gmavals[k])] = np.mean(Aus)
        Austd[tuple(gmavals[k])] = np.std(Aus)
        Auall[tuple(gmavals[k])] = Aus
        rtcall[tuple(gmavals[k])] = np.mean(brtcs)
        rtcstd[tuple(gmavals[k])] = np.std(brtcs)
        rtcmore[tuple(gmavals[k])] = brtcs
        ap[tuple(gmavals[k])] = np.mean(aps)
        apstd[tuple(gmavals[k])] = np.std(aps)
        apall[tuple(gmavals[k])] = aps
    
        #rearrange information & store information
        for f in range(len(files)):
            afile = open(os.path.join('./'+str(foldername),files[f]+'.p'), 'wb') 
            pickle.dump(to_store[f], afile) 
            afile.close()
    