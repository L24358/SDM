'''
Adds inhibition after decision is made
'''

import os
import time
import concurrent.futures
import numpy as np
import matplotlib.pyplot as plt
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
import dynalysis.classes as clss 
from math import sqrt
from aadm.sode import RungeKutta4

'''=====Parameters====='''
p = 1
c = 0
mu0 = 0.0156/3
Ne = Ni = 50
N = 2*Ne + 2*Ni
snum, enum = 12, 228
gmaval = [0.5,0.6,0.7,1,1]
nsig, tA, dt, nsigg = 5, 0.002, 0.0001, 2.5
crnumdic = {0:1, 3.2:2, 6.4:3, 12.8:4, 25.6:5, 51.2:6}
rnum = crnumdic[c] 

sup = 0.005
thre = 0.17
delay = 0.15
print('suppression: '+str(sup))
print('threshold: '+str(thre))
print('delay: '+str(delay))
print('Ne='+str(Ne)+', p='+str(p)+', nsig='+str(nsig)+', mu0='+str(mu0)+', nsigg='+str(nsigg))

'''=====Network Setup====='''
#Synaptic weights
gma, gma2, gma3, gma4, alpha = gmaval
params = hpr.get_params_s2(gma, gma2, gma3, gma4, alpha)
gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params
gees, giis = np.array([gees, giis])/(Ne*p-1)
geis, gies, gee, gei, gie, gii = np.array([geis, gies, gee, gei, gie, gii])/(Ne*p)

#Connection matrix
P = np.random.binomial(1,p,size=(N,N))
for i in range(N): P[i][i] = 0
np.random.seed(0) ##Temporary seed

if True: #Probabilistic connection
    Pee = np.random.binomial(1,p,size=(2*Ne,2*Ne))
    Pee = np.tril(Pee) + np.triu(Pee.T, 1)
    Pie = np.random.binomial(1,p,size=(2*Ne,2*Ni)) #Only works on Ne=Ni
    Pie = np.tril(Pie) + np.triu(Pie.T, 1)
    Pei = np.random.binomial(1,p,size=(2*Ni,2*Ne))
    Pei = np.tril(Pei) + np.triu(Pei.T, 1)
    Pii = np.random.binomial(1,p,size=(2*Ni,2*Ni))
    Pii = np.tril(Pii) + np.triu(Pii.T, 1)
    P = np.vstack([np.hstack([Pee,Pie]), np.hstack([Pei,Pii])])
    for i in range(N): P[i][i] = 0

if False: #Uniformly distributed weights
    Pee = np.random.uniform(0.99,1.01,size=(2*Ne,2*Ne))
    Pee = np.tril(Pee) + np.triu(Pee.T, 1)
    Pie = np.random.uniform(0.99,1.01,size=(2*Ne,2*Ni)) #Only works on Ne=Ni
    Pie = np.tril(Pie) + np.triu(Pie.T, 1)
    Pei = np.random.uniform(0.99,1.01,size=(2*Ni,2*Ne))
    Pei = np.tril(Pei) + np.triu(Pei.T, 1)
    Pii = np.random.uniform(0.99,1.01,size=(2*Ni,2*Ni))
    Pii = np.tril(Pii) + np.triu(Pii.T, 1)
    P = np.vstack([np.hstack([Pee,Pie]), np.hstack([Pei,Pii])])
    
for i in range(N): P[i][i] = 0

#Folder
Q = 'sup='+str(sup)+'_thre='+str(thre)+'_delay='+str(delay)+'_'
obran = clss.branch('Ne='+str(Ne)+', p='+str(p)+', nsig='+str(nsig)+', mu0='+str(mu0)+', nsigg='+str(nsigg), os.getcwd())
bran = clss.branch(Q+'Ne='+str(Ne)+', p='+str(p)+', nsig='+str(nsig)+', mu0='+str(mu0)+', nsigg='+str(nsigg), os.getcwd())
bran.mkdir()

#Functions
def unpack(v):
    return v[:Ne], v[Ne:2*Ne], v[2*Ne:2*Ne+Ni], v[2*Ne+Ni:N], v[N:N+2*Ne], v[N+2*Ne:]

def cb(sn1, sn2, sg1, sg2, ws, exist):
    a, b, c, d = ws
    return np.dot(list(a*sn1)+list(b*sn2)+list(-c*sg1)+list(-d*sg2), exist)

def downsample(resy):
    return [np.array([train[t*100] for t in range(100)]) for train in resy]

def rhs(v, t):
    sn1, sn2, sg1, sg2, ep, epg = unpack(v)
        
    fn1s = [hpr.S_N(sn1[i], cb(sn1, sn2, sg1, sg2, [gees,gee,gies,gie], P[i])+Ie+ep[i]) for i in range(Ne)]
    fn2s = [hpr.S_N(sn2[i], cb(sn1, sn2, sg1, sg2, [gee,gees,gie,gies], P[i+Ne])+Ie+ep[i+Ne]) for i in range(Ne)]
    fg1s = [hpr.S_G(sg1[i], cb(sn1, sn2, sg1, sg2, [geis,gei,giis,gii], P[i+2*Ne])+Ii+epg[i]+sup) for i in range(Ni)]
    fg2s = [hpr.S_G(sg2[i], cb(sn1, sn2, sg1, sg2, [gei,geis,gii,giis], P[i+2*Ne+Ni])+Ii+epg[i+Ni]+sup) for i in range(Ni)]
    fes = [(-ep[i] + hpr.noise(0,1)*sqrt(tA*nsig**2))/tA for i in range(2*Ne)]
    fegs = [(-epg[i] + hpr.noise(0,1)*sqrt(tA*nsigg**2))/tA for i in range(2*Ni)]
    return np.array(fn1s+fn2s+fg1s+fg2s+fes+fegs)

def read_init(rp, delay=delay):
    resy = np.load(os.path.join(obran.pathlink,'run'+str(1)+'_tr'+str(rp)+'.npy'))
    s1s, s2s, sg1s, sg2s, _, _ = unpack(resy)
    mean_s1, mean_s2 = np.mean(s1s,axis=0), np.mean(s2s,axis=0)
    rt, win = hpr.accuracy(mean_s1, mean_s2, thre=thre); rt = int((rt+delay)*100)+50
    return resy[:,rt]

def solve(rp):
    np.random.seed(rp+10000)
    solver = RungeKutta4(rhs)
    init = read_init(rp)
    solver.set_initial_conditions(init)
    resy, rest = solver.solve(np.arange(0,1,0.0001))
    resy = np.transpose(resy)
    resy = downsample(resy)
    np.save(os.path.join(bran.pathlink,\
                          'run'+str(rnum)+'_tr'+str(rp)+'.npy'), resy)

    for neu in range(Ne): plt.plot(resy[neu], 'g')
    for neu in range(Ne,2*Ne): plt.plot(resy[neu], 'r')
    plt.savefig(os.path.join(bran.pathlink,'sample.png'), dpi=200); plt.clf()
    
    for neu in range(2*Ne,2*Ne+Ni): plt.plot(resy[neu], 'b')
    for neu in range(2*Ne+Ni,2*Ne+2*Ni): plt.plot(resy[neu], color='hotpink')
    plt.savefig(os.path.join(bran.pathlink,'sample_inh.png'), dpi=200); plt.clf()

if __name__=='__main__':
    
    start = time.perf_counter()
    
    print(snum)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(solve, range(snum, enum))

    finish = time.perf_counter()
    print('Total time: '+str(finish-start))
