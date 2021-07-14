'''Directly updated noise parameters to fit 4e'''

import os
import gc
import time
import pickle
import concurrent.futures
import numpy as np
import matplotlib.pyplot as plt
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
import dynalysis.basics as dbcs
import dynalysis.classes as clss
from math import sqrt
from datetime import date
from aadm.sode import RungeKutta4

'''=====Define parameters here====='''
Jes = [0.25176]*6
Jis = [0.04381]*6*3
Ibs = [0.32691]*6
snums = bcs.rpt(len(Jes), None)

pmnum = int(len(Jes)/6)
enums = np.repeat([1500,1500,1250,500,300,200], pmnum)
cs = np.repeat([0,3.2,6.4,12.8,25.6,51.2], pmnum)
crnumdic = {0:1, 3.2:2, 6.4:3, 12.8:4, 25.6:5, 51.2:6}
nsig, tA, dt = 0.5*10, 0.002, 0.0001
'''================================'''

def sim(k):
    snum = snums[k]; enum = enums[k]
    c = cs[k]; rnum = crnumdic[c]
    Je, Ji, Ib = Jes[k], Jis[k], Ibs[k]
    fname = '2MA_pm='+str(Je)+', '+str(Ji)+', '+str(Ib)
    
    def rhs(v, t):
        s1, s2, n1, n2 = v
        f1 = hpr.S_N(s1, Je*s1-Ji*s2+Ib+hpr.pulse(hpr.Inp1(c), t-0.5)+n1)
        f2 = hpr.S_N(s2, Je*s2-Ji*s1+Ib+hpr.pulse(hpr.Inp2(c), t-0.5)+n2)
        f3 = (-n1 + hpr.noise(0,1)*sqrt(tA*nsig**2))/tA
        f4 = (-n2 + hpr.noise(0,1)*sqrt(tA*nsig**2))/tA
        return np.array([f1,f2,f3,f4])
    
    bran = clss.branch(fname, bcs.datapath()); bran.mkdir()
    nbran = clss.branch(fname, bcs.datapath2()); nbran.mkdir()
    if snum==None:
        snum = bcs.get_max(bran.pathlink, rnum)+1
    print(snum)
    
    pdic = {'Je':Je, 'Ji':Ji, 'Ib':Ib, 'c':c, 'nsig':nsig, 'dt':dt}
    to_save_list = [key+': '+str(pdic[key]) for key in pdic]
    today = date.today()
    dbcs.output_line(os.path.join(bran.pathlink, 'run'+str(rnum)+'_param.dat'), '\n'+today.strftime("%d/%m/%Y"))
    dbcs.output_single_col(os.path.join(bran.pathlink, 'run'+str(rnum)+'_param.dat'), to_save_list)
    
    check = True
    for rp in range(snum, enum):
        np.random.seed(rp)
        solver = RungeKutta4(rhs)
        solver.set_initial_conditions(list(hpr.get_EQ_u1(hpr.get_params_u1(0.3))[:-1])+[0,0])
        resy, rest = solver.solve(np.arange(0,3,dt))
        
        if check:
            for i in range(2): plt.plot(rest, resy[:,i])
            plt.savefig(fname+'.png'); plt.clf()
            check = False
        
        ntrace = np.transpose(np.asarray(resy))
        dtrace = vsn.downsample(ntrace)
        nfile = os.path.join(nbran.pathlink,'run'+str(rnum)+'_tr'+str(rp)+'.npy')
        dfile = os.path.join(bran.pathlink,'run'+str(rnum)+'_tr'+str(rp)+'.npy')
        np.save(nfile, ntrace)
        np.save(dfile, dtrace)
        gc.collect()
        gc.collect()
        
    print('rnum'+str(rnum)+' done!')

if __name__=='__main__':
    start = time.perf_counter()
    
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(sim, range(len(Jes)))

    finish = time.perf_counter()
    print('Total time: '+str(finish-start))