'''
Simulate time-dependent linear reduced models in parallel.

* Notes:
- Here, since bg!=0 (i.e. it is not a fair trial), the stimulus input begins at t=0.
    (In other simulation trials, stimulus input starts at t=50 ms for ease of verifying stability of original attractor.)
'''

import os
import gc
import time
import pickle
import numpy as np
import matplotlib.pyplot as plt
import concurrent.futures
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
import dynalysis.basics as dbcs
import dynalysis.classes as clss
from math import sqrt
from datetime import date
from aadm.sode import RungeKutta4

'''=====Define parameters here====='''
agvals = list(np.load(os.path.join(bcs.datapath(),'files','agvals.npy'))) # a_L values to loop over
agvals = [agvals[x] for x in np.arange(2,38,7)]*6
bgvals = list(np.load(os.path.join(bcs.datapath(),'files','bgvals.npy'))) # b_L values to loop over
bgvals = [bgvals[x] for x in np.arange(2,38,7)]*6
gmaval = [0.5,0.6,0.7,1,1] # parameters for the full model in which a_L, b_L were extracted from
snums = bcs.rpt(len(agvals), None)

pmnum = int(len(agvals)/6)
enums = np.repeat([1500,1500,1500,1000,300,200], pmnum) # number of trials to simulate for each coherence value
cs = np.repeat([0,3.2,6.4,12.8,25.6,51.2], pmnum) # coherence values
crnumdic = {0:1, 3.2:2, 6.4:3, 12.8:4, 25.6:5, 51.2:6}
nsig, tA, dt = 0.5*10, 0.002, 0.0001 # noise and time steps
'''================================'''

def sim_seed(k):

    # extract simulation parameters for this trial
    snum = snums[k]; enum = enums[k]
    c = cs[k]; rnum = crnumdic[c]
    gma, gma2, gma3, gma4, alpha = gmaval
    params = hpr.get_params_s2(gma, gma2, gma3, gma4, alpha)
    gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params
    a, b = vsn.get_ab_SM(params)
    ag, bg = agvals[k], bgvals[k]
    
    # ODE for time-dependent linear reduced model
    m, n = gies+gie, gies-gie
    def rhs_abv2(v, t):
        s1, s2, n1, n2 = v
        f1 = hpr.S_N(s1, gees*s1+gee*s2-m*(a*(s1+s2)+b)/2-n*(ag*(s1-s2)/2+bg)+Ie+hpr.pulse(hpr.Inp1(c), t)+n1)
        f2 = hpr.S_N(s2, gees*s2+gee*s1-m*(a*(s1+s2)+b)/2+n*(ag*(s1-s2)/2+bg)+Ie+hpr.pulse(hpr.Inp2(c), t)+n2)
        f4 = (-n1 + hpr.noise(0,1)*sqrt(tA*nsig**2))/tA
        f5 = (-n2 + hpr.noise(0,1)*sqrt(tA*nsig**2))/tA
        return np.array([f1,f2,f4,f5])
    
    # define path, folders, starting index
    func = rhs_abv2; init = list(hpr.get_EQ_ab_s2(params, a, b))+[0,0]
    fname = 'agval='+str(round(agvals[k],5))
    bran = clss.branch(fname, bcs.datapath()); bran.mkdir()
    if snum==None:
        snum = bcs.get_max(bran.pathlink, rnum)+1
    print(snum)
    
    # save parameters
    pdic = {'gees':gees, 'gies':gies, 'geis':geis, 'giis':giis,\
            'gee':gee, 'gie':gie, 'gei':gei, 'gii':gii,\
            'Ii':Ii, 'Ie':Ie, 'gma':gma, 'gma2':gma2, 'gma3':gma3,\
            'c':c, 'nsig':nsig, 'dt':dt}
    to_save_list = [key+': '+str(pdic[key]) for key in pdic]
    today = date.today()
    dbcs.output_line(os.path.join(bran.pathlink, 'run'+str(rnum)+'_param.dat'), '\n'+today.strftime("%d/%m/%Y"))
    dbcs.output_single_col(os.path.join(bran.pathlink, 'run'+str(rnum)+'_param.dat'), to_save_list)
    
    # main
    check = True # check simulation result (for one trial)
    for rp in range(snum, enum):
        np.random.seed(rp)
        solver = RungeKutta4(func)
        solver.set_initial_conditions(init)
        resy, rest = solver.solve(np.arange(0,3,dt)) 
        
        if check:
            for i in range(2): plt.plot(rest, resy[:,i])
            plt.savefig(fname+'.png'); plt.clf()
            check = False
        
        ntrace = np.transpose(np.asarray(resy))
        dtrace = vsn.downsample(ntrace)
        dfile = os.path.join(bran.pathlink,'run'+str(rnum)+'_tr'+str(rp)+'.npy')
        np.save(dfile, dtrace)
        gc.collect()
        
    print('rnum'+str(rnum)+' done!')

if __name__=='__main__':
    start = time.perf_counter()
    
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(sim_seed, range(len(agvals)))

    finish = time.perf_counter()
    print('Total time: '+str(finish-start))
    
    
    
    