import os
import pickle
import numpy as np
import aadm.basics as bcs
import aadm.helper as hpr
import dynalysis.classes as clss
from scipy.optimize import curve_fit
from aadm.sode import RungeKutta4

'''=====basics====='''

def downsample(resy):
    return [np.array([train[t*100] for t in range(300)]) for train in resy]

'''=====Simulation====='''

def sim_SM_nonoise(ss, params,d=0, c=0):
    gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params
    def rhs(v, t):
        s1, s2, sG1, sG2 = v
        f1 = hpr.S_N(s1, gees*s1+gee*s2-gies*sG1-gie*sG2+Ie+hpr.pulse(hpr.Inp1(c), t))
        f2 = hpr.S_N(s2, gees*s2+gee*s1-gies*sG2-gie*sG1+Ie+hpr.pulse(hpr.Inp2(c), t))
        f3 = hpr.S_G(sG1, geis*s1+gei*s2-giis*sG1-gii*sG2+Ii)
        f4 = hpr.S_G(sG2, geis*s2+gei*s1-giis*sG2-gii*sG1+Ii)
        return np.array([f1,f2,f3,f4])

    init = [ss+d, ss] + list(hpr.get_ssSG_s1(ss+d, ss, params))
    solver = RungeKutta4(rhs)
    solver.set_initial_conditions(init)
    resy, rest = solver.solve(np.arange(0,3,0.0001)) 
    return np.transpose(np.asarray(resy))
    
def sim_2M_nonoise(ss, params, a, b, ag, d=0, c=0):
    gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params
    m, n = gies+gie, gies-gie
    def rhs(v, t):
        s1, s2 = v
        f1 = hpr.S_N(s1, gees*s1+gee*s2-m*(a*(s1+s2)+b)/2-n*(ag*(s1-s2)+0)/2+Ie+hpr.pulse(hpr.Inp1(c), t))
        f2 = hpr.S_N(s2, gees*s2+gee*s1-m*(a*(s1+s2)+b)/2+n*(ag*(s1-s2)+0)/2+Ie+hpr.pulse(hpr.Inp2(c), t))
        return np.array([f1,f2])

    init = [ss+d, ss]
    solver = RungeKutta4(rhs)
    solver.set_initial_conditions(init)
    resy, rest = solver.solve(np.arange(0,3,0.0001)) 
    return np.transpose(np.asarray(resy))
    
def sim_2MQ_nonoise(ss, params, a, b, ag, bg, d=0, c=0):
    gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params
    m, n = gies+gie, gies-gie
    def rhs(v, t):
        s1, s2 = v
        f1 = hpr.S_N(s1, gees*s1+gee*s2-m*(a*(s1+s2)+b)/2-n*(ag*((s1-s2)/2)**2+bg*((s1-s2)/2)+0)+Ie+hpr.pulse(hpr.Inp1(c), t))
        f2 = hpr.S_N(s2, gees*s2+gee*s1-m*(a*(s1+s2)+b)/2+n*(ag*((s1-s2)/2)**2+bg*((s1-s2)/2)+0)+Ie+hpr.pulse(hpr.Inp2(c), t))
        return np.array([f1,f2])

    init = [ss+d, ss]
    solver = RungeKutta4(rhs)
    solver.set_initial_conditions(init)
    resy, rest = solver.solve(np.arange(0,3,0.0001)) 
    return np.transpose(np.asarray(resy))
    
def get_ab_SM(params):
    ss = hpr.get_EQ_s1(params)[0]
    s1, s2, sg1, sg2 = sim_SM_nonoise(ss, params)
    splus = (s1+s2); sgplus = (sg1+sg2)
    popt, pcov = curve_fit(hpr.affine, splus, sgplus)
    return popt
    
def get_abg_SM(params, ss=0.25, d=0.01):
    s1, s2, sg1, sg2 = downsample(sim_SM_nonoise(ss, params, d = d))
    sminus = (s1-s2)/2; sgminus = (sg1-sg2)/2
    rt = int(hpr.accuracy(s1, s2, t0=0)[0]*10000)
    popt, pcov = curve_fit(hpr.affine, sminus[:rt], sgminus[:rt])
    return popt
    
def get_abg_SM2(params, ss=0.25, d=0.01):
    s1, s2, sg1, sg2 = downsample(sim_SM_nonoise(ss, params, d = d))
    sminus = (s1-s2)/2; sgminus = (sg1-sg2)/2
    rt = int(hpr.accuracy(s1, s2, t0=0)[0]*10000)
    popt, pcov = curve_fit(hpr.Quad, sminus[:rt], sgminus[:rt])
    return popt

'''=====Fetch data====='''

def inquire_trace_sode(run, fname, rnum, seed, mother=bcs.datapath2(), snum=0, rdm=False, down=True):
    if seed !=None: np.random.seed(seed)
    b = clss.branch(fname, mother); os.chdir(b.pathlink)
    #print(b.pathlink)
    enum = bcs.get_max(b.pathlink, run)
    if enum < rnum: rnum = enum+1
    if rdm: rps = np.random.choice(range(snum, enum), rnum)
    else: rps = range(snum, rnum)
    traces = []
    print(rps) ##
    for rp in rps:
        sim = pickle.load(open('run'+str(run)+'_tr'+str(rp)+'.p', 'rb'))
        resy = np.transpose(np.asarray(sim['res.y']))
        if down: traces.append(downsample(resy))
        else: traces.append(resy)
    os.chdir(mother)
    return traces
    
def inquire_trace_npy(run, fname, rnum, seed=None, mother=bcs.datapath(), snum=0, rdm=False):
    if seed !=None: np.random.seed(seed)
    b = clss.branch(fname, mother); os.chdir(b.pathlink)
    enum = bcs.get_max(b.pathlink, run)
    if enum < rnum: rnum = enum+1
    if rdm: rps = np.random.choice(range(snum, enum), rnum)
    else: rps = range(snum, rnum)
    traces = []
    for rp in rps:
        resy = np.load('run'+str(run)+'_tr'+str(rp)+'.npy')
        traces.append(resy)
    os.chdir(mother)
    return traces
    
def rhs_Quad(t,v,params, a, b, ag, bg,c=0, cg=0):
    gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params 
    m, n = gies+gie, gies-gie
    s1, s2 = v
    f1 = hpr.S_N(s1, gees*s1+gee*s2-m*(a*(s1+s2)+b)/2-n*(ag*((s1-s2)/2)**2+bg*((s1-s2)/2)+cg)+Ie+hpr.pulse(hpr.Inp1(c), t-0.5))
    f2 = hpr.S_N(s2, gees*s2+gee*s1-m*(a*(s1+s2)+b)/2+n*(ag*((s1-s2)/2)**2+bg*((s1-s2)/2)+cg)+Ie+hpr.pulse(hpr.Inp2(c), t-0.5))
    return [f1,f2]
    
def rhs_lin(t,v,params, a, b, ag, c=0, bg=0):
    gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params 
    m, n = gies+gie, gies-gie
    s1, s2 = v
    f1 = hpr.S_N(s1, gees*s1+gee*s2-m*(a*(s1+s2)+b)/2-n*(ag*(s1-s2)/2+bg)+Ie+hpr.pulse(hpr.Inp1(c), t-0.5))
    f2 = hpr.S_N(s2, gees*s2+gee*s1-m*(a*(s1+s2)+b)/2+n*(ag*(s1-s2)/2+bg)+Ie+hpr.pulse(hpr.Inp2(c), t-0.5))
    return [f1,f2]
