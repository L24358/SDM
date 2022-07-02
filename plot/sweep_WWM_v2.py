'''
Plot parameter sweep of WWM as 2D slices.
'''

import os
import gc
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import aadm.basics as bcs
import aadm.helper as hpr
import dynalysis.classes as clss
from scipy.optimize import fsolve
from itertools import product

inp1 = inp2 = 0

def same(lst, x):
    '''Determines if x belongs to one of the equilibrium points already found in lst.'''
    for item in lst:
        deter = np.any([round(x[i],3)==round(item[i],3) for i in range(2)])
        if deter: return True
    return False

def deter_root_count(params):
    '''Determines the number of roots.'''
    J1, J2, Ie = params
    def nullcline_UM(z):
        s1 = z[0]; s2 = z[1]
        F = np.empty((2))
        F[0] = -s1/0.1 + (1-s1)*0.641*bcs.H(J1*s1-J2*s2+Ie+inp1)
        F[1] = -s2/0.1 + (1-s2)*0.641*bcs.H(J1*s2-J2*s1+Ie+inp2)
        return F
    
    all_res = []
    for s10 in np.arange(0,1,0.01):
            res = fsolve(nullcline_UM, [s10,1-s10])
            if not same(all_res, res): all_res.append(res)
            
    all_res = np.asarray(all_res)
    return all_res[:,0], all_res[:,1]

def get_coors(x, y):
    dic = {}
    for i in range(len(x)): dic[x[i]] = y[i]
    nx = sorted(x); ny = [dic[j] for j in nx]
    return list(zip(nx, ny))

def check_distance(coors):
    vel = hpr.get_velocity(coors)
    mag = [np.linalg.norm(v) for v in vel]
    if min(mag) < 0.2: return False
    return True

def check_region(coors):
    for c in coors:
        if sum(c) > 0.7: return False
    return True
    
b_pics = clss.branch('solns', os.getcwd()); b_pics.mkdir()
pm_J1 = np.arange(0.26,0.27,0.001)[1:]
pm_J2 = np.arange(0.01,0.26,0.005)[1:]
pm_Ie = np.arange(0.31,0.33,0.001)[8:-2]

plot = False; check = False
note = []
pmlist = list(product(pm_J1, pm_J2))
for Ie in pm_Ie[4:8]:
    matrix = []
    for J1 in pm_J1:
        temp = []
        for J2 in pm_J2:
            x, y = deter_root_count([J1, J2, Ie])
            coors = get_coors(x, y)
            criterion = check_distance(coors)
            criterion2 = check_region(coors)
            if len(x) == 5 and criterion and criterion2:
                temp.append(1)
                note.append([J1,J2,Ie])
                if plot:
                    plt.plot(x, y, 'k.')
                    fname = 'J1='+str(J1)+'_J2='+str(J2)+'_Ie='+str(Ie)+'.png'
                    plt.savefig(os.path.join(b_pics.pathlink, fname)); plt.clf()
            else:
                temp.append(0)
                if check and len(x)!=5:
                    plt.plot(x, y, 'k.')
                    fname = 'J1='+str(J1)+'_J2='+str(J2)+'_Ie='+str(Ie)+'.png'
                    plt.savefig(os.path.join(b_pics.pathlink, fname)); plt.clf()
        matrix.append(temp)
    plt.close('all')
    ax = sns.heatmap(matrix, cmap=['whitesmoke','orange'])
    ax.set_xticks(np.arange(0,50,10))
    ax.set_yticks(np.arange(0,10,2))
    # ax.tick_params(left=False, bottom=False)
    # ax.set(xticklabels=[]); ax.set(yticklabels=[])
    plt.savefig('Ie='+str(Ie)+'.png', dpi=300); plt.clf()      
gc.collect()
