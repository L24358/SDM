import numpy as np
import matplotlib.pyplot as plt
import aadm.basics as bcs
from scipy.optimize import fsolve
from joblib import Parallel, delayed

from __main__ import pdic
params = ['gees','gee','geis','gei','gies','gie','giis','gii','Ie','Ii']
default = {}
for pm in params: default[pm] = 1
default.update(pdic)
for pm in params: exec(pm+' = default["'+pm+'"]')

def bunch(x0, tpe='ori'):
    
    func = bcs.phi if tpe == 'ori' else bcs.phiL
    def nullcline_UM(z):
        s1 = z[0]; sg1 = z[1]
        F = np.empty((2))
        F[0] = -s1/0.1 + (1-s1)*0.641*bcs.H(gees*s1+gee*s2-gie*sg1+Ie)
        F[1] = -sg1 + func(gei*s1+gei*s2-gii*sg1+Ii)*0.005
        return F
    
    err = []; x = []; y = []
    for s2 in np.arange(0, 1, 0.01):
        sol = Parallel(n_jobs=10)(delayed(fsolve)(nullcline_UM, x0=p) for p in x0)
        for s in sol:
            er = abs(nullcline_UM(s)); err += list(er)
            if er[0] < 0.0000001: x.append(s2); y.append(s)
    y = np.transpose(y)
    return x, y, err
 