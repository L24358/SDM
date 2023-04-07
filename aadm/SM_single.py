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
        s1 = z[0]; sg1 = z[1]; sg2 = z[2]
        F = np.empty((3))
        F[0] = -s1/0.1 + (1-s1)*0.641*bcs.H(gees*s1+gee*s2-gies*sg1-gie*sg2+Ie)
        F[1] = -sg1 + func(geis*s1+gei*s2-giis*sg1-gii*sg2+Ii)*0.005
        F[2] = -sg2 + func(geis*s2+gei*s1-giis*sg2-gii*sg1+Ii)*0.005
        return F
    
    err = []; x = []; y = []
    for s2 in list(np.arange(0, 0.1, 0.005)):#+list(np.arange(0.1, 1, 0.05)):
        sol = Parallel(n_jobs=10)(delayed(fsolve)(nullcline_UM, x0=p) for p in x0)
        err = [abs(nullcline_UM(s))[0] for s in sol]
        y.append(sol[err.index(min(err))]); x.append(s2)
    y = np.transpose(y)
    return x, y, err
