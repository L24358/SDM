import numpy as np
import matplotlib.pyplot as plt
import aadm.basics as bcs
from scipy.optimize import fsolve
from joblib import Parallel, delayed

from __main__ import pdic
params = ['giis','gii','Ii','sg2_max']
default = {}
for pm in params: default[pm] = 1
default.update(pdic)
for pm in params: exec(pm+' = default["'+pm+'"]')

def bunch(x0, tpe='ori'):
    
    func = bcs.phi if tpe == 'ori' else bcs.phiL
    def nullcline_UM(z):
        sg1 = z[0]
        F = np.empty((1))
        F[0] = -sg1 + func(-giis*sg1-gii*sg2+Ii)*0.005
        return F
    
    err = []; x = []; y = []
    for sg2 in np.arange(0, sg2_max, 0.01):
        sol = Parallel(n_jobs=10)(delayed(fsolve)(nullcline_UM, x0=p) for p in x0)
        for s in sol:
            er = abs(nullcline_UM(s)); err += list(er)
            if er[0] < 0.0000001: x.append(sg2); y.append(s)
    y = np.transpose(y)
    return x, y, err
