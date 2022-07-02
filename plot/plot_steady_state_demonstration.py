import os
import numpy as np
import matplotlib.pyplot as plt
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
from mpl_toolkits import mplot3d
from scipy.optimize import curve_fit

mother = os.getcwd()
gmavals = [[0.3,0.6,0.6,1,1], [0.3,0.6,0.7,1,1], [0.3,0.6,0.8,1,1], [0.3,0.8,1.25,1,1],\
           [0.5,0.6,0.7,1,1], [0.5,0.6,0.8,1,1], [0.3,0.6,0.5,1,1], [0.5,1,1,0,2]]

def sim_sminus(params, ss=0.25, d=0.01):
    s1, s2, sg1, sg2 = vsn.downsample(vsn.sim_SM_nonoise(ss, params, d = d))
    splus = (s1+s2)/2; sgplus = (sg1+sg2)/2
    sminus = (s1-s2)/2; sgminus = (sg1-sg2)/2
    return splus, sgplus, sminus, sgminus

x, y, c = [], [], []
for k in range(len(gmavals)):
    gma, gma2, gma3, gma4, alpha = gmavals[k]
    params = hpr.get_params_s2(gma, gma2, gma3, gma4, alpha)
    fname = str(gma)+', '+str(gma2)+', '+str(gma3)+', '+str(gma4)+', '+str(alpha)
    
    gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params
    sp_inp = 307.5*geis*(1+gma2)
    sp_rel = 200+307.5*giis*(1+gma4)
    sm_inp = 307.5*geis*(1-gma2)
    sm_rel = 200+307.5*giis*(1-gma4)
    
    splus, sgplus, sminus, sgminus = sim_sminus(params)
    ssplus = (sp_inp*splus-88.5+307.5)/sp_rel
    ssminus = (sm_inp*sminus)/sm_rel
    dp = sgplus - ssplus
    dm = sgminus - ssminus
    
    





