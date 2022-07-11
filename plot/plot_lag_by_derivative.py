'''
Plot deviation (lag) as a function of the derivative of splus, sminus.
'''

import os
import numpy as np
import matplotlib.pyplot as plt
import aadm.basics as bcs
import aadm.helper as hpr
import aadm.vision as vsn
from mpl_toolkits import mplot3d
from scipy.optimize import curve_fit

mother = os.getcwd()

def sim_sminus(params, ss=0.25, d=0.1):
    s1, s2, sg1, sg2 = vsn.sim_SM_nonoise(ss, params, d = d)
    splus = (s1+s2)/2; sgplus = (sg1+sg2)/2
    sminus = (s1-s2)/2; sgminus = (sg1-sg2)/2
    return splus, sgplus, sminus, sgminus, s1, s2, sg1, sg2

x, y, c = [], [], []
gma, gma2, gma3, gma4, alpha = [0.5,0.6,0.7,1,1]
params = hpr.get_params_s2(gma, gma2, gma3, gma4, alpha)

gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params
sp_inp = 307.5*geis*(1+gma2)
sp_rel = 200+307.5*giis*(1+gma4)
sm_inp = 307.5*geis*(1-gma2)
sm_rel = 200+307.5*giis*(1-gma4)

splus, sgplus, sminus, sgminus, s1, s2, sg1, sg2 = sim_sminus(params, ss=0.25)
s1, s2, sg1, sg2 = vsn.downsample([s1, s2, sg1, sg2])
ssplus = (sp_inp*splus-88.5+307.5)/sp_rel
ssminus = (sm_inp*sminus)/sm_rel
dp = sgplus - ssplus
dm = sgminus - ssminus

dm, dp, splus_down, sminus_down = vsn.downsample([dm, dp, splus, sminus])

dp_plot = -dp[:49]
dm_plot = -dm[:49]
splus_plot = splus_down[1:][:50]-splus_down[:-1][:50]
sminus_plot = sminus_down[1:][:50]-sminus_down[:-1][:50]

fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
ax.plot(dp_plot, 'b', label='deviation, +') #/np.std(dp_plot)
ax.plot(dm_plot, 'r', label='deviation, -') #/np.std(dm_plot)

ax.plot(splus_plot*(sp_inp/sp_rel**2)*100, color='cornflowerblue', linestyle='--', label='derivative, +')
ax.plot(sminus_plot*(sm_inp/sm_rel**2)*100, color='lightcoral', linestyle='--', label='derivative, -')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('t (ms)'); plt.ylabel('Normalized S'); plt.legend()
plt.tight_layout()

plt.savefig(os.path.join(mother, 'lag_as_derivative.png'), dpi=300)





