import numpy as np
import aadm.helper as hpr

def path1(d):
    '''
    Returns:
        - Array of coordinates (S_E1, S_E2, S_I1, S_I2)
        - Array of dx values
    '''
    steps = np.arange(0, 0.7, d)
    zeros = np.zeros(len(steps))
    

params = hpr.get_params_s2(0.5, 0.6, 0.7, 1, 1, giis=1.2, Ii=1)
path, dx = path1()

integral = 0
for i in range(len(path)):
    f = hpr.rhs_s1(0.6, path[i], params, 0, n1=0, n2=0)
    integral += np.array(f).dot(dx[i])