# -*- coding: utf-8 -*-
"""
Plot velocity magnitude as a contour plot.

* Notes:
- https://jakevdp.github.io/PythonDataScienceHandbook/04.04-density-and-contour-plots.html
"""

import numpy as np
import matplotlib.pyplot as plt
import aadm.helper as hpr

x = np.linspace(0, 0.6, 100)
y = np.linspace(0, 0.6, 100)
X, Y = np.meshgrid(x, y)
gma, gma2, gma3, gma4, alpha = [0.5,0.6,0.7,1,1]
params = hpr.get_params_s2(gma, gma2, gma3, gma4, alpha)

def f(X, Y):
    Z = np.zeros(X.shape)
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            s1, s2 = X[i][j], Y[i][j]
            sg = hpr.get_ssSG_s1(s1, s2, params)
            v = [s1, s2] + list(sg)
            f1, f2, _, _ = hpr.rhs_s1(0.6, v, params, 0)
            Z[i][j] = np.linalg.norm([f1, f2])
    return Z

Z = f(X, Y)
plt.contourf(X, Y, Z, 20, cmap='RdGy')
plt.colorbar();