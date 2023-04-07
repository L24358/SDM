import os
import numpy as np
from itertools import product

def H(x, a=270, b=108, d=0.154): return (a*x-b)/(1-2.71828**(-d*(a*x-b)))

def phi(x, a=615, b=177, d=0.087): return (a*x-b)/(1-2.71828**(-d*(a*x-b)))/2

def phiL(x, a=615):
    aeq, beq = a/2, -0.4425/0.005
    return aeq*x + beq

def init_x0(n, d=0.1, strt=0, end=1):
    l = np.arange(strt, end, d)
    return [np.array(item) for item in product(l, repeat=n)]

def datapath():
    return 'C:\\Users\\belle\\Documents\\03_SDM\\Sim07 Nullclines-3'
    
def datapath2():
    return 'D:\\15_DM\\Sim07 Nullclines-3'
    
def graphpath():
    return 'C:\\Users\\belle\\Documents\\03_SDM\\Sim09 Plotting'

def get_max(path, run):
    folders = os.listdir(path); maxx = -1
    for fold in folders:
        if ('run'+str(run) in fold) and ('.npy' in fold):
            new = int(fold.split('_')[1].split('.')[0][2:])
            if new > maxx: maxx = new
    return maxx
    
def remove_nan(x, y):
    newx, newy = [], []
    for i in range(len(x)):
        if np.isnan(x[i]) or np.isnan(y[i]): pass
        else: newx.append(x[i]); newy.append(y[i])
    return newx, newy
    
def default():
    return 0.2609, 0.0497, 0.3255
    
def rename(folder, irun, frun):
    ori = os.getcwd()
    os.chdir(folder)
    files = os.listdir(folder)
    for f in files:
        if 'run'+str(irun) in f:
            name = 'run'+str(frun)+f[4:]
            os.rename(f, name)
    os.chdir(ori)
    
def rn(num):
    if round(num) == round(num, 1): return int(num)
    elif round(num,1) == round(num,2): return round(num,1)
    return round(num,2)
    
def rpt(N, val): return [val]*N


# function to convert to subscript
def get_sub(x):
    normal = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-=()"
    sub_s = "ₐ₈CDₑբGₕᵢⱼₖₗₘₙₒₚQᵣₛₜᵤᵥwₓᵧZₐ♭꜀ᑯₑբ₉ₕᵢⱼₖₗₘₙₒₚ૧ᵣₛₜᵤᵥwₓᵧ₂₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎"
    res = x.maketrans(''.join(normal), ''.join(sub_s))
    return x.translate(res)