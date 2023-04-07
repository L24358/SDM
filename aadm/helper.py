import os
import pickle
import numpy as np
import aadm.basics as bcs
import dynalysis.basics as dbcs
import dynalysis.classes as clss
from scipy.optimize import fsolve, curve_fit
from joblib import Parallel, delayed
from scipy.signal import find_peaks
from scipy.integrate import solve_ivp
from math import sqrt

'''Basic Functions'''

def pulse(A, t): return A*np.heaviside(t, 1)

def pulse2(A, t1, t2): return A*np.heaviside(t1, 1)*np.heaviside(t2, 1)

def Inp1(c, mu0=0.0156): return mu0*(1+c/100)

def Inp2(c, mu0=0.0156): return mu0*(1-c/100)

def S_N(s, Isyn): return -s/0.1 + (1-s)*0.641*bcs.H(Isyn)

def S_G(s, Isyn): return (-s/0.005 + bcs.phi(Isyn))

def noise(mu, sigma): return np.random.normal(mu, sigma)

def corrnoise(mu, sigma): return np.random.multivariate_normal(mu, sigma)

def affine(x, a, b): return a*x + b

def cubic(x, a, b, c, d): return a*x**3 + b*x**2 + c*x + d

def reaction_time(y, t0=50, thre=0.5):
    r = y[t0:]
    for t in range(len(r)):
        if r[t] > thre: return t*0.01
    return 9999999

def accuracy(y1, y2, t0=50, thre=0.5):
    rt1 = reaction_time(y1, t0=t0, thre=thre)
    rt2 = reaction_time(y2, t0=t0, thre=thre)
    rt = min(rt1, rt2)
    if rt < 10 and rt1 == rt: return rt, 1
    elif rt < 10 and rt2 == rt: return rt, 2
    else: return None, None
    
def largest(l, judge):
    judge = list(judge)
    one, two = sorted(judge, reverse=True)[:2]
    ind1, ind2 = judge.index(one), judge.index(two)
    return [l[ind1], l[ind2]]

def get_flow(corrs):
    dists = []
    for t in range(len(corrs)-1):
        dist = np.linalg.norm(np.array(corrs[t+1])-np.array(corrs[t]))
        dists.append(dist)
    return dists

def get_peaks(dist, prom=0.001):
    pks, pdic = find_peaks(dist, prominence=prom, height=0.002)
    if len(pks)==2: return pks, pdic['peak_heights']
    elif len(pks)<2: return get_peaks(dist, prom=prom/2)
    else:
        pns = pdic['prominences']
        pks, phs = largest(pks, pns), largest(pdic['peak_heights'], pns)
        return pks, phs
    
def get_valley(dist, pks):
    try:
        dist = list(dist)
        return dist.index(min(dist)), min(dist)
    except:
        print(pks)
    
def to_zero(dist):
    for d in range(len(dist)):
        if dist[d] < 0.002: return d
    return len(dist)

def get_velocity(coors):
    return [np.array(coors[c+1])-np.array(coors[c]) for c in range(len(coors)-1)]
    
def get_coors(x, y):
    dic = {}
    for i in range(len(x)): dic[x[i]] = y[i]
    nx = sorted(x); ny = [dic[j] for j in nx]
    return list(zip(nx, ny))
    
def check_distance(coors):
    vel = get_velocity(coors)
    mag = [np.linalg.norm(v) for v in vel]
    if min(mag) < 0.03: return False
    return True

def check_region(coors):
    for c in coors:
        if sum(c) > 0.7: return False
    return True
    
'''Get params, EQ, nullclines, ssSG, flow, Jacobian'''

def get_params_u1(gma, gii=1.2, Ii=1):
    def param_map(v):
        a = 1.5375; b = -0.4425
        ges, gie, inpe = v
        f1 = ges - (a*gie**2)/(1 + a*gii) - 0.2609
        f2 = gma*ges - (a*gie**2)/(1 + a*gii) + 0.0497
        f3 = inpe - (gie*a)/(1 + a*gii)*Ii - 0.3255 - gie*b/(1 + a*gii)
        return [f1,f2,f3]
    
    gees, gie, Ie = fsolve(param_map, [0.5,0.5,0.5])
    return gees, gees*gma, gie, gie, gii, Ie, Ii

def get_params_u2(gma, gei, gii, Ii=1):
    def param_map(v):
        a = 1.5375; b = -0.4425
        ges, gie, inpe = v
        f1 = ges - (a*gie*gei)/(1 + a*gii) - 0.2609
        f2 = gma*ges - (a*gie*gei)/(1 + a*gii) + 0.0497
        f3 = inpe - (gie*a)/(1 + a*gii)*Ii - 0.3255 - gie*b/(1 + a*gii)
        return [f1,f2,f3]
    
    gees, gie, Ie = fsolve(param_map, [0.5,0.5,0.5])
    return gees, gees*gma, gei, gie, gii, Ie, Ii
    
def get_params_u3(gma, alpha, gii, Ii, Je=0.2609, Ji=0.0497, Ib=0.3255):
    def param_map(v):
        a = 1.5375; b = -0.4425
        ges, gie, inpe = v
        f1 = ges - (a*gie*gie/alpha)/(1 + a*gii) - Je
        f2 = gma*ges - (a*gie*gie/alpha)/(1 + a*gii) + Ji
        f3 = inpe - (gie*a)/(1 + a*gii)*Ii - Ib - gie*b/(1 + a*gii)
        return [f1,f2,f3]
    
    gees, gie, Ie = fsolve(param_map, [0.5,0.5,0.5])
    return gees, gees*gma, gie/alpha, gie, gii, Ie, Ii

def get_params_s2(gma, gma2, gma3, gma4, alpha, giis=1.2, Ii=1):
    gii = giis*gma4
    def param_map(v):
        a = 1.5375; b = -0.4425
        di = 1 + a*giis - (a*gii)**2/(1 + a*giis)
        gees, gies, inpe = v
        f1 = gees - gies/di*(a*gies*alpha - (a**2*gii*gies*gma2*alpha)/(1 + a*giis)) - \
            gies*gma3/di*(a*gies*gma2*alpha - (a**2*gii*gies*alpha)/(1 + a*giis)) - 0.2609
        f2 = gees*gma - gies/di*(a*gies*gma2*alpha - (a**2*gii*gies*alpha)/(1 + a*giis)) - \
            gies*gma3/di*(a*gies*alpha - (a**2*gii*gies*gma2*alpha)/(1 + a*giis)) + 0.0497
        f3 = (-gies - gies*gma3)/di*(b + a*Ii)*(1 - a*gii/(1 + a*giis)) + \
        inpe - 0.3255
        return [f1,f2,f3]
    
    gees, gies, Ie = fsolve(param_map, [0.5,0.5,0.5])
    geis = gies*alpha
    gee = gees*gma
    gei = geis*gma2
    gie = gies*gma3
    return gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii
    
def get_params_s3(gma, gma2, gma3, gma4, alpha, giis, Ii, Je, Ji, Ib):
    gii = giis*gma4
    def param_map(v):
        a = 1.5375; b = -0.4425
        di = 1 + a*giis - (a*gii)**2/(1 + a*giis)
        gees, gies, inpe = v
        f1 = gees - gies/di*(a*gies*alpha - (a**2*gii*gies*gma2*alpha)/(1 + a*giis)) - \
            gies*gma3/di*(a*gies*gma2*alpha - (a**2*gii*gies*alpha)/(1 + a*giis)) - Je
        f2 = gees*gma - gies/di*(a*gies*gma2*alpha - (a**2*gii*gies*alpha)/(1 + a*giis)) - \
            gies*gma3/di*(a*gies*alpha - (a**2*gii*gies*gma2*alpha)/(1 + a*giis)) + Ji
        f3 = (-gies - gies*gma3)/di*(b + a*Ii)*(1 - a*gii/(1 + a*giis)) + \
        inpe - Ib
        return [f1,f2,f3]
    
    gees, gies, Ie = fsolve(param_map, [0.5,0.5,0.5])
    geis = gies*alpha
    gee = gees*gma
    gei = geis*gma2
    gie = gies*gma3
    return gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii

def adjust_params(params, fee, fei, fie, fii):
    gees, gee, gei, gie, gii, Ie, Ii = params
    gees = gees*fee; gee=gee*fee
    gei = gei*fei; gie = gie*fie
    gii = gii*fii
    return gees, gee, gei, gie, gii, Ie, Ii

def get_nullcline_u1(c, sG, t, tpe, params, minn=0, maxx=1, step=0.01):
    gees, gee, gei, gie, gii, Ie, Ii = params
    x0 = bcs.init_x0(1, d=0.2)
    func = Inp1 if tpe==1 else Inp2
    def rhs(v):
        s1 = v
        f1 = S_N(s1, gees*s1+gee*s2-gie*sG+Ie+pulse(func(c), t-0.5))
        return f1
    
    err = []; x = []; y = []
    for s2 in np.arange(minn, maxx, step):
        sol = Parallel(n_jobs=10)(delayed(fsolve)(rhs, x0=p) for p in x0)
        for s in sol:
            er = [abs(i) for i in rhs(s)]; err += list(er)
            if er[0] < 0.0000001: x.append(s2); y.append(s)
    y = np.transpose(y)
    return x, y, err
    
def get_nullcline_u2(c, t, x0, params, tpe, tpe2='ori', minn=0, maxx=1, step=0.01, mu0=0.0156):
    gees, gee, gei, gie, gii, Ie, Ii = params
    func = Inp1 if tpe==1 else Inp2
    func2 = bcs.phi if tpe2 == 'ori' else bcs.phiL
    def nullcline_UM(z):
        s1 = z[0]; sg1 = z[1]
        F = np.empty((2))
        F[0] = -s1/0.1 + (1-s1)*0.641*bcs.H(gees*s1+gee*s2-gie*sg1+Ie+pulse(func(c, mu0=mu0), t-0.5))
        F[1] = -sg1 + func2(gei*s1+gei*s2-gii*sg1+Ii)*0.005
        return F
    
    err = []; x = []; y = []
    for s2 in np.arange(minn, maxx, step):
        sol = Parallel(n_jobs=10)(delayed(fsolve)(nullcline_UM, x0=p) for p in x0)
        for s in sol:
            er = abs(nullcline_UM(s)); err += list(er)
            if er[0] < 0.0000001: x.append(s2); y.append(s)
    y = np.transpose(y)
    return x, y, err
    
def get_nullcline_two(c, t, x0, params, tpe, tpe2='ori', minn=0, maxx=1, step=0.01):
    Je, Ji, Ib = params
    func = Inp1 if tpe==1 else Inp2
    def nullcline_UM(z):
        s1 = z[0]
        F = np.empty((1))
        F[0] = -s1/0.1 + (1-s1)*0.641*bcs.H(Je*s1-Ji*s2+Ib+pulse(func(c), t-0.5))
        return F
    
    err = []; x = []; y = []
    for s2 in np.arange(minn, maxx, step):
        sol = Parallel(n_jobs=10)(delayed(fsolve)(nullcline_UM, x0=p) for p in x0)
        for s in sol:
            er = abs(nullcline_UM(s)); err += list(er)
            if er[0] < 0.0000001: x.append(s2); y.append(s)
    y = np.transpose(y)
    return x, y, err

def get_nullcline_s2(c, t, x0, params, tpe, tpe2='ori', minn=0, maxx=1, step=0.01):
    gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params
    func = Inp1 if tpe==1 else Inp2
    func2 = bcs.phi if tpe2 == 'ori' else bcs.phiL
    def nullcline_UM(z):
        s1 = z[0]; sg1 = z[1]; sg2 = z[2]
        F = np.empty((3))
        F[0] = -s1/0.1 + (1-s1)*0.641*bcs.H(gees*s1+gee*s2-gies*sg1-gie*sg2+Ie+pulse(func(c), t-0.5))
        F[1] = -sg1 + func2(geis*s1+gei*s2-giis*sg1-gii*sg2+Ii)*0.005
        F[2] = -sg2 + func2(geis*s2+gei*s1-giis*sg2-gii*sg1+Ii)*0.005
        return F
    
    err = []; x = []; y = []
    for s2 in np.arange(minn, maxx, step):
        sol = Parallel(n_jobs=10)(delayed(fsolve)(nullcline_UM, x0=p) for p in x0)
        for s in sol:
            er = abs(nullcline_UM(s)); err += list(er)
            if er[0] < 0.0000001: x.append(s2); y.append(s)
    y = np.transpose(y)
    return x, y, err
    
def sort_nullcline(s1, s2):
    dic ={}
    for i in range(len(s2)): dic[s2[i]] = s1[i]
    new_s2 = sorted(s2)
    new_s1 = [dic[new_s2[i]] for i in range(len(s1))]
    return new_s1, new_s2

def sim_UM(params, init=None, c=0, nsig=0.02, tA=0.002, snum=0, t0=0.5):
    gees, gee, gei, gie, gii, Ie, Ii = params
    def rhs(t, v):
        s1, s2, sG, n1, n2 = v
        f1 = S_N(s1, gees*s1+gee*s2-gie*sG+Ie+pulse(Inp1(c), t-t0)+n1)
        f2 = S_N(s2, gees*s2+gee*s1-gie*sG+Ie+pulse(Inp2(c), t-t0)+n2)
        f3 = S_G(sG, gei*s1+gei*s2-gii*sG+Ii)
        f4 = (-n1 + noise(0,1)*sqrt(tA*nsig**2))/tA
        f5 = (-n2 + noise(0,1)*sqrt(tA*nsig**2))/tA
        return [f1,f2,f3,f4,f5]
    
    if init==None: init = list(get_EQ_u1(params))+[0,0]
    res = solve_ivp(rhs, (0,3), init, t_eval=np.arange(0,3,0.01))
    return res.y

def get_EQ_u1(params):
    gees, gee, gei, gie, gii, Ie, Ii = params
    def rhs(v):
        s1, s2, sG = v
        f1 = S_N(s1, gees*s1+gee*s2-gie*sG+Ie)
        f2 = S_N(s2, gees*s2+gee*s1-gie*sG+Ie)
        f3 = S_G(sG, gei*s1+gei*s2-gii*sG+Ii)
        return [f1,f2,f3]
    ssvec = fsolve(rhs, [0.1, 0.1, 0.1])
    return ssvec

def get_EQ_s1(params):
    gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params
    def rhs(v):
        s1, s2, sG1, sG2 = v
        f1 = S_N(s1, gees*s1+gee*s2-gies*sG1-gie*sG2+Ie)
        f2 = S_N(s2, gees*s2+gee*s1-gies*sG2-gie*sG1+Ie)
        f3 = S_G(sG1, geis*s1+gei*s2-giis*sG1-gii*sG2+Ii)
        f4 = S_G(sG2, geis*s2+gei*s1-giis*sG2-gii*sG1+Ii)
        return [f1,f2,f3,f4]
    ssvec = fsolve(rhs, [0.1, 0.1, 0.1, 0.1])
    return ssvec

def get_EQ_ab(params, a, b):
    gees, gee, gei, gie, gii, Ie, Ii = params
    def rhs(v):
        s1, s2 = v
        f1 = S_N(s1, gees*s1+gee*s2-gie*(a*(s1+s2)/2+b)+Ie)
        f2 = S_N(s2, gees*s2+gee*s1-gie*(a*(s1+s2)/2+b)+Ie)
        return [f1,f2]
    ssvec = fsolve(rhs, [0.1, 0.1])
    return ssvec

def get_EQ_ab_s1(params, a, b):
    gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params
    m = gies+gie
    def rhs(v):
        s1, s2, sm = v
        f1 = S_N(s1, gees*s1+gee*s2-m*(a*(s1+s2)+b)/2+Ie)
        f2 = S_N(s2, gees*s2+gee*s1-m*(a*(s1+s2)+b)/2+Ie)
        f3 = 615/2*((geis-gei)*(s1-s2)/2-(giis-gii)*sm)-sm/0.005
        return [f1,f2,f3]
    ssvec = fsolve(rhs, [0.1, 0.1, 0.1])
    return ssvec
    
def get_EQ_ab_s2(params, a, b):
    gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params
    m = gies+gie
    def rhs(v):
        s1, s2 = v
        f1 = S_N(s1, gees*s1+gee*s2-m*(a*(s1+s2)+b)/2+Ie)
        f2 = S_N(s2, gees*s2+gee*s1-m*(a*(s1+s2)+b)/2+Ie)
        return [f1,f2]
    ssvec = fsolve(rhs, [0.1, 0.1])
    return ssvec
    
def get_EQ_ab_s3(params, a, b, ag, init=[0.1,0.1], mu=0, bg=0):
    gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params
    m, n = gies+gie, gies-gie
    def rhs(v):
        s1, s2 = v
        f1 = S_N(s1, gees*s1+gee*s2-m*(a*(s1+s2)+b)/2-n*(ag*(s1-s2)/2+bg)+Ie+mu)
        f2 = S_N(s2, gees*s2+gee*s1-m*(a*(s1+s2)+b)/2+n*(ag*(s1-s2)/2+bg)+Ie+mu)
        return [f1,f2]
    ssvec = fsolve(rhs, init)
    return ssvec
    
def get_EQ_ab_s4(params, a, b, ag, bg, init=[0.1,0.1], mu=0):
    gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params
    m, n = gies+gie, gies-gie
    def rhs(v):
        s1, s2 = v
        f1 = S_N(s1, gees*s1+gee*s2-m*(a*(s1+s2)+b)/2-n*(ag*((s1-s2)/2)**2+bg*((s1-s2)/2))+Ie+mu)
        f2 = S_N(s2, gees*s2+gee*s1-m*(a*(s1+s2)+b)/2+n*(ag*((s1-s2)/2)**2+bg*((s1-s2)/2))+Ie+mu)
        return [f1,f2]
    ssvec = fsolve(rhs, init)
    return ssvec

def get_ssSG(s1, s2, params):
    gees, gee, gei, gie, gii, Ie, Ii = params
    def rhs(v):
        f3 = S_G(v, gei*s1+gei*s2-gii*v+Ii)
        return f3
    ssvec = fsolve(rhs, [0.1])
    return ssvec

def get_ssSG_s1(s1, s2, params):
    gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params
    def rhs(v):
        sG1, sG2 = v
        f3 = S_G(sG1, geis*s1+gei*s2-giis*sG1-gii*sG2+Ii)
        f4 = S_G(sG2, geis*s2+gei*s1-giis*sG2-gii*sG1+Ii)
        return f3, f4
    ssvec = fsolve(rhs, [0.1, 0.1])
    return ssvec

def rhs(t, v, params, c, n1=0, n2=0):
    gees, gee, gei, gie, gii, Ie, Ii = params
    s1, s2, sG = v
    f1 = S_N(s1, gees*s1+gee*s2-gie*sG+Ie+pulse(Inp1(c), t-0.5)+n1)
    f2 = S_N(s2, gees*s2+gee*s1-gie*sG+Ie+pulse(Inp2(c), t-0.5)+n2)
    f3 = S_G(sG, gei*s1+gei*s2-gii*sG+Ii)
    return [f1,f2,f3]
    

def rhs_s1(t, v, params, c, n1=0, n2=0):
    gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params
    s1, s2, sG1, sG2 = v
    f1 = S_N(s1, gees*s1+gee*s2-gies*sG1-gie*sG2+Ie+pulse(Inp1(c), t-0.5)+n1)
    f2 = S_N(s2, gees*s2+gee*s1-gies*sG2-gie*sG1+Ie+pulse(Inp2(c), t-0.5)+n2)
    f3 = S_G(sG1, geis*s1+gei*s2-giis*sG1-gii*sG2+Ii)
    f4 = S_G(sG2, geis*s2+gei*s1-giis*sG2-gii*sG1+Ii)
    return [f1,f2,f3,f4]

def get_Jacobian(params, c):
    def rhs_v2(v):
        s1, s2, sG = v
        f1 = S_N(s1, gees*s1+gee*s2-gie*sG+Ie+Inp1(c))
        f2 = S_N(s2, gees*s2+gee*s1-gie*sG+Ie+Inp2(c))
        f3 = S_G(sG, gei*s1+gei*s2-gii*sG+Ii)
        return [f1,f2,f3]
    gees, gee, gei, gie, gii, Ie, Ii = params
    ss1, ss2, _ = fsolve(rhs_v2, [0.1, 0.1, 0.1])
    Hs1 = ss1/(0.641*0.1*(1-ss1))
    Hs2 = ss2/(0.641*0.1*(1-ss2))
    
    def modH(x): return bcs.H(x)-Hs
    Hs = Hs1
    xs1 = fsolve(modH, [1])[0]
    Hs = Hs2
    xs2 = fsolve(modH, [1])[0]
    
    def derH(x, a=270, b=108, d=0.154):
        exp = 2.71828**(-d*(a*x-b))
        return a/(1-exp) - a*d*(a*x-b)*exp/(1-exp)**2
    ders1 = derH(xs1)
    ders2 = derH(xs2)
    
    J = [[-1/0.1-0.641*Hs1+(1-ss1)*0.641*ders1*gees, (1-ss1)*0.641*ders1*gee, -(1-ss1)*0.641*ders1*gie],\
         [(1-ss2)*0.641*ders2*gee, -1/0.1-0.641*Hs2+(1-ss2)*0.641*ders2*gees, -(1-ss2)*0.641*ders2*gie],\
         [615/2*gei, 615/2*gei, -1/0.005-615/2*gii]]
    P = [[1,1,0], [1,-1,0], [0,0,1]]
    Q = np.matmul(np.matmul(P, J), np.linalg.inv(P))
    w, v = np.linalg.eig(Q)
    
    return w, np.transpose(v)

def get_Jacobian_s1(params, c):
    def rhs_v2(v):
        s1, s2, sG1, sG2 = v
        f1 = S_N(s1, gees*s1+gee*s2-gies*sG1-gie*sG2+Ie+Inp1(c))
        f2 = S_N(s2, gees*s2+gee*s1-gies*sG2-gie*sG1+Ie+Inp2(c))
        f3 = S_G(sG1, geis*s1+gei*s2-giis*sG1-gii*sG2+Ii)
        f4 = S_G(sG2, geis*s2+gei*s1-giis*sG2-gii*sG1+Ii)
        return [f1,f2,f3,f4]
    gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params
    ss1, ss2, _, _ = fsolve(rhs_v2, [0.1, 0.1, 0.1, 0.1])
    Hs1 = ss1/(0.641*0.1*(1-ss1))
    Hs2 = ss2/(0.641*0.1*(1-ss2))
     
    def modH(x): return bcs.H(x)-Hs
    Hs = Hs1
    xs1 = fsolve(modH, [1])[0]
    Hs = Hs2
    xs2 = fsolve(modH, [1])[0]
    
    def derH(x, a=270, b=108, d=0.154):
        exp = 2.71828**(-d*(a*x-b))
        return a/(1-exp) - a*d*(a*x-b)*exp/(1-exp)**2
    
    ders1 = derH(xs1)
    ders2 = derH(xs2)
    
    J = [[-1/0.1-0.641*Hs1+(1-ss1)*0.641*ders1*gees, (1-ss1)*0.641*ders1*gee, -(1-ss1)*0.641*ders1*gies, -(1-ss1)*0.641*ders1*gie],\
         [(1-ss2)*0.641*ders2*gee, -1/0.1-0.641*Hs2+(1-ss2)*0.641*ders2*gees, -(1-ss2)*0.641*ders2*gie, -(1-ss2)*0.641*ders2*gies],\
         [615/2*geis, 615/2*gei, -1/0.005-615/2*giis, -615/2*gii],\
         [615/2*gei, 615/2*geis, -615/2*gii, -1/0.005-615/2*giis]]
    P = [[1,1,0,0], [1,-1,0,0], [0,0,1,1], [0,0,1,-1]]
    Q = np.matmul(np.matmul(P, J), np.linalg.inv(P))
    w, v = np.linalg.eig(Q)
    return w, np.transpose(v)

def get_slope_intercept(params):
    s1, s2, sg, _, _ = sim_UM(params, nsig=0)
    return curve_fit(affine, (s1+s2)/2, sg)[0]
    
def same(lst, x):
    for item in lst:
        deter = np.any([round(x[i],3)==round(item[i],3) for i in range(2)])
        if deter: return True
    return False
    
def deter_root_count(params, inp1=0, inp2=0):
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

'''Get data related'''

def inquire_trace(run, fname, rnum, seed, mother=bcs.datapath(), snum=0, rdm=False):
    if seed !=None: np.random.seed(seed)
    b = clss.branch(fname, mother); os.chdir(b.pathlink)
    print(b.pathlink)
    enum = bcs.get_max(b.pathlink, run)
    if enum < rnum: rnum = enum+1
    if rdm: rps = np.random.choice(range(snum, enum), rnum)
    else: rps = range(snum, rnum)
    traces = []
    for rp in rps:
        sim = pickle.load(open('run'+str(run)+'_tr'+str(rp)+'.p', 'rb'))
        resy = sim['res.y']; traces.append(resy)
    os.chdir(mother)
    return traces

def slowpoke(sim, params):
    gees, gee, gei, gie, gii, Ie, Ii = params
    diffs = []; sgs, ssgs = [], []
    for rp in range(len(sim)):
        s1, s2, sg = sim[rp][0], sim[rp][1], sim[rp][2]
        ssg_vec = [get_ssSG(s1[t], s2[t], params)[0] for t in range(len(s1))]
        diffs.append(np.array(sg)-np.array(ssg_vec))
        sgs.append(np.array(sg)); ssgs.append(np.array(ssg_vec))
    return np.mean(diffs, axis=0), np.mean(sgs, axis=0), np.mean(ssgs, axis=0), diffs

def slowpoke_s1(sim, params):
    gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params
    diffs1, diffs2 = [], []; ssgs = []
    for rp in range(len(sim)):
        s1, s2, sg1, sg2 = sim[rp][0], sim[rp][1], sim[rp][2], sim[rp][3]
        ssg1 = [get_ssSG_s1(s1[t], s2[t], params)[0] for t in range(len(s1))]
        ssg2 = [get_ssSG_s1(s1[t], s2[t], params)[1] for t in range(len(s1))]
        diffs1.append(np.array(sg1)-np.array(ssg1))
        diffs2.append(np.array(sg2)-np.array(ssg2))
    return np.mean(diffs1, axis=0), np.mean(diffs2, axis=0), np.mean(ssgs, axis=0), diffs1, diffs2
    
def get_flow_bunch(sim):
    vels = []
    for rp in range(len(sim)):
        s1, s2 = sim[rp][0], sim[rp][1]
        dist = get_flow(list(zip(s1,s2)))
        vels.append(dist)
    return vels
    
def get_period(sim, strt=0.5):
    T1s, T2s, T3s, T4s = [], [], [], []
    lc1, lc2, lc3 = [], [], []
    for rp in range(len(sim)):
        s1, s2 = sim[rp][0], sim[rp][1]
        dist = get_flow(list(zip(s1,s2)))
        try:
            pks, phs = get_peaks(dist); pk1, pk2 = sorted(pks)
            val, pvs = get_valley(dist[pk1:pk2], pks)
            T1, T2, T3, T4 = pk1, val, pk2-(pk1+val), to_zero(dist[pk2:])
            T1s.append(T1/100-strt); T2s.append(T2/100); T3s.append(T3/100); T4s.append(T4/100)
            lc1.append(s1[pk1]+s2[pk1]); lc2.append(s1[pk1+val]+s2[pk1+val]); lc3.append(s1[pk2]+s2[pk2])
        except: pass
    return np.mean(T1s), np.mean(T2s), np.mean(T3s), np.mean(T4s), np.mean(lc1)/2, np.mean(lc2)/2, np.mean(lc3)/2
    
def get_period_complete(sim, strt=0.5):
    T1s, T2s, T3s, T4s = [], [], [], []
    lc1, lc2, lc3 = [], [], []
    for rp in range(len(sim)):
        s1, s2 = sim[rp][0], sim[rp][1]
        dist = get_flow(list(zip(s1,s2)))
        try:
            pks, phs = get_peaks(dist); pk1, pk2 = sorted(pks)
            val, pvs = get_valley(dist[pk1:pk2], pks)
            T1, T2, T3, T4 = pk1, val, pk2-(pk1+val), to_zero(dist[pk2:])
            T1s.append(T1/100-strt); T2s.append(T2/100); T3s.append(T3/100); T4s.append(T4/100)
            lc1.append(s1[pk1]+s2[pk1]); lc2.append(s1[pk1+val]+s2[pk1+val]); lc3.append(s1[pk2]+s2[pk2])
        except: pass
    return T1s, T2s, T3s, T4s, np.array(lc1)/2, np.array(lc2)/2, np.array(lc3)/2
    
def get_period_complete_withNone(sim, strt=0.5):
    print('This function returns result with None.')
    print('The locations are not divided by 2.')
    T1s, T2s, T3s, T4s = [], [], [], []
    lc1, lc2, lc3 = [], [], []
    for rp in range(len(sim)):
        s1, s2 = sim[rp][0], sim[rp][1]
        dist = get_flow(list(zip(s1,s2)))
        try:
            pks, phs = get_peaks(dist); pk1, pk2 = sorted(pks)
            val, pvs = get_valley(dist[pk1:pk2], pks)
            T1, T2, T3, T4 = pk1, val, pk2-(pk1+val), to_zero(dist[pk2:])
            T1s.append(T1/100-strt); T2s.append(T2/100); T3s.append(T3/100); T4s.append(T4/100)
            lc1.append(s1[pk1]+s2[pk1]); lc2.append(s1[pk1+val]+s2[pk1+val]); lc3.append(s1[pk2]+s2[pk2])
        except:
            T1s.append(None); T2s.append(None); T3s.append(None); T4s.append(None)
            lc1.append(None); lc2.append(None); lc3.append(None)
    return T1s, T2s, T3s, T4s, np.array(lc1), np.array(lc2), np.array(lc3)

def get_period_ratio(sim, strt=0.5):
    T1s, T2s, T3s = [], [], []
    for rp in range(len(sim)):
        s1, s2, sg1, sg2, _, _ = sim[rp]
        dist = get_flow(list(zip(s1,s2)))
        try:
            pks, phs = get_peaks(dist); pk1, pk2 = sorted(pks)
            val, pvs = get_valley(dist[pk1:pk2], pks)
            T0 = 0
            T1 = int(pk1+strt*100) #accumulated
            T2 = int(val+pk1+strt*100) #accumulated
            T3 = int(pk2+strt*100) #accumulated
            splus = (s1+s2)/2; sminus = (s1-s2)/2
            sign = np.sign(np.mean(sminus))
            T1s.append(np.mean((sminus/splus)[T0:T1])*sign)
            T2s.append(np.mean((sminus/splus)[T1:T2])*sign)
            T3s.append(np.mean((sminus/splus)[T2:T3])*sign)
        except: pass
    return np.mean(T1s), np.mean(T2s), np.mean(T3s)

def get_period_ratio_fixed(sim, strt=0.5):
    T1s, T2s, T3s = [], [], []
    for rp in range(len(sim)):
        s1, s2, sg1, sg2, _, _ = sim[rp]
        dist = get_flow(list(zip(s1,s2)))
        try:
            pks, phs = get_peaks(dist); pk1, pk2 = sorted(pks)
            val, pvs = get_valley(dist[pk1:pk2], pks)
            T0 = 0
            T2 = int(30+strt*100) #accumulated
            T3 = int(70+strt*100) #accumulated
            splus = (s1+s2)/2; sminus = (s1-s2)/2
            sign = np.sign(np.mean(sminus))
            T2s.append(np.mean((sminus/splus)[T0:T2])*sign)
            T3s.append(np.mean((sminus/splus)[T2:T3])*sign)
        except: pass
    return np.mean(T2s), np.mean(T3s)

def get_peak_index(sim):
    pk1s, pk2s, vals = [], [], []
    for rp in range(len(sim)):
        s1, s2 = sim[rp][0], sim[rp][1]
        dist = get_flow(list(zip(s1,s2)))
        try:
            pks, phs = get_peaks(dist); pk1, pk2 = sorted(pks)
            val, pvs = get_valley(dist[pk1:pk2], pks)
            pk1s.append(pk1); pk2s.append(pk2); vals.append(pk1+val)
        except: pk1s.append(None); pk2s.append(None); vals.append(None)
    return pk1s, pk2s, vals

def psych_func(sims, penalty=2.5, thre=0.5, suppress=False):
    accs, rtcs, rtws = [], [], []
    accstd, rtcstd, rtwstd = [], [], []
    for sim in sims:
        time_corr, time_wrong = [], []; winner = []
        for rp in range(len(sim)):
            s1, s2 = sim[rp][0], sim[rp][1]
            rt, win = accuracy(s1, s2, thre=thre)
            if rt != None:
                winner.append(win)
                if win==1: time_corr.append(rt)
                else: time_wrong.append(rt)
            else:
                winner.append(np.random.choice([1,2],1)[0])
                if win==1: time_corr.append(penalty)
                else: time_wrong.append(penalty)
                if not suppress: print('Undecided!')
        
        acc = winner.count(1)/len(winner); accs.append(acc)
        rtc = np.mean(time_corr); rtcs.append(rtc)
        rtw = np.mean(time_wrong); rtws.append(rtw)
        accstd.append(np.std(winner)/np.sqrt(len(winner)))
        rtcstd.append(np.std(time_corr)/np.sqrt(len(time_corr)))
        rtwstd.append(np.std(time_wrong)/np.sqrt(len(time_wrong)))
    return accs, rtcs, rtws, accstd, rtcstd, rtwstd
    
def psych_func_complete(sims, penalty=2.5, thre=0.5, t0=50):
    accs, rtcs = [], []
    for sim in sims:
        time_corr = []; winner = []
        for rp in range(len(sim)):
            s1, s2 = sim[rp][0], sim[rp][1]
            rt, win = accuracy(s1, s2, thre=thre, t0=t0)
            if rt != None:
                winner.append(win)
                time_corr.append(rt)
            else:
                winner.append(np.random.choice([1,2],1)[0])
                time_corr.append(penalty)
                print('Undecided!')
        
        accs.append(winner)
        rtcs.append(time_corr)
    return accs, rtcs
    
def psych_func_complete2(sims, penalty=2.5, thre=0.5, t0=50):
    accs, rtcs, rtws = [], [], []
    for sim in sims:
        time_corr = []; winner = []; time_wrong = []
        for rp in range(len(sim)):
            s1, s2 = sim[rp][0], sim[rp][1]
            rt, win = accuracy(s1, s2, thre=thre, t0=t0)
            if rt != None:
                if win==1: time_corr.append(rt)
                else: time_wrong.append(rt)
                winner.append(win)
            else:
                winner.append(np.random.choice([1,2],1)[0])
                if win==1: time_corr.append(penalty)
                else: time_wrong.append(penalty)
                print('Undecided!')
        
        accs.append(winner)
        rtcs.append(time_corr)
        rtws.append(time_wrong)
    return accs, rtcs, rtws
    
def get_location_dic(traces, thre=0.05):
    lcdic = {}; c = 0
    for trace in traces:
        s1, s2, sg, _, _ = trace
        splus, sminus = (s1+s2)/2, (s1-s2)/2
        for t in range(len(s1)):
            if abs(sminus[t]) > thre:
                loc = splus[t]
                T1, T2, T3, T4, lc1, lc2, lc3 = get_period([trace])
                evaluate = np.any([np.isnan(x) for x in [T1, T2, T3, T4]])
                if not evaluate:
                    t0, cT1, cT2, cT3, cT4 = 50, 50+int(T1*100), 50+int(T1*100)+int(T2*100),\
                        50+int(T1*100)+int(T2*100)+int(T3*100),\
                        50+int(T1*100)+int(T2*100)+int(T3*100)+int(T4*100)
                    try:
                        a2, b2 = curve_fit(affine, splus[cT1:cT2], sg[cT1:cT2])[0]
                        a3, b3 = curve_fit(affine, splus[cT2:cT3], sg[cT2:cT3])[0]
                        a4, b4 = curve_fit(affine, splus[cT3:cT4], sg[cT3:cT4])[0]
                        lcdic[(loc,c)] = [(a2,b2),(a3,b3),(a4,b4)]; c += 1
                    except:
                        pass
                break
    return lcdic
   
def sim_SM_nonoise(ss, params,d=0, c=0):
    gees, gee, geis, gei, gies, gie, giis, gii, Ie, Ii = params
    def rhs(t, v):
        s1, s2, sG1, sG2 = v
        f1 = S_N(s1, gees*s1+gee*s2-gies*sG1-gie*sG2+Ie+pulse(Inp1(c), t))
        f2 = S_N(s2, gees*s2+gee*s1-gies*sG2-gie*sG1+Ie+pulse(Inp2(c), t))
        f3 = S_G(sG1, geis*s1+gei*s2-giis*sG1-gii*sG2+Ii)
        f4 = S_G(sG2, geis*s2+gei*s1-giis*sG2-gii*sG1+Ii)
        return [f1,f2,f3,f4]

    #ss = list(get_EQ_s1(params))[0]
    init = [ss, ss+d] + list(get_ssSG_s1(ss, ss+d, params))
    res = solve_ivp(rhs, (0,3), init, t_eval=np.arange(0,3,0.01))
    return res.y
    
def get_ab_UM(params):
    s1, s2, sg, n1, n2 = sim_UM(params, nsig=0); splus = (s1+s2)/2
    popt, pcov = curve_fit(affine, splus, sg)
    return popt
    
def get_ab_SM(params):
    ss = get_EQ_s1(params)[0]
    s1, s2, sg1, sg2 = sim_SM_nonoise(ss, params)
    splus = (s1+s2)/2; sgplus = (sg1+sg2)/2
    popt, pcov = curve_fit(affine, splus, sgplus)
    return popt
    
def get_abg_SM(params, ss=0.25, d=0.01):
    s1, s2, sg1, sg2 = sim_SM_nonoise(ss, params, d = d)
    sminus = (s1-s2)/2; sgminus = (sg1-sg2)/2
    rt = int(accuracy(s1, s2, t0=0)[0]*100)
    popt, pcov = curve_fit(affine, sminus[:rt], sgminus[:rt])
    return popt
    
def psyched(c, alpha, beta):
    return 0.5 + 0.5*(1-np.exp(-(c/alpha)**beta))
    
def triang(c):
    return -2/51.2**2*c + 2/51.2
    
def Quad(T, a, b, c):
    return [a*t**2+b*t+c for t in T]