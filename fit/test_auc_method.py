import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from scipy.integrate import simps

dist1 = np.random.normal(1,5,size=1000)
dist2 = np.random.normal(-1,5,size=1000)

x = list(dist1)+list(dist2)
y = [1]*len(dist1)+[0]*len(dist2)
x = np.array(x).reshape((len(x),1))
clf = LogisticRegression(solver="liblinear", random_state=0).fit(x, y)
auc_scipy = roc_auc_score(y, clf.predict_proba(x)[:, 1])
fpr, tpr, thresholds = roc_curve(y, clf.predict_proba(x)[:, 1])
print(auc_scipy)

def get_tpr_fpr(dist1, dist2, thre):
    TP = len([x for x in dist1 if x > thre]); FN = len(dist1)-TP
    FP = len([x for x in dist2 if x > thre]); TN = len(dist2)-FP
    tpr = TP/(TP+FN)
    fpr = FP/(FP+TN)
    return tpr, fpr

def auc(dist1, dist2, ndot):
    strt = min(list(dist1)+list(dist2))
    end = max(list(dist1)+list(dist2))
    coors = []
    for thre in np.linspace(strt, end, ndot):
        tpr, fpr = get_tpr_fpr(dist1, dist2, thre)
        coors.append((fpr,tpr))
        
    return coors

coors = auc(dist1, dist2, 20)
xc, yc = zip(*coors)
xc, yc = np.array(xc), np.array(yc)
plt.plot(xc, yc)
s = np.trapz(yc, xc)
print(s)