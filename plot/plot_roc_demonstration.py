'''
Demonstrative plot of response distribution and roc.
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

x = np.linspace(-5,5)
fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
plt.plot(x, norm.pdf(x, loc=1, scale=1), 'g')
plt.plot(x, norm.pdf(x, loc=-1, scale=1), 'r')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('Response'); plt.ylabel('pdf')
plt.tight_layout()
plt.savefig('response_distribution.png', dpi=300); plt.clf()

fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
x = np.arange(0,1,0.1)
y = -np.square(x-1)+1 
plt.plot(x, y, 'k')
plt.plot(x, y+np.random.normal(0,0.05,size=len(x)), 'b.')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel('False positive rate'); plt.ylabel('True positive rate')
plt.tight_layout()
plt.savefig('roc.png', dpi=300); plt.clf()