'''
Plot the standard deviation of the psychometric function.
'''

import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import aadm.basics as bcs

'''=====Define parameters here====='''
gmavals = [[0.3,0.6,1.7,1,1], [0.3,0.6,0.6,1,1],\
		   [0.3,1.7,0.6,1,1], [0.3,0.8,0.8,1,1],\
		   [0.3,1.25,0.8,1,1], \
		   [0.3,1,1,0,2], [0.1,1,1,0,2],\
		   [0.5,1,1,0,2], [0.75,1,1,0,2], [0.83,1,1,0,2],\
		   [0.9,1,1,0,2], [0.5,0.6,0.6,1,1], [0.75,0.8,0.8,1,1],\
           [0.5,0.6,1.7,1,1], [0.75,0.6,1.7,1,1]]
'''================================'''

infofile = os.path.join(bcs.datapath(), 'files', 'fullacc_4c.p')
accdic = pickle.load(open(infofile, 'rb'))

for k in range(len(gmavals)):
    data = np.asarray(accdic[tuple(gmavals[k])])
    psych = [np.mean(2-np.array(row)) for row in data]
    std = [np.std(2-np.array(row))/np.sqrt(len(row)-1) for row in data]
    plt.errorbar(range(6), psych, yerr=std)