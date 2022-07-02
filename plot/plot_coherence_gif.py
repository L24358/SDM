'''
Plot gif that illustrates the 2AFC experiment considered, and the concept of coherence.
'''
import os
import imageio
import numpy as np
import matplotlib.pyplot as plt
import dynalysis.classes as clss

N = 150
b = clss.branch('gif', os.getcwd()); b.mkdir()

for p in [0.5, 0.7, 0.9]:

    right_centers = np.random.uniform(low=-5, high=5, size=(int(N*p),2))
    left_centers = np.random.uniform(low=-5, high=5, size=(int(N*(1-p)),2))
    
    fig, ax = plt.subplots(figsize=(10,10))
    fig.patch.set_facecolor('black')
    plt.xlim((-5,5))
    plt.ylim((-5,5))
    plt.axis('off')
    ax.text(-2.5, -0.5, r"$\bf{p = " + '{:.1f}'.format(p) + "}$", color='white', fontsize=80)
    plt.savefig(os.path.join(b.pathlink, 'p='+str(p)+'.png'), dpi=50, facecolor='k', transparent=False)
    
    for t in range(30):
        right_centers = np.transpose([right_centers[:,0]+0.1+np.random.normal(0,0.05,size=int(N*p)),\
                                      right_centers[:,1]]+np.random.normal(0,0.05,size=int(N*p)))
        left_centers = np.transpose([left_centers[:,0]-0.1+np.random.normal(0,0.05,size=int(N*(1-p))),\
                                      left_centers[:,1]]+np.random.normal(0,0.05,size=int(N*(1-p))))
        
        fig, ax = plt.subplots(figsize=(10,10))
        fig.patch.set_facecolor('black')
        
        for center in left_centers:
            circle = plt.Circle(center, 0.1, color='white')
            ax.add_patch(circle)
            
        for center in right_centers:
            circle = plt.Circle(center, 0.1, color='white')
            ax.add_patch(circle)
        
        plt.xlim((-5,5))
        plt.ylim((-5,5))
        plt.axis('off')
        ax.text(1.5, 4, r"$\bf{p = " + '{:.1f}'.format(p) + "}$", color='white', fontsize=40)
        plt.savefig(os.path.join(b.pathlink, 'p='+str(p)+'_t='+str(t)+'.png'), dpi=50, facecolor='k', transparent=False)
        plt.close('all')
 
images = []       
for p in [0.5, 0.7, 0.9]:
    for rp in range(10):
        images.append(imageio.imread(os.path.join(b.pathlink, 'p='+str(p)+'.png')))
    for t in range(30):
        images.append(imageio.imread(os.path.join(b.pathlink, 'p='+str(p)+'_t='+str(t)+'.png')))
imageio.mimsave('coherence.gif', images, duration=0.1)
            
    