# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 14:58:29 2018

@author: wangyf
"""

'''
plot PCA
'''
import matplotlib.pyplot as plt 
import matplotlib

font = {'size'   : 20}

matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 1.5
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['xtick.major.width'] = 2
matplotlib.rcParams['ytick.major.size'] = 8
matplotlib.rcParams['ytick.major.width'] = 2
'''
Plot 1 initial Data distribution
'''
cnt = 0
fig, ax = plt.subplots(figsize= (6,5))
for site, col in zip(('top', 'bridge', 'hollow'),
            ('red', 'green', 'blue')):
    indices = np.where(np.array(sitetype_list) == site)[0]
    plt.scatter(X[:,cnt][indices],
                y[indices],
                label=site,
                facecolor = col, 
                alpha = 0.5,
                s  = 100)

ax.set_xlabel('Cluster Size N')
ax.set_ylabel('DFT CO  Binding Energy (eV)')  
plt.xticks(np.arange(0, 30, 5))
plt.yticks(np.arange(-3.0, -1.0, 0.5))
#plt.legend(bbox_to_anchor = (1.02, 1),loc= 'upper left', frameon=False)
plt.show() 


'''
Plot 2 Parity Plot
'''
yobj = y
ypred = y_pca
fig, ax = plt.subplots(figsize= (6,5))
for site, col in zip(('top', 'bridge', 'hollow'),
                ('red', 'green', 'blue')):
        indices = np.where(np.array(sitetype_list) == site)[0]
        ax.scatter(ypred[indices],
                    yobj[indices],
                    label=site,
                    facecolor = col, 
                    alpha = 0.5,
                    s  = 100)
ax.plot([yobj.min(), yobj.max()], [yobj.min(), yobj.max()], 'k--', lw=2)
ax.set_xlabel('Predicted CO Binding Energy (eV)')
plt.xticks(np.arange(-3.0, -1.0, 0.5))
plt.yticks(np.arange(-3.0, -1.0, 0.5))
plt.legend(loc= 'lower right', frameon=False)
plt.show()