# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 13:58:39 2018

@author: wangyf
"""

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
Plot 1 Cluster Size vs DFT cluster energy
'''


fig, ax = plt.subplots(figsize= (6,5))
ax.scatter(NPd_list, Ec,s=120, facecolors='none', edgecolors='b')
ax.set_xlabel('Cluster Size N')
ax.set_ylabel('DFT Cluster Energy (eV)')
plt.axis([1, 25, -25, 1])
plt.xticks(np.arange(0, 30, 5))
plt.yticks(np.arange(-25, 5, 5))
plt.show()

'''
Plot 2 Parity Plot
'''
y_predict_all = lasso_cv.predict(X)

    
fig, ax = plt.subplots(figsize= (6,5))
ax.scatter(y_predict_all, Ec, s=120, facecolors='none', edgecolors='r')
plt.xlabel("Predicted Cluster Energy (eV)")
#plt.ylabel("DFT Cluster Energy (eV)")

lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
]

# now plot both limits against eachother
ax.plot(lims, lims, 'k--', alpha=0.75, zorder=0, lw=2)
plt.axis([-25,1, -25, 1])
plt.xticks(np.arange(-25, 5, 5))
plt.yticks(np.arange(-25, 5, 5))
plt.show()