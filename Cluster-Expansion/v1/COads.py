# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 15:01:12 2018

@author: wangyf
"""
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib


pi =np.matrix([[ 1, 0, 0,  0, 0, 0,  0, 0, 0],  #1
                [ 2, 0, 0,  0, 0, 0,  0, 0, 0], #2
                [ 2, 1, 0,  0, 0, 0,  0, 0, 0], #3
                [ 2, 2, 0,  0, 0, 0,  0, 0, 0], #4-2d
                [ 3, 0, 0,  1, 0, 0,  0, 0, 0], #4-3d
                [ 2, 2, 1,  0, 0, 0,  0, 0, 0], #5-2d
                [ 3, 0, 0,  1, 1, 0,  0, 0, 0], #5-3d
                [ 3, 0, 0,  1, 2, 0,  0, 0, 0], #6
                [ 3, 1, 1,  1, 1, 0,  0, 0, 0], #7
                [ 2, 1, 1,  1, 1, 1,  0, 0, 0], #8
                [ 3, 2, 1,  1, 1, 1,  1, 0, 0], #9
                [ 3, 0, 0,  3, 2, 0,  2, 0, 0], #10
                [ 3, 1, 2,  2, 1, 2,  1, 2, 0], #11
                [ 3, 1, 1,  1, 2, 1,  1, 1, 0], #12
                [ 3, 1, 1,  1, 2, 1,  1, 1, 0], #13
                [ 3, 1, 1,  1, 2, 1,  1, 1, 0], #14
                [ 3, 1, 1,  1, 2, 1,  1, 2, 0], #15
                [ 3, 1, 1,  1, 2, 1,  1, 2, 0], #16
                [ 3, 1, 1,  1, 2, 1,  1, 2, 0], #17
                [ 3, 1, 1,  1, 2, 1,  1, 2, 0], #18
                [ 3, 1, 2,  1, 2, 1,  1, 2, 0], #19
                [ 3, 1, 2,  1, 2, 1,  1, 2, 0], #20
                [ 2, 2, 3,  1, 5, 1,  1, 2, 0], #20
                [ 1, 0, 0,  3, 0, 0,  3, 3, 0]]) #20
ref = -2.4 #binding energy of CO onto Pd1
Ec_DFT =  np.array([-2.4, -2.98, -3.06, -2.94, -2.58, -3.03, -2.23, -2.28, -2.20, -2.36, 
                    -2.03, -2.31, -2.29, -2.43, -2.38, -2.31, -2.25, -2.29, -2.33, -2.33,
                    -2.08, -2.00, -1.91, -1.43])
Ec = Ec_DFT - ref
ns = len(Ec)
J = np.linalg.lstsq(pi, Ec)[0]
MSE = np.sum(np.power((np.dot(pi,J) - Ec),2))/ns 

Ec_model = np.array(np.dot(pi,J))[0] + ref

'''
Parity plot of regression
'''

plt.figure(figsize=(20,20))

fig, ax = plt.subplots()
ax.scatter(Ec_DFT,Ec_model, s=60, facecolors='none', edgecolors='r')
ax.axis([-3.2, -1.2, -3.2, -1.2])
plt.xlabel("Eads DFT (eV)")
plt.ylabel("Eads Model Prediction(eV)")


lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
]

# now plot both limits against eachother
ax.plot(lims, lims, 'k--', alpha=0.75, zorder=0)
ax.set_xlim(lims)
ax.set_ylim(lims)

font = {'family' : 'normal',
        'size'   : 15}

matplotlib.rc('font', **font)
plt.show()


