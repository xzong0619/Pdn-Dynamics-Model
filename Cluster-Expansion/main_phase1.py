# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 10:44:28 2018

@author: wangyf
"""

#%%
'''
main
'''
from Lattice_fun import *  

mother = [(0,0), (1,0), (1/2, 3**0.5/2), (3/2, 3**0.5/2), (0, 3**0.5),
          (1/2, -3**0.5/2), (3/2, -3**0.5/2), (0, -3**0.5),
          (-1,0), (-1/2, -3**0.5/2), (-3/2, -3**0.5/2),
          (-1/2, 3**0.5/2), (-3/2, 3**0.5/2)]
config = [[0],
          [0,1],
          [0,1,2],
          [0,1,2,3],
          [0,1,2,3,5],
          [0,1,2,3,5,6],
          [0,1,2,3,5,9],
          [0,1,2,3,5,11],
          [0,1,2,5,9,8,11],
          [0,1,2,3,5,9,7],
          [0,1,2,5,8,11,12],
          [0,1,2,8,9,11,12]]
'''
Energy of clusters in eV
'''

Ec = [0, -0.45, -1.5, -2.56, -3.41, -4.49, -4.38, -4.32, -6.07, -4.9, -5.11, -5.12]


     
empty = 'grey'
filled = 'r'
occ = [empty, filled]

Gm = gmothers(mother, occ)
    
        
'''
Creat 12 configurations
'''
ns = len(config)  # number of configurations 
Gsv = [] # list of configurations 
for si in range(ns):
    son = config[si]
    Gs = gconfigurations(mother,son, occ)
    Gsv.append(Gs)


#%%
'''
Creat 6 clusters
'''

cmother = [(0,0), (1,0), (1/2, 3**0.5/2), (3/2, 3**0.5/2), (0, 3**0.5)]
cconfig = [[0], [0,1], [0,3], [1,4],[0,1,2],[1,2,4], [0,1,3]]


nc = len(cconfig) # number of clusers
Gcv = [] # list of clusters
for si in range(nc):
    cson = cconfig[si] 
    Gc = gclusters(cmother,cson, occ)
    Gcv.append(Gc)       
        
#%% Stattistical analysis
'''
creat pi matrix
size of number of configuration * numbers of clusters
'''

J, pi =  get_J(Ec,Gsv,Gcv, occ)      
MSE = np.sum(np.power((np.dot(pi,J) - Ec),2))/ns 



#%%
'''
Calculate Leave One Out CV score
'''
Ec_predict = np.zeros(ns)


for i in range(ns):
    '''
    i is the index to be removed from the list
    '''
    Ec_LOO =  LeaveOneOut(Ec,i)
    Gsv_LOO = LeaveOneOut(Gsv,i)
    
    J_LOO = get_J(Ec_LOO, Gsv_LOO, Gcv, occ)[0]
    Ec_predict[i] = np.dot(pi[i],J_LOO)
    
CV = np.sum(np.power((Ec_predict - Ec),2))/ns    

    