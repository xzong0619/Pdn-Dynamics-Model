# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 15:56:13 2018

@author: wangyf
"""

from structure_constants import *
import lattice_functions as lf
import pickle
import pandas as pd
import numpy as np

empty = 'grey'
filled = 'r'
occ = [empty, filled]  
    
[Gm, Gsv, Gcv1, Gcv2, Gcv3] = pickle.load(open("clusters_NN1.p", "rb"))


nCOs = len(COindex)
cn = lf.coordination(occ)

nCN1 = []
nCN2 = []
nGCN = []

nCeCN1 = []
nCeCN2 = []
nCeGCN = []

nZ  = []
for i in range(nCOs):
    
    cn.get_CNs(Gsv[COindex[i]],COsites[i])
    cn.get_z(Gsv[COindex[i]],COsites[i])
    
    nCN1.append(cn.CN1)
    nCN2.append(cn.CN2)
    nGCN.append(cn.GCN)
    
    nCeCN1.append(cn.CeCN1)
    nCeCN2.append(cn.CeCN2)
    nCeGCN.append(cn.CeGCN)
    
    nZ.append(cn.z)

df = pd.DataFrame({"CN1": nCN1,
                   "CN2": nCN2,
                   "GCN": nGCN,
                   "CeCN1": nCeCN1,
                   "CeCN2": nCeCN2,
                   "CeGCN": nCeGCN,
                   "z": nZ})
'''
Use clipboard methd
'''

#df = pd.read_clipboard()
X = df.iloc[:,:].values
np.save('Xmatrix', X, allow_pickle = True)
