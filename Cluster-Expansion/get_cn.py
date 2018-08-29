# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 15:56:13 2018

@author: wangyf
"""

from structure_constants import *
import lattice_functions as lf
import pickle

empty = 'grey'
filled = 'r'
occ = [empty, filled]  
    
[Gm, Gsv, Gcv1, Gcv2, Gcv3] = pickle.load(open("clusters_NN1.p", "rb"))


'''
G = Gsv[25] #No.26 Pd10a [0,1,2,6,11,10,14,29,28,27]
COsites =  [2, 14, 27]

cn = lf.coordination(occ)
cn.get_CNs(G,COsites)

CN1 = cn.CN1
CN2 = cn.CN2
GCN = cn.GCN

CeCN1 = cn.CeCN1
CeCN2 = cn.CeCN2
CeGCN = cn.CeGCN
'''

'''
Iterate and get all CNs for 24 CO adsorbed structures
'''
nCOs = len(COindex)
cn = lf.coordination(occ)

nCN1 = []
nCN2 = []
nGCN = []

nCeCN1 = []
nCeCN2 = []
nCeGCN = []

for i in range(nCOs):
    
    cn.get_CNs(Gsv[COindex[i]],COsites[i])
    nCN1.append(cn.CN1)
    nCN2.append(cn.CN2)
    nGCN.append(cn.GCN)
    
    nCeCN1.append(cn.CeCN1)
    nCeCN2.append(cn.CeCN2)
    nCeGCN.append(cn.CeGCN)
    
    
    