# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 22:24:45 2018

@author: yifan
"""
'''
Make the plot predicted CO binding energy
'''

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt 
import matplotlib
from ase import Atoms

matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['axes.linewidth'] = 1.5
matplotlib.rcParams['xtick.major.size'] = 8
matplotlib.rcParams['xtick.major.width'] = 2
matplotlib.rcParams['ytick.major.size'] = 8
matplotlib.rcParams['ytick.major.width'] = 2


def find_bridge_pairs(Pd_pairs, atoms):
    
    bridge_pairs = []
    for pair in Pd_pairs:
        Pd_Pd = atoms.get_distances([pair[0]], [pair[1]])
        if np.logical_and(Pd_Pd>=NN1[0], Pd_Pd<=NN1[1]):
            bridge_pairs.append(list(pair))
    
    return bridge_pairs
            
def find_hollow_triples(Pd_triples , atoms):
    
    hollow_triples = []
    for triple in Pd_triples:
        Pd_Pd1 = atoms.get_distances(triple[0], [triple[1], triple[2]])
        Pd_Pd2 = atoms.get_distances([triple[1]], [triple[2]])
        flag1 = np.logical_and(Pd_Pd1>=NN1[0], Pd_Pd1<=NN1[1])
        flag2 = np.logical_and(Pd_Pd2>=NN1[0], Pd_Pd2<=NN1[1])
      
        if np.all(list(flag1)+list(flag2)):
            hollow_triples.append(list(triple))
    
    return hollow_triples


'''
Plot the lattice first
'''

unit = 1
NN1 = unit*np.array([0.9, 1.1])

pts = np.array([(0,0), (1,0), (2,0), (3,0),
                (1/2, 3**0.5/2), (1/2+1,3**0.5/2), (1/2+2,3**0.5/2),
                (1,3**0.5), (2,3**0.5), (1.5, 3**0.5/2*3)])
        
pts_3D = np.column_stack((pts, np.zeros(len(pts))))
atoms_plane = Atoms('H'+str(len(pts)), positions = pts_3D)
#view(atoms_plane)

pts_index = list(range(0, len(pts)))
t_index = []
for i in pts_index: t_index.append([i])

bridge_sites = []     
pairs  = list(combinations(pts_index,2))
b_index = find_bridge_pairs(pairs, atoms_plane)

        
hollow_sites = [] 
triples  = list(combinations(pts_index,3))
h_index = find_hollow_triples(triples, atoms_plane)
    
   
sites_index = t_index + b_index + h_index

sites_pos = np.zeros((len(sites_index), 2))
    
for i,site in enumerate(sites_index):
    pos = []
    for j in site: pos.append(atoms_plane[j].position)
    print(pos)
    sites_pos[i,:] = np.mean(pos, axis = 0)[:2] 
        
        
for sitetype, col in zip(('top', 'bridge', 'hollow'),
                        ('red', 'green', 'blue')):
    indices = np.where(np.array(sitetype_list) == sitetype)[0]
    plt.scatter(y_pcg[:,cnt][indices],
             range= (1.5, 4),
             label=site,
             color = col,
             alpha=0.5)
plt.xlabel()
plt.ylabel('Count')
plt.xlim((1.5,4))
plt.ylim((0,10))
plt.legend(bbox_to_anchor = (1.02, 1),loc= 'upper left', frameon=False)
plt.tight_layout()
plt.show()
