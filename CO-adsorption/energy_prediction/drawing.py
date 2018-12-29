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

 # Loop through all sites
def find_neighbor_list(sites_pos):
    neighbor_list = []
    n_sites = len(sites_pos)
    for site_1 in range(n_sites):
        for site_2 in range(n_sites):
        
            c1 = sites_pos[site_1,:]
            c2 = sites_pos[site_2,:]
            d = np.linalg.norm( c1 - c2 )
    
            if np.logical_and(d>=NN1[0], d<=NN1[1]):
                neighbor_list.append([site_1, site_2])
    return neighbor_list

def plot_neighbors(sites_pos, neighbor_list, c):

    for pair in neighbor_list: # neighbors
        p1 = sites_pos[pair[0],:]
        p2 = sites_pos[pair[1],:]
        plt.plot([p1[0], p2[0]], [p1[1], p2[1]], '-', color = c, linewidth=1.5,  alpha=0.5)


'''
Plot the lattice first
'''

unit = 0.27 #In nm, 1 angstrom = 0.1 nm 
NN1 = unit * np.array([0.9, 1.1])

pts = unit *  np.array([(0,0), (1,0), (2,0), (3,0),
                        (1/2, 3**0.5/2), (1/2+1,3**0.5/2), (1/2+2,3**0.5/2),
                        (1,3**0.5), (2,3**0.5), (1.5, 3**0.5/2*3)])
    
n_pts = len(pts)        
pts_3D = np.column_stack((pts, np.zeros(n_pts)))
atoms_plane = Atoms('H'+str(len(pts)), positions = pts_3D)
#view(atoms_plane)

pts_index = list(range(0, n_pts))
t_index = []
for i in pts_index: t_index.append([i])
n_top = len(t_index)

bridge_sites = []     
pairs  = list(combinations(pts_index,2))
b_index = find_bridge_pairs(pairs, atoms_plane)
n_bridge = len(b_index)
        
hollow_sites = [] 
triples  = list(combinations(pts_index,3))
h_index = find_hollow_triples(triples, atoms_plane)
n_hollow = len(h_index)    
   
sites_index = t_index + b_index + h_index
n_sites = len(sites_index)

sites_pos = np.zeros((n_sites, 2))
    
for i,site in enumerate(sites_index):
    pos = []
    for j in site: pos.append(atoms_plane[j].position)
    sites_pos[i,:] = np.mean(pos, axis = 0)[:2] 

t_pos =  sites_pos[:n_top]
b_pos = sites_pos[n_top: n_top + n_bridge]
h_pos = sites_pos[n_top + n_bridge:]

all_NN = find_neighbor_list(sites_pos)
t_NN = find_neighbor_list(t_pos)
b_NN = find_neighbor_list(b_pos)
h_NN = find_neighbor_list(h_pos)

plt.figure(figsize=(5, 4))                    
for sitetype, col, mk in zip(('top', 'bridge', 'hollow'),
                        ('red', 'green', 'blue'), 
                        ('o','s' ,'v')):
    indices = np.where(np.array(sitetype_list) == sitetype)[0]
    x_pos = sites_pos[indices][:,0]
    y_pos = sites_pos[indices][:,1]
    plt.scatter(x_pos, y_pos, s= 120,  label=sitetype, marker = mk,  color = col, alpha=1)

#plot_neighbors(sites_pos, all_NN, 'grey')
plot_neighbors(t_pos, t_NN, 'r')
plot_neighbors(b_pos, b_NN, 'green')
plot_neighbors(h_pos, h_NN, 'blue')

x_lim = (sites_pos[:,0].min() - unit/5, sites_pos[:,0].max() + unit/5)   
y_lim = (sites_pos[:,1].min() - unit/5, sites_pos[:,1].max() + unit/5)     
plt.xlabel('x (nm)')
plt.ylabel('y (nm)')
plt.xlim(x_lim)
plt.ylim(y_lim)
plt.legend(bbox_to_anchor = (1.02, 1), loc= 'upper left', frameon=False)
plt.show()

#%%
'''
Plot the intensity of CO adsorption energy
'''
plt.figure(figsize=(5, 4))   
x_pos = sites_pos[:,0]
y_pos = sites_pos[:,1] 
#plot_neighbors(sites_pos, all_NN, 'grey')
plt.scatter(x_pos, y_pos, s= 200, c = y_pcg, edgecolors = 'k', cmap = 'viridis')

x_lim = (sites_pos[:,0].min() - unit/5, sites_pos[:,0].max() + unit/5)   
y_lim = (sites_pos[:,1].min() - unit/5, sites_pos[:,1].max() + unit/5)     
plt.xlabel('x (nm)')
plt.ylabel('y (nm)')
plt.xlim(x_lim)
plt.ylim(y_lim)
#plt.colorbar()
#plt.legend(bbox_to_anchor = (1.02, 1),loc= 'upper left', frameon=False)
plt.show()
