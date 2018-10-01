# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 11:08:30 2018

@author: wangyf
"""
'''
Objective Oriented Version of Lattice Building functon
'''

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import math
from networkx.algorithms import isomorphism as iso
from itertools import combinations

#%%
'''
Basic functions
'''

def two_points_D(A,B):
    
    '''
    Calculates the distance between two points A and B
    '''
    n = len(A)
    s = 0
    for i in range(n):
        s = s+ (A[i]-B[i])**2
    d = math.sqrt(s)
    d = float(format(d, ' .3f')) # round off to three decimal digits
    return d


def drawing(G):
    '''
    takes in graph g, draw it
    '''
    color = nx.get_node_attributes(G,'color')
    pos=nx.get_node_attributes(G,'pos')
    plt.figure() 
    nx.draw(G, pos, with_labels=True, node_color = list(color.values()))
    
    
def LeaveOneOut(A, a):
    '''
    takes in a list A and returns a new list B by leaving ath element out
    '''     
    B = [x for i,x in enumerate(A) if i!=a]
    
    return B 

def add_z(v, z):
    '''
    takes in a np array and add z coordinates to it
    '''
    vd = np.concatenate((v, np.array([z*np.ones(len(v))]).T), axis =1) 

    return vd
    
#%%
'''
clusters object
'''

class clusters():
    
    def __init__(self, occupancy, NN1, draw):
        
        '''
        takes in the occupancy color vector 
        occupancy[0] is the empty color
        occupancy[1] is the filled color
        '''
        
        self.occupancy = occupancy
        self.empty = self.occupancy[0]
        self.filled = self.occupancy[1]
        self.NN1 = NN1
        self.draw = draw
        
    def gmothers(self, mother, dz):
    
        '''
        takes in mother cooridate list 
        returns connected lattice graph
        '''
        draw_mother = self.draw[0]
        self.mother = mother
        self.nm = len(mother)
        self.dz = dz
        Gm = nx.Graph()
        
        for i in range(self.nm):
            Gm.add_node(i, pos = mother[i][:2], z = str(int(mother[i][2]/self.dz)), color = self.empty)
        
        
        self.edge = []
        self.edge_d = []
        self.edge_z = []
        
        # Add all egdes and calculate the edge distance
        for i in range(self.nm):
            for j in np.arange(i+1,self.nm):
                self.edge.append((i,j))
                self.edge_d.append(two_points_D(mother[i],mother[j]))
                self.edge_z.append(str(int(mother[i][2]/self.dz))+str(int(mother[j][2]/self.dz)))
                
                
        self.ne = len(self.edge)
        for i in range(self.ne): 
            if self.NN1: # only draw 1st Nearest Neighbors 
                if self.edge_d[i] <= 1: 
                    Gm.add_edges_from([self.edge[i]], z = self.edge_z[i], length = self.edge_d[i])
            else:
                Gm.add_edges_from([self.edge[i]],  z = self.edge_z[i], length = self.edge_d[i])
        
        if draw_mother:
            drawing(Gm)
            plt.title('%d lattice points' %self.nm)
            
        return Gm
    
    
    def gconfigurations(self, son):
        
        '''
        takes in mother coordinate list and son's index number and occupancy vector
        returns the shaded son graph
        '''   
        draw_config = self.draw[1]
        ns = len(son)
        Gs = nx.Graph()

        for i in range(self.nm):
            Gs.add_node(i, pos = self.mother[i][:2], z = str(int(self.mother[i][2]/self.dz)), color = self.empty)

        for i in range(self.ne): 
            if self.NN1: # only draw 1st Nearest Neighbors 
                if self.edge_d[i] <= 1:
                    Gs.add_edges_from([self.edge[i]], z = self.edge_z[i], length = self.edge_d[i])
            else:
                Gs.add_edges_from([self.edge[i]], z = self.edge_z[i], length = self.edge_d[i])

        for si in range(ns):
            Gs.node[son[si]]['color'] = self.filled
        
        if draw_config:
            drawing(Gs)
            plt.title('Pd %d' %ns)
        
        return Gs
    

    def gclusters(self, cmother, cson):
    
        '''
        takes in clusters 
        return cluster graph objective
        '''
        draw_clusters = self.draw[2]
        Gc = nx.Graph()
        cns = len(cson)
        
        for i in range(cns):
            c = cson[i]
            Gc.add_node(i, pos = cmother[c][:2], z = str(int(cmother[c][2]/self.dz)), color = self.filled)
            
        cedge = []
        cedge_d = []
        cedge_z = []
        
        for i in range(cns):        
            for j in np.arange(i+1,cns):
                c  = cson[i]
                d = cson[j]
                cedge.append((i,j))
                cedge_d.append(two_points_D(cmother[c],cmother[d])) 
                cedge_z.append(str(int(cmother[c][2]/self.dz))+str(int(cmother[d][2]/self.dz)))
        
        cne = len(cedge)
        for i in range(cne):
           Gc.add_edges_from([cedge[i]], z = cedge_z[i], length = cedge_d[i])
           
        if draw_clusters:            
            drawing(Gc)
            plt.title('Pd %d' %cns)
            
        return Gc
    

    def get_mother(self, mother, dz):
        '''
        takes in mother coordinates list and 
        add mother attribute to the class
        '''
        
        self.Gm  = self.gmothers(mother, dz)
        
        
    def get_configs(self, config):
        
        '''
        takes in configuration index list
        get a list of configurations as graph objects
        '''
        
        self.Gsv = []
        self.nconfig = len(config)
        
        for si in range(self.nconfig):
            son_i = config[si]
            Gs = self.gconfigurations(son_i)
            self.Gsv.append(Gs)
     
        
    def get_clusters(self, cmother, ccluster):
        
        '''
        takes in cluster coordinates list and cluster index list
        returns a list of clusters as graph objects
        '''
        
        self.nc = len(ccluster) # number of clusers
        self.Gcv = [] # list of clusters
        
        for si in range(self.nc):
            cson = ccluster[si] 
            Gc = self.gclusters(cmother,cson)
            self.Gcv.append(Gc)   
           
        

#%%        
class calculations():
    
    '''
    Perform statistical calculation for cluster expansion
    '''
    
    def __init__(self,occupancy):
        
        '''
        takes in the occupancy color vector 
        occupancy[0] is the empty color
        occupancy[1] is the filled color
        '''
        
        self.occupancy = occupancy
        self.empty = self.occupancy[0]
        self.filled = self.occupancy[1]
        
        
    def get_occupancy(self, G, i):
        
        '''
        Get the occupancy from the graph G for node i 
        Occupied is 1 and unoccupied is 0
        '''
        
        if G.nodes[i]['color'] == self.empty: o = -1
        if G.nodes[i]['color'] == self.filled: o = 1 
        
        return o    

    def get_delta_G(self, Gl, Gs):
        
        '''
        takes in larger graph Gl and smaller graph Gs
        find sub isomorphic graphs of Gs from Gl
        calculate the delta value in pi matrix 
        '''
        '''
        if there are more than 2 nodes in a cluster
        '''
        if len(Gs) > 1:            
            '''
            find subgraphs using edge distance match
            '''
            GMl = iso.GraphMatcher(Gl, Gs, edge_match=iso.numerical_edge_match(['length'],[1.0]))
            '''
            find subgraphs using node layer match
            '''
            GMz= iso.GraphMatcher(Gl, Gs, edge_match= iso.categorical_edge_match(['z'],[1.0])  )
            '''
            list down total number of subgraphs niso GMz||GMl
            '''
            x = [y for y in GMz.subgraph_isomorphisms_iter() if y in GMl.subgraph_isomorphisms_iter()]
            
        else:
            '''
            find subgraphs using node layer match
            '''
            GMn = iso.GraphMatcher(Gl, Gs, node_match= iso.categorical_edge_match(['z'],[1]) )
            x = [y for y in GMn.subgraph_isomorphisms_iter()]
            
            
        niso =len(x)
        '''
        save subgraphs to a list
        '''
        subg = list()
        for i in range(niso):    
            subg.append(tuple(x[i].keys()))
        
        '''
        save product into a list 
        and caclulate the sum divid by total number of subgraphs
        '''
        subi = []
        subs = []
        for i in range(niso):
            subi.append([])
            for j in range(len(subg[i])):
                subi[i].append(self.get_occupancy(Gl,subg[i][j]))   
            subs.append(np.product(subi[i]))
        delta = np.sum(subs)/niso
        
        
        return delta
    
    def get_delta_l(self, Gl, Gs):
    
        '''
        takes in larger graph Gl and smaller graph Gs
        find sub isomorphic graphs of Gs from Gl
        calculate the delta value in pi matrix 
        '''
        '''
        if there are more than 2 nodes in a cluster
        '''
        niso =len(Gs)
        
        '''
        save product into a list 
        and caclulate the sum divid by total number of subgraphs
        '''
        subi = []
        subs = []
        for i in range(niso):
            subi.append([])
            for j in range(len(Gs[i])):
                subi[i].append(self.get_occupancy(Gl,Gs[i][j]))   
            subs.append(np.product(subi[i]))
        delta = np.sum(subs)/niso
    
        return delta
    
    def get_pi_matrix(self, G1v, G2v):
        '''
        The function that gets 
            
            configuration graphs, G1v
            cluster graphs, G2v
        and returns the interaction correlation matrix pi
        '''
        n1 = len(G1v)
        n2 = len(G2v)
        pi = np.zeros((n1,n2))
        progress = 0
        
        for i in range(n1):
            for j in range(n2):
                pi[i][j] = self.get_delta_l(G1v[i],G2v[j])
                
                progress = progress + 1
                per = progress/n1/n2 *100
                print('%.2f %% done!' %per)
                        
        self.pi = pi
        
        return pi
    
    def get_J(self, Ev):
        '''
        The function input energy of configurations, Ev
        Returns cluster energy J from linear regression
        '''
        self.Ev = np.array(Ev)
        J = np.linalg.lstsq(self.pi, self.Ev)[0]
        self.J = J
        
        return J
    
    def get_MSE(self):
        '''
        Returns MSE of prediction and real cluster energy
        '''
        ns = len(self.Ev)    
        MSE = np.sum(np.power((np.dot(self.pi,self.J) - self.Ev),2))/ns
        
        self.MSE = MSE
        
        return MSE

#%%
class subgraphs():
    '''
    generate subgraph list with the nodes numbers under the mother graph
    '''
    
    def __init__(self, mother, dz):
        
       self.index= np.arange(len(mother)) # generate the index of nodes
       self.mother = mother
       self.dz = dz
       
    @staticmethod
    def layer_tuple(mother, dz, ci):
        
        '''
        takes in a combo of index and returns tuple of layers they are in 
        '''
    
        n = len(ci)
        index = []
        for i in range(n):
            index.append(ci[i])
            
        layers = []
        
        for i in range(n):
            layers.append(int(mother[index[i]][2]/dz))
            
        layers= tuple(layers)
            
        return layers
    
    @staticmethod   
    def distance_tuple(mother, ci):
        '''
        takes in a combo of index and returns sorted distance between nodes
        '''
        n = len(ci)
        index = []
        
        for i in range(n):
            index.append(ci[i])
        
        combo = list(combinations(index,2))
        ncombo = len(combo) #0 for 1 node, 2 for 2 nodes, 3 for 3 nodes
          
        distances = []
        
        for i in range(ncombo):
            pt1 = mother[combo[i][0]]
            pt2 = mother[combo[i][1]]
            distances.append(two_points_D(pt1, pt2))
            
        distances = tuple(sorted(distances))
        
        return distances
    
    @staticmethod
    def unique_combo(combo, indices_list):
    
        Gv_list = []
        nclusters = len(indices_list)
        
        for i in range(nclusters):
            Gv_list.append([])
            niso = len(indices_list[i])
            for j in range(niso):
                Gv_list[i].append(combo[indices_list[i][j]])
        
        return Gv_list


    def get_s(self, n_atoms):
        
        '''
        Input number of nodes in a subgraph
        Generate combinations among the nodes
        '''
        self.n_atoms = n_atoms
        
        combo = list(combinations(self.index, self.n_atoms))
        ncombo  = len(combo)
        
        
        '''
        generate the inform2tion list
        store the sorted distance of nodes in tuple 1
        + the layer each node is in in tuple 2
        '''
        
        info = [] 
        
        for i in range(ncombo):
            ci  = combo[i]
            
            distances = self.distance_tuple(self.mother, ci)
            layers = self.layer_tuple(self.mother, self.dz, ci)
            
            info.append((distances, layers))
        
        info_set = list(set(info))
        
        index_list =[]
        
        for i in info_set:
            index_list.append(info.index(i))
            
        index_list.sort() # sort the list and take out those indices
        
        s_np = np.array(combo)[index_list]
        
        '''
        convert 2D np array to list
        '''
        
        s_list = []
        for i in range(s_np.shape[0]):
            s_list.append(list(s_np[i]))
            
            
        return s_list
    
    def get_s2(self, n_atoms):
        
        '''
        Input number of nodes in a subgraph
        Generate combinations among the nodes
        '''
        self.n_atoms = n_atoms
        
        combo = list(combinations(self.index, self.n_atoms))
        ncombo  = len(combo)
        
        
        '''
        generate the information list
        store the sorted distance of nodes in tuple 1
        + the layer each node is in in tuple 2
        '''
        
        info = [] 
        
        for i in range(ncombo):
            ci  = combo[i]
            
            distances = self.distance_tuple(self.mother, ci)
            layers = self.layer_tuple(self.mother, self.dz, ci)
            
            info.append((distances, layers))
        
        info_set = list(set(info))
        #print(info_set)
        
        index_list =[]
        indices_list = []
        
        for i in info_set:
            index_list.append(info.index(i))
        
        index_list.sort() # sort the list and take out those indices
            
        for i in index_list:
            indices_list.append([a for a, x in enumerate(info) if x == info[i]])
            
        Gcv_list = self.unique_combo(combo, indices_list)    
           
        return Gcv_list

#%%
class coordination():

    '''
    calculate the coordination number (CN1, CN2) 
    and the general coordination number (GCN)
    '''
    
    def __init__(self,occupancy):
        
        '''
        takes in the occupancy color vector 
        occupancy[0] is the empty color
        occupancy[1] is the filled color
        '''
        
        self.occupancy = occupancy
        self.empty = self.occupancy[0]
        self.filled = self.occupancy[1]
    
    def num_1NN(self, G, i):
    
        '''
        G is a networkx graph
        i is the index of node in G
        
        returns number of nearest neighbors of node i 
        and a list of neighbor index number
        '''
        n_1NN = 0   
        list_1NN = []
        if G.nodes[i]['color'] == self.filled:  # check if the node is occupied 
            for j in list(G.neighbors(i)): #iterate through 1st NN
                if G.nodes[j]['color'] == self.filled: #check if the node is occupied
                    n_1NN = n_1NN + 1
                    list_1NN.append(j)
        else:
            print('No atoms detected at this position') 
        return n_1NN, list_1NN
    
    def num_2NN(self, G,i):
        
        '''
        G is a networkx graph
        i is the index of node in G where CO adsorbs
        
        returns a list of numbers of 2nd nearest neighbors of node i for each 1NN
        and a 2D list of 2NN index numbers
        '''
        n_2NN = []
        list_2NN = []
        if G.nodes[i]['color'] == self.filled: # check if the node is occupied 
            for j in G.neighbors(i):      # iterate through 1st NN
                if G.nodes[j]['color'] == self.filled: # check if the node is occupied 
                    n_2NN.append(self.num_1NN(G,j)[0]) # Add number of 2NN for 1NNs
                    list_2NN.append(self.num_1NN(G,j)[1]) # Add neighbor index number
        else:
            print('No atoms detected at this position')            
        return n_2NN, list_2NN
    
    def cal_CN1(self, G,COsites):
        
        '''
        G is a networkx graph
        COsites are the index of adsorption site
        
        returns 1st CN number
        '''
        
        CN1 = []
        sitetype = len(COsites)
        for i in range(sitetype):
            CN1.append(self.num_1NN(G, COsites[i])[0])
            
        '''
        take arithmaric mean for CN1 for bridge and hollow sites
        '''
        
        CN1 = np.mean(np.array(CN1))   
        
        return CN1
    
    
    def cal_CN2(self, G,COsites):
           
        '''
        G is a networkx graph
        COsites are the index of adsorption site
        
        returns 2nd CN number
        '''
        
        
        list_CN2 = []
        CN2 = []
        
        sitetype = len(COsites)
        
        for i in range(sitetype):
            list_CN2.append(self.num_2NN(G, COsites[i])[0])
            
        '''
        sum up 2NN numbers for each 1NN
        '''
        
        for i in range(sitetype):
            if len(list_CN2[i]) == 0: CN2.append(0)
            else: CN2.append(np.sum(np.array(list_CN2[i])))
        
        '''
        take arithmaric mean for CN2 for bridge and hollow sites
        '''
        
        CN2 = np.mean(np.array(CN2)) 
        
        return CN2
    
    def cal_GCN(self, G, COsites):
        
        '''
        G is a networkx graph
        COsites are the index of adsorption site
        
        returns general coordination number
        '''
        
        GCN = []
        
        sitetype = len(COsites)
        list_1NN = []
        '''
        find all avaiable 1NN index 
        '''
        for i in range(sitetype):
            list_1NN = list_1NN + self.num_1NN(G, COsites[i])[1]
        
        '''
        Use set to avoid double counting
        '''    
        list_1NN = list(set(list_1NN))
        '''
        Get CN for these 1NN nodes
        '''
        for i in list_1NN:
            GCN.append(self.num_1NN(G,i)[0])
    
        '''
        Set weight based on Pd(111)
        '''
        
        if len(COsites) == 1: weight = 12
        if len(COsites) == 2: weight = 18
        if len(COsites) == 3: weight = 22
        
        GCN = np.sum(np.array(GCN))/weight
    
        return GCN
    
    
    def num_Ce1NN(self, G, i):
        
        '''
        G is a networkx graph
        i is the index of node in G
        
        returns the flag of whether the atom is next to a Ce atom
        '''
        
        nCe = 0 
        if G.nodes[i]['color'] == self.filled: # check if the node is occupied
            if G.nodes[i]['z'] == '1':  # check if the node is in base layer
                nCe = 1 # 1 means the atom is in contact with 3 Ce atoms underneath
        return nCe
            
    def num_Ce2NN(self, G, i):
        
        '''
        G is a networkx graph
        i is the index of node in G
        
        returns the number of 2nd nearest Ce neighbors of node i 
        and a list of atom adjacent to Ce base layer
        '''
    
        n_2NN = 0   
        list_2NN = []
        if G.nodes[i]['color'] == self.filled: # check if the node is occupied
            for j in list(G.neighbors(i)): # iterate through 1st NN
                n_2NN = n_2NN + self.num_Ce1NN(G,j)  # check if the node is next to Ce
                if self.num_Ce1NN(G,j): list_2NN.append(j) # Append the index to a list
        else:
            print('No atoms detected at this position') 
        return n_2NN, list_2NN
    
    def cal_CeCN1(self, G, COsites):
        
        '''
        G is a networkx graph
        COsites are the index of adsorption site
        
        returns coordination number of Ce
        '''
        
        CN1 = []
        sitetype = len(COsites)
        
        for i in range(sitetype):
            CN1.append(self.num_Ce1NN(G, COsites[i]))
        
        '''
        take arithmaric mean for CN2 for bridge and hollow sites
        '''
        
        CN1 = np.mean(np.array(CN1)) * 3 # each Pd atom is coordinate by 3 Ce
        
        return CN1
        
    def cal_CeCN2(self, G, COsites):
        
        '''
        G is a networkx graph
        COsites are the index of adsorption site
        
        returns 2nd coordination number of Ce
        '''
        
        CN2 = []
        sitetype = len(COsites)
        
        for i in range(sitetype):
            CN2.append(self.num_Ce2NN(G, COsites[i])[0])
        
        '''
        take arithmaric mean for CN2 for bridge and hollow sites
        '''
        
        CN2 = np.mean(np.array(CN2)) *3 
           
        return CN2
    
    def cal_CeGCN(self, G, COsites):
    
        '''
        G is a networkx graph
        COsites are the index of adsorption site
        
        returns general coordination number of Ce
        '''
        
        GCN = []
        
        sitetype = len(COsites)
        list_1NN = []
        
        '''
        find all avaiable 1NN next to Ce index 
        '''
        
        for i in range(sitetype):
            list_1NN = list_1NN + self.num_Ce2NN(G, COsites[i])[1]
            
        list_1NN = list(set(list_1NN))
        
        '''
        Check if Ce is around for 1NN nodes
        '''
        
        for i in list_1NN:
            GCN.append(self.num_Ce1NN(G,i))
        
        
        if len(COsites) == 1: weight = 3
        if len(COsites) == 2: weight = 5
        if len(COsites) == 3: weight = 6
        
        GCN = np.sum(np.array(GCN))/weight * 3
    
        return GCN
    
    def get_CNs(self, G, COsites):
        
        '''
        Take in configuration G
        CO adsorption configuration index list 
        and CO sites index list
        add properties to the self object
        '''
        self.CN1 = self.cal_CN1(G,COsites)
        self.CN2 = self.cal_CN2(G,COsites)
        self.GCN = self.cal_GCN(G,COsites)

        self.CeCN1 = self.cal_CeCN1(G,COsites)
        self.CeCN2 = self.cal_CeCN2(G,COsites)
        self.CeGCN = self.cal_CeGCN(G,COsites)
   
    def get_z(self, G, COsites):
        
        '''
        Take in configuration G
        CO adsorption configuration index list 
        and CO sites index list
        add average layer number to the self object
        '''
        list_z = []
        sitetype = len(COsites)
        
        for i in range(sitetype):
            list_z.append(int(G.nodes[COsites[i]]['z']))
            
            
        self.z = np.mean(np.array(list_z))    
            
            
            
            
            
            
            
            
            
            