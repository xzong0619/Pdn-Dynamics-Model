# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 15:54:14 2018

@author: wangyf
"""

GG = nx.Graph()
GG.add_nodes_from([1,2,3])

nx.draw(GG, node_color = ['green','red','blue'])