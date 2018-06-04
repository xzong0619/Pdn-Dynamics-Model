# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 14:07:25 2018

@author: wangyf
"""

import Cluster_structure as cs

n = 3
old_edge = [(1, 2), (1, 3) ,(2, 3)]
old_wt_i = []

cs.plot_graph([old_edge])


n = 4

old_edge = [(1, 2), (1, 3) ,(2, 3)]
old_wt_i = []
Pd4 = cs.Pdn(n, old_edge, old_wt_i)
Pd4.new_graph()
ni_u = Pd4.ni_u
new_H = Pd4.H
new_wt_i = Pd4.wt_i
print('Pd%d done' %n)
cs.plot_graph(Pd4.H)

n = 5
Pd5 = cs.Parallel_growth(n, ni_u, new_H, new_wt_i)
Pd5.unique()
ni_u = Pd5.mi_u
new_H = Pd5.H_u
new_wt_i = Pd5.wt_i_u
print('Pd%d done' %n)
cs.plot_graph(Pd5.H_u)

n = 6
Pd6 = cs.Parallel_growth(n, ni_u, new_H, new_wt_i)
Pd6.unique()
ni_u = Pd6.mi_u
new_H = Pd6.H_u
new_wt_i = Pd6.wt_i_u
print('Pd%d done' %n)
cs.plot_graph(Pd6.H_u)

n = 7
Pd7 = cs.Parallel_growth(n, ni_u, new_H, new_wt_i)
Pd7.unique()
ni_u = Pd7.mi_u
new_H = Pd7.H_u
new_wt_i = Pd7.wt_i_u
print('Pd%d done' %n)
cs.plot_graph(Pd7.H_u)

n = 8
Pd8 = cs.Parallel_growth(n, ni_u, new_H, new_wt_i)
Pd8.unique()
ni_u = Pd8.mi_u
new_H = Pd8.H_u
new_wt_i = Pd8.wt_i_u
print('Pd%d done' %n)

n = 9
Pd9 = cs.Parallel_growth(n, ni_u, new_H, new_wt_i)
Pd9.unique()
ni_u = Pd9.mi_u
new_H = Pd9.H_u
new_wt_i = Pd9.wt_i_u
print('Pd%d done' %n)

'''
n = 10
Pd10 = cs.Parallel_growth(n, ni_u, new_H, new_wt_i)
Pd10.unique()
ni_u = Pd10.mi_u
new_H = Pd10.H_u
new_wt_i = Pd10.wt_i_u
print('Pd%d done' %n)

n = 11
Pd11 = cs.Parallel_growth(n, ni_u, new_H, new_wt_i)
Pd11.unique()
ni_u = Pd11.mi_u
new_H = Pd11.H_u
new_wt_i = Pd11.wt_i_u
print('Pd%d done' %n)
'''