# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 17:59:55 2018

@author: wangyf
"""
import matplotlib as mat
import matplotlib.pyplot as plt
import numpy as np


# up to 24 colors stored
colors_pool = [ 'b','yellowgreen', 'gold', 'lightcoral',
               'lightskyblue', 'darkgreen', 'orange','salmon',
               'powderblue','olivedrab', 'burlywood',  'indianred', 
               'steelblue', 'lawngreen', 'y', 'hotpink',
               'slategrey', 'palegreen', 'sandybrown', 'tomato',
               'darkviolet', 'lightgreen', 'tan','maroon']
        
        

def PlotOptions():

    '''
    Set matplotlib options
    '''

    mat.rcParams['mathtext.default'] = 'regular'
    mat.rcParams['text.latex.unicode'] = 'False'
    mat.rcParams['legend.numpoints'] = 1
    mat.rcParams['lines.linewidth'] = 2
    mat.rcParams['lines.markersize'] = 12
            

def PlotTimeSeries(x_series, y_series, xlab = 'Time (s)', ylab = '', xlimit = [], series_labels = [], fname = '', logscale = False):
    
    '''
    Plot multiple series against time
    '''
    
    PlotOptions()
    plt.figure(figsize=(10,8))
    
    for i in range (len(y_series)):
        plt.plot(x_series, y_series[i], color = colors_pool[i], linewidth = 3.0)
    
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.xlabel(xlab, size=24)
    plt.ylabel(ylab, size=24)
	
    if not xlimit == []:
        plt.xlim(xlimit)
	
    if not series_labels == []:
        plt.legend(series_labels, loc=1, prop={'size':20}, frameon=False)
    plt.tight_layout()
    
    if logscale:
        plt.yscale('log')
    
    if fname == '':
        plt.show()
    else:
        plt.savefig(fname)
        plt.close()
        
def PlotPie(x_series, series_labels = [], fname = '') :
    
    '''
    Plot pie chart showing the end state population distribution
    '''
    PlotOptions()
    plt.figure(figsize=(8,8))
    mat.rcParams['font.size'] = 20
    n_s = len(x_series)
    colors = colors_pool[:n_s]
    plt.axis('equal')
    
    explode = []
    for i in range(n_s): 
        
        if x_series[i] == 0:
            del colors[i]
            del series_labels[i]
            nx_series = np.delete(x_series, i)
            
    for i in range(len(nx_series)):     
            if nx_series[i]/sum(nx_series) >= 0.1:
                explode.append(0.0)
            else:
                explode.append(0.1)
                
    plt.pie(nx_series, explode = explode, labels = series_labels, colors = colors, autopct = '%1.1f%%', shadow = True, startangle = 140)
    
    if fname == '':
        plt.show()
    else:
        plt.savefig(fname)
        plt.close()
