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
        plt.plot(x_series, y_series[i], #marker = 'o', markerfacecolor = None, 
                 color = colors_pool[i], linewidth = 3.0)
    
    #plt.xticks(size=20)
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
    mat.rcParams['font.size'] = 12
    n_s = len(x_series)
    colors = colors_pool[:n_s]
    plt.axis('equal')
    
    explode = []
    ncolors = []
    nseries_labels = []
    nx_series = []
    
    for i in range(n_s): 
        
        if not x_series[i] == 0:
            ncolors.append(colors[i])
            nseries_labels.append(series_labels[i])
            nx_series.append(x_series[i])
    val = 0        
    for i in range(len(nx_series)):    

            if nx_series[i]/sum(nx_series) < 0.1:
                val = val + 0.15     
            explode.append(val)
        
    
    #plt.pie(nx_series, explode = explode, labels = nseries_labels, colors = ncolors, autopct = '%1.1f%%', shadow = False, startangle = 140)
    plt.pie(nx_series, explode = explode, colors = ncolors, startangle = 140, radius = 0.6, autopct = '%1.2f%%', pctdistance = 1.2)
    plt.legend(labels = nseries_labels, loc=5, prop={'size':11}, frameon=False)
    if fname == '':
        plt.show()
    else:
        plt.savefig(fname)
        plt.close()
