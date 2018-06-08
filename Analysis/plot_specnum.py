# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 17:59:55 2018

@author: wangyf
"""
import matplotlib as mat
import matplotlib.pyplot as plt


        
        

def PlotOptions():

    '''
    Set matplotlib options
    '''

    mat.rcParams['mathtext.default'] = 'regular'
    mat.rcParams['text.latex.unicode'] = 'False'
    mat.rcParams['legend.numpoints'] = 1
    mat.rcParams['lines.linewidth'] = 2
    mat.rcParams['lines.markersize'] = 12
            

def PlotTimeSeries(x_series, y_series, xlab = 'Time (s)', ylab = '', series_labels = [], fname = '', logscale = False):
    
    '''
    Plot multiple series against time
    '''
    
    PlotOptions()
    plt.figure(figsize=(10,8))
    
    for i in range (len(y_series)):
        plt.plot(x_series, y_series[i])
    
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.xlabel(xlab, size=24)
    plt.ylabel(ylab, size=24)
    
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
        
        
        
