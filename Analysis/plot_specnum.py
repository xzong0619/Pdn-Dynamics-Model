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
    
    plt.xticks(size=20)
    if max(x_series) < 0.01:
        plt.ticklabel_format(style = 'sci', axis = 'x', scilimits= (0, 0))
        mat.rc('font', size = 20 )
    plt.yticks(size=20)
    plt.xlabel(xlab, size=24)
    plt.ylabel(ylab, size=24)
	
    if not xlimit == []:
        plt.xlim(xlimit)
	
    if not series_labels == []:
        plt.legend(series_labels, 
                   bbox_to_anchor = (1.02,1),loc= 'upper left',
                   prop={'size':20}, frameon=False)
    plt.tight_layout()
    
    
    if logscale:
        plt.yscale('log')
    
    if fname == '':
        plt.show()
    else:
        plt.savefig(fname, bbox_inches = "tight")
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
    plt.legend(labels = nseries_labels, 
               bbox_to_anchor = (0.8, 0.8),loc= 'upper left',
               prop={'size':11}, frameon=False)
    if fname == '':
        plt.show()
    else:
        plt.savefig(fname, bbox_inches = "tight")
        plt.close()


def PlotFreqs(event_freqs_vec, fname = ''):
    
    '''
    Plot a bar graph of elementary step frequencies versus time - output in elem_step_freqs.png in the directory with the Zacros run

    '''
    PlotOptions()
    plt.figure(figsize=(10,8))
    
    width = 0.2
    ind = 0
    yvals = []
    ylabels = []
    bar_vals = []
    
    store_ind = 0       # index of where the data is stored
    n_rxn =  int(len(event_freqs_vec)/2)
    
    for r in range(n_rxn):
   
        fwd_rate = event_freqs_vec[store_ind]
        store_ind += 1
        bwd_rate = event_freqs_vec[store_ind]   
        store_ind += 1
        
        if fwd_rate + bwd_rate > 0:
            
            net_freq = abs(fwd_rate - bwd_rate)
 
            if fwd_rate > 0:              
                plt.barh(ind-0.4, fwd_rate, width, color='r', log = True)
                bar_vals.append(fwd_rate)
            if bwd_rate > 0:
                plt.barh(ind-0.6, bwd_rate, width, color='b', log = True)
                bar_vals.append(bwd_rate)
            if net_freq > 0:
                plt.barh(ind-0.8, net_freq, width, color='g', log = True)
                bar_vals.append(net_freq)
            ylabels.append(str(r+1))
            yvals.append(ind-0.6)
            ind = ind - 1

    bar_vals = np.array(bar_vals)
    
    log_bar_vals = np.log10(bar_vals)
    xmin = 10**np.floor(np.min(log_bar_vals))
    xmax = 10**np.ceil(np.max(log_bar_vals))
            
    plt.xticks(size=20)
    plt.yticks(size=18)
    plt.xlabel('Frequency',size=20)
    plt.yticks(yvals, ylabels)
    r_patch = mat.patches.Patch(color = 'red', label = 'fwd')
    b_patch = mat.patches.Patch(color = 'blue', label = 'rev')
    g_patch = mat.patches.Patch(color = 'green', label = 'net')
    
    plt.legend(handles = [r_patch, b_patch, g_patch ],
               bbox_to_anchor = (1.05, 1),loc= 'upper left',
               prop={'size':18},frameon=False)
    
    plt.xlim([xmin, xmax])        
    plt.tight_layout()
    
    
    if fname == '':
        plt.show()
    else:
        plt.savefig(fname, bbox_inches = "tight")
        plt.close()
    