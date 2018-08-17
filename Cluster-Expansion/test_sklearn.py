# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 13:10:08 2018

@author: wangyf
"""

from sklearn.model_selection import LeaveOneOut

X = [1,2,3,4]
loo = LeaveOneOut()

for train,test in loo.split(X):
    print("%s %s" % (train, test))