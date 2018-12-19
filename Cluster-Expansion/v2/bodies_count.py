# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 12:00:10 2018

@author: wangyf
"""
bodies_count = np.array([1, 3, 32])
bodies = ['Points','Pairs', 'Triplets']
base_line = 0
x_pos = np.arange(len(bodies_count))
opacity = 1

plt.figure(figsize=(4,3))
rects1 = plt.bar(x_pos, bodies_count - base_line, #yerr=std_train,
                alpha=opacity, color='lightblue',
                label='Train')
#plt.ylim([-1,18])
plt.xticks(x_pos, bodies)
plt.xlabel('Clusters')
plt.ylabel('Count')


