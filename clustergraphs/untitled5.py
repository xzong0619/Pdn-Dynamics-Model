# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 13:37:09 2018

@author: wangyf
"""

def bubbleSort(arr):
    n = len(arr)
 
    # Traverse through all array elements
    for i in range(n):
 
        # Last i elements are already in place
        for j in range(1, n-i):
 
            # traverse the array from 0 to n-i-1
            # Swap if the element found is greater
            # than the next element
            if arr[j] > arr[j-1] :
                arr[j-1], arr[j] = arr[j], arr[j-1]
                
                
# Driver code to test above
arr = [3333, 1111111, 0.1, 2, 3]
 
bubbleSort(arr)
 
print ("Sorted array is:")
for i in range(len(arr)):
    print ("%d" %arr[i])