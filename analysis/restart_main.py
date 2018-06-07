# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 13:49:38 2018

@author: wangyf
"""
import sys
sys.path.append(r"C:\Users\wangyf\Documents\GitHub\Pdn-Dynamics-Model\Generator")
sys.path.append(r"C:\Users\wangyf\Documents\GitHub\Pdn-Dynamics-Model\Analysis")
from IO_Pdn import *
from process_specnum import *

'''
User Input section (reaction conditions):
--------------------------------------------------------------
'''
nc = [None,  # neighboring structure
      '1-2',
      '1-2 1-3 2-3',
      '1-2 1-3 2-3 4-2 4-3 ',
      '1-2 1-3 2-3',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 ',
      '1-2 1-3 2-3 4-2 4-3 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-3 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-4 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-4 6-2 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 ',
      '1-2 1-3 2-3 4-2 4-3 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-3 7-6 7-3 1-7 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-4 7-6 7-5 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-3 7-6 7-5 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-4 6-2 7-6 7-2 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-3 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-5 6-4 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 6-4 6-2 ',
      '1-2 1-3 2-3 4-2 4-3 5-4 5-3 ']

Base_path = os.getcwd()
input_dir = os.path.join(Base_path, 'zacros_inputs')
output_dir = os.path.join(Base_path, 'outputs')

State =  StateOut(nc)
end_state = restart(input_dir)

s_name = end_state.s_name
s_n = end_state.s_n
s_count = end_state.s_count


State.WriteIn(output_dir, s_name, s_n, s_count)