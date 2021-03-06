#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 14:48:22 2018

@author: wangyifan
"""

'''
Test_GA version two
Test on Pd20 
'''

import numpy as np
import pickle
from ase.visualize import view

try:
	from mpi4py import MPI
except:
	COMM = None
	rank = 0
else:
	COMM = MPI.COMM_WORLD
	rank = COMM.rank
    

from deap import base
from deap import creator
from deap import tools

from GA_functions import get_time
[Gcv_nonzero, J_nonzero, 
 intercept, MSE_test, MSE_train,
 pi_nonzero,y] =  pickle.load(open("lasso.p", "rb"))


import GA_functions as GA
import lattice_functions as lf



#%%
#Genetic Hyperparameters
ngoal = 21 #5, 10, 15, 20
nodes = 36 #lattice node size
n = 100 #Size of population
ngen = 2 #Number of generations
cxpb = 0.8 #The probability of mating two individuals
mutpb = 0.1 #The probability of mutating an individual
k = n
tournsize = 10

score_weights = (-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0) #tuple for min-1.0, max+1.0

print('{}  Core {}  Reading files'.format(get_time(), rank))

#Create the fitness object
creator.create("FitnessMin", base.Fitness, weights = score_weights)
#Create an individual
creator.create("Individual", list, fitness = creator.FitnessMin)

#Create the toolbox
toolbox = base.Toolbox()

toolbox.register("attr_binary", GA.occupancy)
toolbox.register("individual", tools.initRepeat, creator.Individual,
                 toolbox.attr_binary, n=nodes)

toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb = mutpb)
toolbox.register("select", tools.selTournament, k = k, tournsize = tournsize)
toolbox.register("evaluate", GA.evaluate, ngoal= ngoal)

population = GA.make_initial_population(COMM, toolbox, n)
population = GA.evaluate_population(COMM, toolbox, population)
history = []


#%%
ind1 = toolbox.individual()
fit1 = toolbox.evaluate(ind1)
pi1 = GA.evaluate_pi(ind1) 

ind2 = toolbox.individual()
fit2 = toolbox.evaluate(ind2)
pi2 = GA.evaluate_pi(ind2) 
#%%
sim_mat = np.zeros(shape = (n, n))
for generation in range(ngen):
    GA.print_generation_number(COMM, generation)
    offspring = GA.generate_offspring(COMM, toolbox, population, cxpb)
    offspring = GA.mutate_offspring(COMM, toolbox, offspring)
    offspring = GA.evaluate_population(COMM, toolbox, offspring)
    population = GA.make_next_population(COMM, population, offspring)
    GA.calculate_statistics(COMM, generation, population)
    GA.find_best_individual(COMM, population)
    history = GA.write_history(population, history)
    
    
#%%
n_hof  = 20
ihof,E_hof  = GA.hall_of_fame(COMM, history, 20)
(best_ind, best_fitness, best_pi, best_config) = GA.winner_details(COMM, ihof)
lf.drawing(best_config[0])
best_G = GA.individual_config(best_ind)
GA.ase_object(best_ind)
ind_list = list(np.nonzero(best_ind)[0])
pickle.dump([ihof, E_hof], open('pd_'+str(ngoal) + '.p','wb'))


