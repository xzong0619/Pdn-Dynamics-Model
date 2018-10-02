# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 10:17:00 2018

@author: wangyf
"""

# test the genetic algorithm 
import numpy as np


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
from test_lasso import Gcv_nonzero
from test_lasso_krige import newpoints
import GA_functions as GA
import lattice_functions as lf
from structure_constants import mother, dz


#%%
'''
Convert individual occupancy vector to configuration in graph
'''
empty = 'grey'
filled = 'r'
occ = [empty, filled]

'''
only draw 1st nearest neighbors?
'''
NN1 = 1
'''
Draw mother/conifgurations/clusters?
'''
draw = [0, 0, 0]


Clusters = lf.clusters(occ, NN1, draw)
Clusters.get_mother(mother, dz)
    
    

#%%
#Genetic Hyperparameters
n = 100 #Size of population
ngen = 100 #Number of generations
cxpb = 1.0 #The probability of mating two individuals
mutpb = 0.001 #The probability of mutating an individual
k = n
tournsize = 5


print('{}  Core {}  Reading files'.format(get_time(), rank))

#Create the fitness object
creator.create("FitnessMin", base.Fitness, weights = (-1.0,))
#Create an individual
creator.create("Individual", list, fitness = creator.FitnessMin)

#Create the toolbox
toolbox = base.Toolbox()

toolbox.register("attr_binary", GA.occupancy)
toolbox.register("individual", tools.initRepeat, creator.Individual,
                 toolbox.attr_binary, n=n)

toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb = mutpb)
toolbox.register("select", tools.selTournament, k = k, tournsize = tournsize)
toolbox.register("evaluate", GA.evaluate, Clusters = Clusters, occ = occ, Gcv =  Gcv_nonzero, pi_true = newpoints)

population = GA.make_initial_population(COMM, toolbox, n)
population = GA.evaluate_population(COMM, toolbox, population)

ind1 = toolbox.individual()
fit1 = toolbox.evaluate(ind1)
 
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
    GA.save_population(COMM, population, 'population{}.txt'.format(generation))
    