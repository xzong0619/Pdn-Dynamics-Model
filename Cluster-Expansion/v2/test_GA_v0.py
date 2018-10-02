# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 11:51:22 2018

@author: wangyf
"""

# test the genetic algorithm 
import numpy as np
import random
from sklearn.metrics import mean_squared_error
from deap import base
from deap import creator
from deap import tools


import lattice_functions as lf
from structure_constants import mother, dz
from test_lasso import Gcv_nonzero
from test_lasso_krige import newpoints

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

def occupancy():
    '''
    Creat occupancy for a node in the configuration, value between 0 or 1
    '''
    occ = random.randint(0, 1)
    return occ
    
    
def individual_config(individual):
    
    ind_list = list(np.nonzero(individual)[0])
    Clusters.get_configs([ind_list])
    Gsv = Clusters.Gsv
    
    return Gsv

def evaluate(individual, Gcv, pi_true):
      
    Gsv = individual_config(individual)
    Cal = lf.calculations(occ)
    pi_pred =  Cal.get_pi_matrix(Gsv ,Gcv) 

    fitness = mean_squared_error(pi_pred, pi_true) 
    # possible to put lower energy clusters as fitness
    return (fitness,)

#Genetic Hyperparameters
n = 100 #Size of population
ngen = 100 #Number of generations
cxpb = 1.0 #The probability of mating two individuals
mutpb = 0.001 #The probability of mutating an individual
k = n
tournsize = 5

IND_SIZE = 36

#Create the fitness object
creator.create("FitnessMin", base.Fitness, weights = (-1.0,))
#Create an individual
creator.create("Individual", list, fitness = creator.FitnessMin)

toolbox = base.Toolbox()
toolbox.register("attr_binary", occupancy)
toolbox.register("individual", tools.initRepeat, creator.Individual,
                 toolbox.attr_binary, n=IND_SIZE)


toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb = mutpb)
toolbox.register("select", tools.selTournament, k = k, tournsize = tournsize)
toolbox.register("evaluate", evaluate, Gcv =  Gcv_nonzero, pi_true = newpoints)

pop = toolbox.population(n = n)

ind1 = toolbox.individual()
print(ind1)               # [0.86..., 0.27..., 0.70..., 0.03..., 0.87...]



ind_fit = evaluate(ind1, Gcv_nonzero, newpoints)


#%%
for g in range(ngen):
    # Select the next generation individuals
    offspring = toolbox.select(pop)
    # Clone the selected individuals
    offspring = [toolbox.clone(individual) for individual in offspring]

    # Apply crossover on the offspring
    for child1, child2 in zip(offspring[::2], offspring[1::2]):
        if random.random() < cxpb:
            toolbox.mate(child1, child2)
            del child1.fitness.values
            del child2.fitness.values

    # Apply mutation on the offspring
    for mutant in offspring:
        if random.random() < mutpb:
            toolbox.mutate(mutant)
            del mutant.fitness.values

    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
    fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit

    # The population is entirely replaced by the offspring
    pop[:] = offspring
    

