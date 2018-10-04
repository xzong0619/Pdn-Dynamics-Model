# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 16:14:01 2018

@author: wangyf
"""

import random
import numpy as np 
import pickle
try:
    from mpi4py import MPI
except:
    pass

from datetime import datetime
from sklearn.metrics import mean_squared_error
import lattice_functions as lf
from numpy.linalg import norm
from test_connectivity_fitness import connect_score
#%%
def get_rank(COMM = None):
    if COMM is None:
        return 0
    else:
        return COMM.rank

def get_size(COMM = None):
    if COMM is None:
        return 1
    else:
        return COMM.size
    
def get_time():
    return datetime.now()   


def save_population(COMM = None, population = None, file_name = None):
	if COMM is None:
		rank = 0
	else:
		rank = COMM.rank

	if rank == 0: 
		with open(file_name, 'wb') as f_ptr:
			pickle.dump(population, f_ptr)
            
#%%
def occupancy():
    '''
    Creat occupancy for a node in the configuration, value between 0 or 1
    '''
    occ = random.randint(0, 1)
    return occ
    
    
def individual_config(individual, Clusters):
    
    ind_list = list(np.nonzero(individual)[0])
    Clusters.get_configs([ind_list])
    Gsv = Clusters.Gsv
    
    return Gsv

def evaluate_pi(individual, Clusters, Gcv):
    
    occ = Clusters.occupancy
    Gsv = individual_config(individual, Clusters)
    Cal = lf.calculations(occ)
    pi_pred =  Cal.get_pi_matrix(Gsv ,Gcv)
    
    return pi_pred
    
def evaluate(individual, Clusters, Gcv, pi_true, J, intercept):
    
    occ = Clusters.occupancy
    Gsv = individual_config(individual, Clusters)
    Cal = lf.calculations(occ)
    pi_pred =  Cal.get_pi_matrix(Gsv ,Gcv) 
    fitness1 = mean_squared_error(pi_pred, pi_true) 
    fitness2 = norm(pi_pred-pi_true, ord = np.inf)
    fitness3 = connect_score(individual)
    fitness4 = (np.dot(pi_pred, J) + intercept)[0]
    # possible to put lower energy clusters as fitness
    return (fitness1,fitness2,fitness3,fitness4)


def make_initial_population(COMM = None, toolbox = None, n = None):
    rank = get_rank(COMM)
    if rank == 0:
        print('\t{}  Core {}  Building initial population'.format(get_time(), rank))
        population = toolbox.population(n = n)
    else:
        population = None
    return population

def evaluate_population(COMM = None, toolbox = None, population = None):
    rank = get_rank(COMM)
    size = get_size(COMM)
    if rank == 0:
            jobs_split = np.array_split(range(len(population)), size)
            population_split = []
            for jobs in jobs_split:
                x = []
                for job in jobs:
                    x.append(population[job])
                population_split.append(x)
            print('\t{}  Core {}  Distributing individuals'.format(get_time(), rank))
    else:
        population_split = None
        jobs_split = None

    if COMM is None:
        population_mpi = population
        jobs_mpi = range(len(population))
    else:
        population_mpi = COMM.scatter(population_split, root = 0)
        jobs_mpi = COMM.scatter(jobs_split, root = 0)

    #Evaluate fitness
    fitnesses_mpi = {}
    for i, individual_mpi in zip(jobs_mpi, population_mpi):
        fitnesses_mpi[i] = toolbox.evaluate(individual_mpi)
    print('\t{}  Core {}  Finished evaluating individuals'.format(get_time(), rank))
    if COMM is None:
        fitnesses_list = [fitnesses_mpi]
    else:
        fitnesses_list = MPI.COMM_WORLD.gather(fitnesses_mpi, root = 0)
    if rank == 0:
        print('\t{}  Core {}  Assigning fitness to population.'.format(get_time(), rank))
        for fitnesses_dict in fitnesses_list:
            for i, fitness in fitnesses_dict.items():
                population[i].fitness.values = fitness
    else:
        population = None
    return population

def generate_offspring(COMM = None, toolbox = None, population = None, cxpb = None):
    rank = get_rank(COMM)
    if rank == 0:
        print( '\t{}  Core {}  Generating offspring'.format(get_time(), rank))
        print( '# Offspring before selection: {}'.format(len(population)))
        offspring = toolbox.select(population)
        print( '# Offspring after selection: {}'.format(len(offspring)))
        offspring = [toolbox.clone(individual) for individual in offspring]

        #Apply crossover and mutation on the offspring
        #print '\tMaking children'
        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            if random.random() < cxpb:
                toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values
    else:
        offspring = None
    return offspring

def mutate_offspring(COMM = None, toolbox = None, population = None):
    rank = get_rank(COMM)
    if rank == 0:
        print( '\t{}  Core {}  Mutating offspring'.format(get_time(), rank))
        #print '\tApplying mutations'
        for mutant in population:
            toolbox.mutate(mutant)
            del mutant.fitness.values
    else:
        population = None
    return population

def make_next_population(COMM = None, population = None, offspring = None):
    rank = get_rank(COMM)
    if rank == 0:
        print( '\t{}  Core {}  Generating new offspring'.format(get_time(), rank))
        population[:] = offspring
    return population


def calculate_statistics(COMM = None, generation = None, population = None):
    rank = get_rank(COMM)
    if rank == 0:
        print( '\t{}  Core {}  Calculating statistics'.format(get_time(), rank))
        fitnesses = get_fitnesses(population)
        avg = np.mean(fitnesses)
        sd = np.std(fitnesses)
        min_val = np.min(fitnesses)
        max_val = np.max(fitnesses)
        with open('stats.out', 'a') as f_ptr:
            f_ptr.write('{}\t{}\t{}\t{}\t{}\n'.format(generation, avg, sd, min_val, max_val))

def print_generation_number(COMM = None, generation = None):
    rank = get_rank(COMM)
    if rank == 0:
        print( '{}  Core {}  Generation {}'.format(get_time(), rank, generation))

def find_best_individual(COMM = None, population = None):
    rank = get_rank(COMM)
    if rank == 0:
        fitnesses = get_fitnesses(population)
        i = np.where(fitnesses == min(fitnesses))[0][0]
        print( '\tIndividual with best fitness:')
        print( '\tFitness = {} '.format(population[i].fitness.values))
        #print( '\tCV RMSE = {} eV'.format(np.sqrt(population[i].fitness.values[0])))

def get_fitnesses(population):
    '''
    Add the sum of the fitness values
    '''
    fitness_tuple = abs(np.array([individual.fitness.values for individual in population]))
    fitness_max = np.max(fitness_tuple, axis = 0)
    normalized_fit = fitness_tuple/fitness_max
    normalized_fit[:,-1] = 1-normalized_fit[:,-1]
    
    return normalized_fit.sum(axis = 0)


def final_best_individual(population, Clusters, Gcv):
    
    fitnesses = get_fitnesses(population)
    i = np.where(fitnesses == min(fitnesses))[0][0]
    
    best_ind = population[i]
    best_fitness = best_ind.fitness.values
    best_pi = evaluate_pi(best_ind, Clusters, Gcv)
    best_config = individual_config(best_ind, Clusters)
    
    return(best_ind, best_fitness, best_pi, best_config)
    
    