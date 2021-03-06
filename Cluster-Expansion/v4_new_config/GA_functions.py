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
from test_connectivity_fitness import connect_score_2
from ase.io import read, write
from ase.visualize import view
from ase.data import covalent_radii
from ase import Atoms, Atom
from ase.build import surface
from structure_constants import mother, dz

[Gcv_nonzero, J_nonzero, 
 intercept, MSE_test, MSE_train,
 pi_nonzero,y] =  pickle.load(open("lasso.p", "rb"))

#%%
def initialize_Clusters_object():
    
    empty = 'grey'
    filled = 'r'
    occ = [empty, filled]
    
    '''
    only draw 1st nearest neighbors?
    '''
    NN1 = 0
    '''
    Draw mother/conifgurations/clusters?
    '''
    draw = [0, 0, 0]
    
    
    Clusters = lf.clusters(occ, NN1, draw)
    Clusters.get_mother(mother, dz)
    
    return Clusters

Clusters = initialize_Clusters_object()

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

def one_hot_to_index(individual):
    '''
    Convert an individual from one hot encoding to a list index
    '''
    ind_index = list(np.nonzero(individual)[0])
    return ind_index

def index_to_one_hot(ind_index, n_nodes = 36):
    '''
    Convert an individual from a list index to one hot encoding  
    '''
    individual = np.zeros(n_nodes, dtype = int)
    individual[np.array(ind_index)] = 1
    individual = list(individual)
    
    return individual 

    
def ase_object(individual, view_flag = False):
    '''
    take in the 1/0 vector for index
    '''
    ind_index = one_hot_to_index(individual)
    Pdr = covalent_radii[46]
    Or = covalent_radii[8]
    PdPd = Pdr*2
    PdO = Pdr + Or
    Pdnodes = mother[ind_index]
    OO = 3.882
    
    def ceria():
        #Lattice constant
        a = 5.49
        CeO2 = Atoms('Ce4O8', scaled_positions =[ (0., 0., 0.),
                      (0., 0.5, 0.5),
                      (0.5, 0., 0.5),
                      (0.5, 0.5, 0.),
                      (0.75, 0.25, 0.25),
                      (0.25, 0.75, 0.75),
                      (0.75, 0.75, 0.75),
                      (0.25, 0.25, 0.25),
                      (0.25, 0.25, 0.75),
                      (0.75, 0.75, 0.25),
                      (0.25, 0.75, 0.25),
                      (0.75, 0.25, 0.75)],
                      cell = [a,a,a],
        			  pbc = True	)
        #Scales the atomic positions with the unit cell
        #CeO2.set_cell(cell, scale_atoms=True)
        #(1,1,1) is the slab type. There are 2 unit cells along 
        #z direction
        slab = surface(CeO2, (1, 1, 1), 2)
        
        #Repeating the slab 2 unit cells in x and 1 unit cell 
        #in y directions
        slab = slab.repeat((3,3,1))
        slab.center(vacuum=10, axis=2)
        del slab[[atom.index for atom in slab if atom.z>15]]
        
        return slab
    
    support = ceria()
    #view(support)

    origin = support[89].position
    origin[2] = origin[2] +  PdO
    Pdm = Pdnodes*PdPd + origin
    #PdNP = Atoms('Pd36', positions = Pdm)
    nPd = len(Pdnodes)
    for i in range(nPd):
        support.append(Atom('Pd', position = Pdm[i]))
    
    if view_flag: view(support) 
    
    return support
                
#%%
def occupancy():
    '''
    Creat occupancy for a node in the configuration, value between 0 or 1
    '''
    occ = random.randint(0, 1)
    return occ
    
    
def individual_config(individual, Clusters = Clusters ):
    
    ind_index = one_hot_to_index(individual)
    Clusters.get_configs([ind_index])
    Gsv = Clusters.Gsv
    
    return Gsv

def evaluate_pi(individual, Clusters = Clusters, Gcv =  Gcv_nonzero):
    
    occ = Clusters.occupancy
    Gsv = individual_config(individual, Clusters)
    Cal = lf.calculations(occ)
    pi_pred =  Cal.get_pi_matrix(Gsv ,Gcv)
    
    return pi_pred

def predict_E(ind_index, Clusters = Clusters, Gcv =  Gcv_nonzero, J = J_nonzero, intercept = intercept):
    '''
    Predict Energy of the cluster only, take in ind index
    '''
    individual = index_to_one_hot(ind_index)
    occ = Clusters.occupancy
    Gsv = individual_config(individual, Clusters)
    Cal = lf.calculations(occ)
    pi_pred =  Cal.get_pi_matrix(Gsv ,Gcv)
    E_pred = float(np.dot(pi_pred, J) + intercept)
    
    return E_pred
    
    
def evaluate(individual, ngoal,  Clusters = Clusters, Gcv =  Gcv_nonzero, J = J_nonzero, intercept = intercept):
    
    occ = Clusters.occupancy
    Gsv = individual_config(individual, Clusters)
    occ_nodes = one_hot_to_index(individual)
    Cal = lf.calculations(occ)
    pi_pred =  Cal.get_pi_matrix(Gsv ,Gcv)
    n_node = sum(individual)
    
    #fitness1 = mean_squared_error(pi_pred, pi_true)/4
    #fitness2 = norm(pi_pred-pi_true, ord = np.inf)
    #fitness3 = connect_score(individual)
    #fitness4 = (np.dot(pi_pred, J) + intercept)[0]
    # possible to put lower energy clusters as fitness
    #return (fitness1,fitness2,fitness3,fitness4)
    fitness1 = abs(n_node-ngoal)
    fitness2, fitness3, fitness4, fitness5, fitness6 = connect_score_2(occ_nodes)
    fitness7 = float(np.dot(pi_pred, J) + intercept)
    
    if ngoal <= 10:
        return fitness1,  fitness2, fitness3, fitness6, fitness7, fitness4, fitness5
    if ngoal > 10 and ngoal <= 15:
        return fitness1,  fitness2, fitness3, fitness4, fitness7, fitness5, fitness6
    if ngoal > 15:
        return fitness1,  fitness7, fitness2, fitness3, fitness4, fitness5, fitness6


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
        print('\t{} Core {} Before gather.'.format(get_time(), rank))
        fitnesses_list = MPI.COMM_WORLD.gather(fitnesses_mpi, root = 0)
        print('\t{} Core {} After gather.'.format(get_time(), rank))
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

def find_best_individual(COMM = None, population = None, nbest = 1):
    rank = get_rank(COMM)
    if rank == 0:
        fitnesses = get_fitnesses(population)
        fiv = []
        for fi in range(fitnesses.shape[1]):
            fiv.append(fitnesses[:,fi])
           
        ind = np.lexsort((fiv[-1], fiv[-2], fiv[-3], fiv[-4], fiv[-5], fiv[-6], fiv[-7]))
        best_i = ind[0]
        
        print( '\tIndividual with best fitness:')
        print( '\tFitness = {} '.format(population[best_i].fitness.values))
        #print( '\tCV RMSE = {} eV'.format(np.sqrt(population[i].fitness.values[0])))
    return best_i


def find_best_individuals(COMM = None, population = None, nbest = 1):
    
    rank = get_rank(COMM)
    if rank == 0:
        fitnesses = get_fitnesses(population)
        fiv = []
        for fi in range(fitnesses.shape[1]):
            fiv.append(fitnesses[:,fi])
           
        ind = np.lexsort((fiv[-1], fiv[-2], fiv[-3], fiv[-4], fiv[-5], fiv[-6] , fiv[-7]))
        best_index = np.array([ind[0]])
        counter = 0 #Index of indiviuals
        count = 1 #count of qualifing individuals
    
        while count < nbest:
            best_pre = population[best_index[-1]]
            best_next = population[ind[counter]]
            if not best_pre.fitness.values == best_next.fitness.values:
                best_index = np.append(best_index, ind[counter])
                count = count + 1
                
            counter = counter + 1
       
    return best_index
            
    
def get_fitnesses(population = None):
    '''
    normalize the fitness values
    '''
    fitness_tuple = np.array([individual.fitness.values for individual in population])
    
    #fitness_max = np.max(fitness_tuple, axis = 0)
    #normalized_fit = fitness_tuple/fitness_max
    #normalized_fit[:,-1] = 1-normalized_fit[:,-1]
    
    return fitness_tuple

def write_history(population, history):
    for pi in population:
        history.append(pi)
    
    return history
    
def hall_of_fame(COMM = None, history = None, nbest = None):
    
    print( '\nHall of Fame top {}:'.format(nbest))
    ihof = find_best_individuals(COMM, history, nbest)
    hof = []
    E_hof = []
    for ih in ihof:
        hof.append(history[ih])
        E_hof.append(min(history[ih].fitness.values))
    for i in range(nbest):
        print( '\tIndividual {} with best fitness:  Fitness = {}'.format(i, hof[i].fitness.values))
    
    return hof, E_hof 



def winner_details(COMM = None, population = None, Clusters = Clusters, Gcv =  Gcv_nonzero):
    
    print( '\nWinner Details:')
    i = find_best_individual(COMM, population)
    
    best_ind = population[i]
    best_fitness = best_ind.fitness.values
    best_pi = evaluate_pi(best_ind, Clusters, Gcv)
    best_config = individual_config(best_ind, Clusters)
    
    return(best_ind, best_fitness, best_pi, best_config)
    
    