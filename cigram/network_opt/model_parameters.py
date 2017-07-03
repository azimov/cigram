'''
CiGRAM - Gaussian Random grAph Model

Copyright (C) 2014  James Gilbert

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
more details: http://www.gnu.org/licenses/gpl.html
'''
from __future__ import division
import networkx as nx
import json
import multiprocessing
import numpy as np
import random
import os
from cigram.network_opt.fitness_measures import summary_stat_fitness
import time
import logging
logger = logging.getLogger('inspyred.ec')

try:
        import inspyred
except AttributeError:
        # Python 2.6 backwards compaitability
        from ordereddict import OrderedDict
        import collections
        collections.OrderedDict = OrderedDict
        import inspyred



def pso(topology=inspyred.swarm.topologies.star_topology):
    '''
    create PSO optimiser instance
    '''
    ea = inspyred.swarm.PSO(random) 
    ea.terminator = inspyred.ec.terminators.evaluation_termination
    ea.topology = topology
    return ea

def dea():
    ea = inspyred.ec.DEA(random)
    ea.terminator = inspyred.ec.terminators.evaluation_termination
    return ea


def sim_anneal():
    ea = inspyred.ec.SA(random)
    ea.terminator = inspyred.ec.terminators.evaluation_termination
    return ea

model_store = {}

def graph_fit_evaluator(candidates, args):
    """
    Compare the distance between target and real network
    This is a generic function that runs the get_fittness method found in the network model class.
    """
    model = model_store[args["name"]]	
    fit = []
    for candidate in candidates:
            dist = []
            for _ in range(args["candidate_replicates"]):
                    distance = model.get_fitness(candidate)
                    dist.append(distance)
            fit.append(np.mean(dist))

    return fit

def graph_fit_evaluator_avg(candidates, args):
    """
    Compare the distance between target and real network with distance that uses averages over num replicates
    This is a generic function that runs the get_fittness method found in the network model class.
    """
    model = model_store[args["name"]]	
    fit = []
    for candidate in candidates:
            fit.append(model.get_fitness_avg(candidate, args["candidate_replicates"]))
    
    return fit


def time_logger_observer(population, num_generations, num_evaluations, args):
    """
    Save each generation of the optimi
    """
    e_time = time.time() - args["start_time"]
    save_path = args["observer_save_path"]
    model = model_store[args["name"]]
    if "best_fit_history" not in args:
            args["best_fit_history"] = []

    if "observer_times" not in args:
            args["observer_times"] = []

    if "eval_counts" not in args:
            args["eval_counts"] = []
    
    args["eval_counts"].append(num_evaluations)
    # get the population
    results = sorted([(model.get_params(c.candidate), c.fitness) for c in population], key=lambda x: x[1])
            
    best = results[0]
    args["best_fit_history"].append(best[1])
    args["observer_times"].append(e_time)
    observer_results = {
        "stats":{
            "num_generations": num_generations,
            "num_evaluations": num_evaluations,
            "observer_times": args["observer_times"],
            "best_fit": best[1],
            "best_fit_history": args["best_fit_history"],
            "eval_counts":args["eval_counts"]
        },
        "population":results
    }
    
    with open(save_path, "w+") as obs_file:
        json.dump(observer_results, obs_file, indent=4)

def default_observer(population, num_generations, num_evaluations, args):
    """
    Save each generation of the optimi
    """
    save_path = args["observer_save_path"]
    model = model_store[args["name"]]

    # get the population
    results = sorted([(model.get_params(c.candidate), c.fitness) for c in population], key=lambda x: x[1])
    pbest = results[0]
    
    model.set_gbest(pbest)

    observer_results = {
        "stats":{
            "num_generations": num_generations,
            "num_evaluations": num_evaluations,
            "best_fit": model.gbest[1],
            "best":model.gbest[0],
            
            "current_best_fit": pbest[1],
            "current_best":pbest[0],
        }
    }

    with open(save_path, "w+") as obs_file:
            json.dump(observer_results, obs_file, indent=4)

def optimise_model(model, num_cpus=multiprocessing.cpu_count(), pop_size=50, max_evals=10000, num_replicates=5, observer=default_observer, ea=pso(), seed=None, start_time=time.time(), evaluator=graph_fit_evaluator):
    '''
    Generic function to optimise a model to fit a target graph
    
    The objective is to minimise the distance between the target graph and some 
    Use specified model class which provides the interface for function parameters, and ranges used by the optimiser.
    
    Returns sorted list of model parameters with observed fitness. 
    
    Parameters:
            model gram_optimiser object that contains methods required for this process to work
            num_cpus: number of processes to be spawned - defaults to number of cpus
            
            pop_size: swarm size
            max_evals: number of max_evaluations
            num_replicates: number of candidate_replicates
            observer: optional bool to save results as the process is running
            seed: random seed
    '''
    np.random.seed(seed)
    random = np.random
    ea.observer = observer
    
    try:
        os.mkdir(".cigram_cache/")
    except OSError:
        # Dir already exists
        pass
    
    model_store[model.name] = model
    
    optimisation_parameters = dict(
        pop_size=pop_size,
        maximize=False,
        max_evaluations=max_evals,
        candidate_replicates=num_replicates,
        bounder=model.bounder(),
        generator=model.parameter_generator,
        name=model.name,
        observer_save_path=".cigram_cache/{0}.json".format(model.name),
        start_time = start_time
    )

    # This step is largely because the single process code is far easier to debug
    if num_cpus > 1:
        optimisation_parameters["evaluator"] = inspyred.ec.evaluators.parallel_evaluation_mp
        optimisation_parameters["mp_num_cpus"] = num_cpus
        optimisation_parameters["mp_evaluator"] = evaluator
    else:
        optimisation_parameters["evaluator"] = evaluator
            
    
    final_pop = ea.evolve(**optimisation_parameters)
    
    results =  sorted([model.gbest] + [(model.get_params(c.candidate), c.fitness) for c in final_pop], key=lambda x: x[1])
    
    return results


import cma
import multiprocessing as mp
        

def cma_es_optimise_model(model, num_cpus=mp.cpu_count(), pop_size=50, cma_sigma=0.5, max_evals=10000, num_replicates=5, seed=None, start_time=time.time(), evaluator=graph_fit_evaluator):
    
    '''
    Best fit optimisation with CMA-ES
    
     * Figure out how to set number of max_evaluations/time
     * Returning results of best fit networks
    '''

    # Start optimiser with random starting parameters
    es = cma.CMAEvolutionStrategy(model.parameter_generator(np.random, {}), cma_sigma, { "verb_disp":0, "verbose":0})
    
    def fitness_wrap(x, res, lock):
        fit = model.get_fitness(x)
        # Weird copying of variables is a necissary work around
        lock.acquire()
        result = res
        result.append((x, fit))
        res = result
        lock.release()

    
    manager = mp.Manager()
    lock = manager.Lock()
    
    while not es.stop():
        # First, get canidate solutions
        candidates = es.ask()
        # Multithreaded calulcation for each canidate solution
        
        f_values = manager.list()
        
        procs = []
        for x in candidates:

            proc = mp.Process(target=fitness_wrap, args=(x,f_values,lock))
            proc.start()
            procs.append(proc)
        
        for proc in procs:
            proc.join()
        
        # pass results to cma-es algorithm
        es.tell(*zip(*list(f_values)))
        es.disp()
        es.logger.add()
    
    return es