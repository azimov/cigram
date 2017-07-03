from __future__ import division
import inspyred
import numpy as np
from itertools import product

class arpso(inspyred.swarm.PSO):
    """
    Particle swarm optimisation with attractive repulsive behaviour
    
    Works as standard PSO, however, when diversity of the swarm is below a given threshold the algorithm
    repels particles until the diversity is high enough.
    
    This causes the particles to explore more of the search space than they would under the standard implementation
    
    """
    
    
    def __init__(self, random):
        inspyred.ec.EvolutionaryComputation.__init__(self, random)
        self.topology = inspyred.swarm.topologies.star_topology
        self._previous_population = []
        self.selector = self._swarm_selector
        self.replacer = self._swarm_replacer
        self.variator = self._swarm_variator
        self.archiver = self._swarm_archiver
        self.direction = 1
        self.length = None
        
    def _population_diversity(self, population, bounder):
        diversity = 0.0
        #Np = float(len(population))
        
        ### Compute the average candidate of the swarm
        #avg_candidate = np.mean([x.candidate for x in population], axis=0)
        
        
        
        ## Compute the diversity of the swarm
        #for x in population:
            #p_diversity = 0.0
            #for xi, avg_i, lowerb, upperb  in zip(x.candidate, avg_candidate, bounder.lower_bound, bounder.upper_bound):
                #p_diversity += ((xi - lowerb)/(upperb - lowerb) - (avg_i - lowerb)/(upperb - lowerb))**2
            
            #diversity += np.sqrt(p_diversity) 
        
        
        ### Compute the length of the longest search space diagonal
        ##if self.length is None:
            ##self.length = 0.0
            ##varset = range(len(bounder.lower_bound))
            ##for a, b in product(varset, varset):
                ##A = bounder.upper_bound[a] - bounder.lower_bound[a]
                ##B = bounder.upper_bound[b] - bounder.lower_bound[b]
                ##D = np.sqrt(A**2 + B**2)
                ##if D > self.length:
                    ##self.length = D
        
        ##length = self.length
        
        cart_prod = product(population, population)
        distance = []
        for (p, q) in cart_prod:
            d = 0
            for x, y in zip(p.candidate, q.candidate):
                d += (x - y)**2
            
            #print d
            distance += [np.sqrt(d)]
            
        #diversity *= 1/Np
        return distance
    
    def _swarm_variator(self, random, candidates, args):
        inertia = args.setdefault('inertia', 0.5)
        cognitive_rate = args.setdefault('cognitive_rate', 2.1)
        social_rate = args.setdefault('social_rate', 2.1)
        
        diversity_lower_threshold = args.setdefault('diversity_lower_threshold', 1.0)
        diversity_upper_threshold = args.setdefault('diversity_upper_threshold', 2.0)
        
        if len(self.archive) == 0:
            self.archive = self.population[:]
        if len(self._previous_population) == 0:
            self._previous_population = self.population[:]
    
            
        diversity = self._population_diversity(self.population, self.bounder)
        
        if self.direction == 1 and max(diversity) < diversity_lower_threshold:
            self.direction = -1
        elif self.direction == -1 and max(diversity) >= diversity_upper_threshold:
            self.direction = 1
        
        print max(diversity), self.direction
        
        neighbors = self.topology(self._random, self.archive, args)
        offspring = []
        for x, xprev, pbest, hood in zip(self.population, 
                                            self._previous_population, 
                                            self.archive, 
                                            neighbors):
            nbest = max(hood)
            particle = []
            for xi, xpi, pbi, nbi in zip(x.candidate, xprev.candidate, 
                                            pbest.candidate, nbest.candidate):
                value = xi + inertia * (xi - xpi) +  (self.direction * (cognitive_rate * random.random() * (pbi - xi))) +  (self.direction * (social_rate * random.random() * (nbi - xi)))

                particle.append(value)
            particle = self.bounder(particle, args)
            offspring.append(particle)
        return offspring