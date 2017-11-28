"""
Module contains classes for finding best fit GRAM parameters in the parameter selection tool.
"""
from __future__ import division
import networkx as nx
from cigram import cigram_graph
import numpy as np

import inspyred
from collections import OrderedDict

import logging
logger = logging.getLogger('inspyred.ec')


class CigramOptimiser(object):
    
    def __init__(self, g, fit_evaluator,
                 a_bounds=(-5.0, 5.0),
                 sf_bounds=(0.6, 2.0),
                 sr_bounds=(0.6, 2.0),
                 csf_bounds=(0.8, 2.0),
                 csr_bounds=(0.8, 2.0),
                 ek_per_bounds=(0.01, 0.4),
                 po_bounds=(0.01, 0.15),
                 max_n=None,
                 k_bounds=None,
                 random=np.random,
                 min_degree=1,
                 starting_population=None,
                 graph_props=None,
                 seed=None,
                 **kwargs):
        """
        Optimisation class to fit a given type of graph model.
        This code of this class could be generalised to another type of model to adapt the fitting procedure


        :param g: networkx.Graph object to fit
        :param fit_evaluator:
        :param a_bounds:
        :param sf_bounds:
        :param sr_bounds:
        :param csf_bounds:
        :param csr_bounds:
        :param ek_per_bounds:
        :param po_bounds:
        :param max_n:
        :param k_bounds:
        :param random:
        :param min_degree:
        :param starting_population:
        :param graph_props:
        :param seed:
        :param kwargs:
        """

        self.graph = g
        self.fixed_model_params = {
                "n": g.number_of_nodes(),
                "density": nx.density(g),
                "connected": False,
                "min_degree": min_degree
        }
        
        self.gbest = None

        if starting_population is None:
            starting_population = []

        self.starting_population = starting_population
        
        if max_n is not None and max_n < g.number_of_nodes():
                scaling_factor = g.number_of_nodes()/max_n
                self.fixed_model_params["n"] = max_n
                self.fixed_model_params["density"] *= scaling_factor

        if k_bounds is None:
                k_bounds = (1, int(self.fixed_model_params["N"]/20))
        
        # Fix other parameters if we're not generating community strucuture
        if k_bounds == 1:
                csf_bounds = 0
                csr_bounds = 0
                ek_per_bounds = 0
                po_bounds = 0

        # Set the boundaries of the optimisation that map to the parameters.
        self.param_boundaries = OrderedDict()
        self.param_boundaries["sigma_nodes"] = sf_bounds
        self.param_boundaries["sigma_edges"] = sr_bounds
        self.param_boundaries["a"] = a_bounds
        self.param_boundaries["k"] = k_bounds
        self.param_boundaries["community_sigma_f"] = csf_bounds
        self.param_boundaries["community_sigma_r"] = csr_bounds
        self.param_boundaries["ek_per"] = ek_per_bounds
        self.param_boundaries["p_o"] = po_bounds
        
        self.random_param_func = dict([(p, random.uniform) for p in self.param_boundaries.keys()])
        self.random_param_func["k"] = random.randint
        
        for p, bounds in self.param_boundaries.items():
                if type(bounds) is not tuple:
                        self.fixed_model_params[p] = bounds
                        del self.param_boundaries[p]
                elif bounds[0] == bounds[1]:
                        self.fixed_model_params[p] = bounds[0]
                        del self.param_boundaries[p]
        
        if "name" in kwargs:
                self.name = kwargs["name"]
        else:
                self.name = "{0}".format(g.name)
        
        self.fit_evaluator = fit_evaluator
        # Store the target degree distribution/assortativity/clustering
        # can be already specified if user wants to fit e.g. a degree distribution not obtained from a graph
        if graph_props is not None:
                self.graph_props_a = graph_props
        else: 
                self.graph_props_a = self.fit_evaluator.graph_properties(g)
                
        self.seed = seed
    
    def set_gbest(self, pbest):
        """
        Updates the best solution
        :param pbest: solution
        :return: None
        """
        if self.gbest is None or self.gbest[1] > pbest[1]:
                self.gbest = pbest
    
    def set_name(self, name):
        """
        The name to give this optimiser object.
        Useful in the case of multiple optimisations of different models that can be resumed part way through
        """
        self.name = name

    def generate(self, candidate):
        """
        Generate the model given a candidates parameter set.
        Must return a graph object
        """
        params = dict(self.fixed_model_params.items() + self.convert_params(candidate).items())
        
        if self.seed is not None:
                params["seed"] = self.seed
                self.seed += 1

        logger.debug("Generating {0}".format(params))
        g, _, _ = cigram_graph(**params)
        logger.debug("Generated {0}".format(params))
        return g

    def parameter_generator(self):
        """
        Generates a random set of parameters for use by the optimiser
        if the starting population kwarg is used at init, the population is not randomised but starts with specified
        values
        """
        
        if len(self.starting_population):
            cand = self.starting_population.pop()
            pset = [cand[p] for p in self.param_boundaries.keys()]
        else:
            pset = [self.random_param_func[p](*bounds) for p, bounds in self.param_boundaries.items()]
        return pset

    def convert_params(self, candidate):
        """
        A candidate set will be in a given order (specified by the parameter_generator method)
        This function just converts the params to those used by the model generator function.
        """
        keys = self.param_boundaries.keys()
        return dict((keys[i], candidate[i]) for i in range(len(candidate)))

    def get_candidate(self, parameters):
        """
        Given params, convert to a given candidate format
        """
        return [parameters[p] for p in self.param_boundaries]

    def get_params(self, candidate):
        
        return dict(self.fixed_model_params.items() + self.convert_params(candidate).items())

    def get_fitness(self, candidate):
        """
        given a set of parameters, evaluate the fitness
        This function tests the assortativity and distance of the degree distribution (using the KS and log difference
        of max degree)
        Create a subclass of this class to change this function.
        """
        m = self.generate(candidate)
        graph_props_b = self.fit_evaluator.graph_properties(m)
        del m
        return self.fit_evaluator.graph_fitness(self.graph_props_a, graph_props_b)
    
    def get_fitness_avg(self, candidate, num_reps=1):
        """
        given a set of parameters, evaluate the fitness over average range of samples
        This function tests the assortativity and distance of the degree distribution (using the KS and log difference
        of max degree)
        Create a subclass of this class to change this function.
        """
        
        if num_reps < 1:
                num_reps = 1
        graphs = [self.generate(candidate) for _ in range(num_reps)]
        graph_props_b = self.fit_evaluator.graph_properties_avg(graphs)
        return self.fit_evaluator.graph_fitness(self.graph_props_a, graph_props_b)

    def bounder(self):
        """
        Create the bounder object for the model - required by inspyred optimiser
        """
        lower_constraints = tuple([v[0] for v in self.param_boundaries.values()])
        upper_constraints = tuple([v[1] for v in self.param_boundaries.values()])
        
        return inspyred.ec.Bounder(lower_constraints, upper_constraints)
