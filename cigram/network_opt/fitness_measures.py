from __future__ import division
import networkx as nx
import numpy as np
from networkx.algorithms.approximation.clustering_coefficient import average_clustering
from cigram.utils.analysis import ks_distance

class summary_stat_fitness(object):
    '''
    Uses distance between different summary statistics to compute graph similiarity
    '''
    
    def __init__(self,
                max_degree_fit_weight=0.5,
                degree_fit_weight=1.0,
                clustering_fit_weight=0.5,
                assort_fit_weight=0.5,
                **kwargs):
        
        self.mdfw = max_degree_fit_weight
        self.dfw = degree_fit_weight
        self.cfw = clustering_fit_weight
        self.afw = assort_fit_weight
            
    
    def graph_properties(self, G):
        '''
        compute the properties of the graph in question
        
        max degree, ks_dist, degree_assortativity, clustering coefficient
        '''
        degree = G.degree().values()
        props =  dict(
                max_degree=max(degree),
                degree=degree,
                clustering  = average_clustering(G, trials=1000),
                assort = nx.degree_assortativity_coefficient(G)
        )

        return props
    
    def graph_fitness(self, prop_a, prop_b):
        '''
        Summary statistics dissimmilarity
        '''
        max_deg_dist = np.abs(np.log10(prop_a["max_degree"]) - np.log10(prop_b["max_degree"]))
        ks_dist = ks_distance(prop_a["degree"], prop_b["degree"], normed=False, d_min=10)
        assort_dist = np.abs(prop_a["assort"] - prop_b["assort"]) 
        clust_dist = np.abs(prop_a["clustering"] - prop_b["clustering"])

        return (ks_dist + 0.5 * (max_deg_dist + assort_dist + clust_dist))


class no_clustering_summary(object):
    '''
    Graph summary stat similarity without measuring the expensive to compute clustering coefficient
    '''
    def __init__(self,
                    max_degree_fit_weight=0.5,
                    degree_fit_weight=1.0,
                    assort_fit_weight=0.5,
                    **kwargs):
            
        self.mdfw = max_degree_fit_weight
        self.dfw = degree_fit_weight
        self.afw = assort_fit_weight
    
    
    def graph_properties(self, G):
        '''
        compute the properties of the graph in question (nx graph object)
        
        max degree, ks_dist, degree_assortativity,
        '''
        degree = G.degree().values()
        props =  dict(
                degree=degree,
                max_degree = max(degree),
                assort = nx.degree_assortativity_coefficient(G)
        )

        return props
    
    def graph_fitness(self, prop_a, prop_b):
        '''
        Summary statistics dissimmilarity
        '''
        max_deg_dist = np.abs(np.log10(prop_a["max_degree"]) - np.log10(prop_b["max_degree"]))
        ks_dist = ks_distance(prop_a["degree"], prop_b["degree"], normed=False, d_min=10)
        assort_dist = np.abs(prop_a["assort"] - prop_b["assort"]) 
        
        return (ks_dist + 0.5 * (max_deg_dist + assort_dist))


class degree_cdf_l2_norm(object):
    
    def __init__(self, **kwargs):
        pass
    
    
    def graph_properties(self, G):
        '''
        '''
        degree = G.degree()
        
    
        c_degree_diff = [ len([k for k in degree.values() if k <= x] )/float(len(degree)) for x in range(max(degree.values()) + 1)]
        
        return c_degree_diff
    
    
    def graph_fitness(self, degree_diff_a, degree_diff_b):
        '''
        distance between cumulative degree difference distributions
        '''

        if len(degree_diff_a) > len(degree_diff_b):
            degree_diff_b += [1.0] * (len(degree_diff_a) - len(degree_diff_b))
        elif len(degree_diff_a) < len(degree_diff_b):
            degree_diff_a += [1.0] * (len(degree_diff_b) - len(degree_diff_a))
        
        return np.mean(np.abs(np.array(degree_diff_a) - np.array(degree_diff_b))) 


class Fit_degree_difference(object):
    '''
    Graph summary stat similarity without measuring the expensive to compute clustering coefficient
    '''
    def __init__(self, **kwargs):
        pass
    
    
    def graph_properties(self, G):
        '''
        compute the properties of the graph in question (nx graph object)
        
        max degree, ks_dist, degree_assortativity,
        '''
        degree = G.degree()
        
        degree_diff = [np.abs(degree[i] - degree[j]) for i,j in  G.edges() ] 
        c_degree_diff = [ 1.0/len(degree_diff) * len([k for k in degree_diff if k <= x] ) for x in range(max(degree_diff) + 1)]
        
        return c_degree_diff
    
    
    def graph_fitness(self, degree_diff_a, degree_diff_b):
        '''
        distance between cumulative degree difference distributions
        '''

        if len(degree_diff_a) > len(degree_diff_b):
            degree_diff_b += [1.0] * (len(degree_diff_a) - len(degree_diff_b))
        elif len(degree_diff_a) < len(degree_diff_b):
            degree_diff_a += [1.0] * (len(degree_diff_b) - len(degree_diff_a))
        
        
        return np.mean(np.abs(np.array(degree_diff_a) - np.array(degree_diff_b))) 