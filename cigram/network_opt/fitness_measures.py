from __future__ import division
import networkx as nx
import numpy as np
from networkx.algorithms.approximation.clustering_coefficient import average_clustering


def ks_distance(g_degree, m_degree, d_min=1, normed=True):
    """
    Calculate the Kolmogorov-Smirnov distance between two distrubtions

    This test is highly sensitive to the tail of the cumulative distrubtions
    For this reason d_min truncates the values (by default) to 10.

    Params:
    g_degree, m_degree should be degree sequences (e.g. g.degree().values() for networkx)
    """
    cum_vals_g = []
    cum_vals_m = []

    for k in range(int(max(g_degree + m_degree)) + 1):
        cum_vals_g.append(len([i for i in g_degree if i <= k]) / len(g_degree))
        cum_vals_m.append(len([i for i in m_degree if i <= k]) / len(m_degree))

    cum_vals_m = np.array(cum_vals_m[d_min:])
    cum_vals_g = np.array(cum_vals_g[d_min:])

    dists = np.abs(cum_vals_m - cum_vals_g)

    # make the distance uniformly sensitive across the range, useful for fat tails
    if normed:
        valid_d = np.where(dists > 0)
        # No point greater than 0, dists are identical
        if len(dists[valid_d]) == 0:
            return 0.0

        denom = np.sqrt(cum_vals_m[valid_d] * (1 - cum_vals_m[valid_d]))
        dists = dists[valid_d]
        valid_d = np.where(denom > 0)

        ds = dists[valid_d] / denom[valid_d]
        return np.max(ds)

    return np.max(dists)


class SummaryStatFitness(object):
    
    def __init__(self,
                 max_degree_fit_weight=0.5,
                 degree_fit_weight=1.0,
                 clustering_fit_weight=0.5,
                 assort_fit_weight=0.5):
        """
        Uses distance between different summary statistics to compute graph similiarity
        :param max_degree_fit_weight:
        :param degree_fit_weight:
        :param clustering_fit_weight:
        :param assort_fit_weight:
        """
        self.mdfw = max_degree_fit_weight
        self.dfw = degree_fit_weight
        self.cfw = clustering_fit_weight
        self.afw = assort_fit_weight

    @staticmethod
    def graph_properties(g):
        """
        compute the properties of the graph in question
        
        max degree, ks_dist, degree_assortativity, clustering coefficient
        """
        degree = g.degree().values()
        props = dict(
                max_degree=max(degree),
                degree=degree,
                clustering=average_clustering(g, trials=1000),
                assort=nx.degree_assortativity_coefficient(g)
        )

        return props

    @staticmethod
    def graph_fitness(prop_a, prop_b):
        """
        Summary statistics dissimmilarity
        """
        max_deg_dist = np.abs(np.log10(prop_a["max_degree"]) - np.log10(prop_b["max_degree"]))
        ks_dist = ks_distance(prop_a["degree"], prop_b["degree"], normed=False, d_min=10)
        assort_dist = np.abs(prop_a["assort"] - prop_b["assort"]) 
        clust_dist = np.abs(prop_a["clustering"] - prop_b["clustering"])

        return ks_dist + 0.5 * (max_deg_dist + assort_dist + clust_dist)


class NoClusterSummaryFitness(object):

    def __init__(self,
                 max_degree_fit_weight=0.5,
                 degree_fit_weight=1.0,
                 assort_fit_weight=0.5):
        """
        Graph summary stat similarity without measuring the expensive to compute clustering coefficient
        """
        self.mdfw = max_degree_fit_weight
        self.dfw = degree_fit_weight
        self.afw = assort_fit_weight

    @staticmethod
    def graph_properties(g):
        """
        compute the properties of the graph in question (nx graph object)
        
        max degree, ks_dist, degree_assortativity,
        """
        degree = g.degree().values()
        props = dict(
            degree=degree,
            max_degree=max(degree),
            assort=nx.degree_assortativity_coefficient(g)
        )

        return props

    @staticmethod
    def graph_fitness(prop_a, prop_b):
        """
        Summary statistics dissimmilarity
        """
        max_deg_dist = np.abs(np.log10(prop_a["max_degree"]) - np.log10(prop_b["max_degree"]))
        ks_dist = ks_distance(prop_a["degree"], prop_b["degree"], normed=False, d_min=10)
        assort_dist = np.abs(prop_a["assort"] - prop_b["assort"]) 
        
        return ks_dist + 0.5 * (max_deg_dist + assort_dist)


class DegreeDifference(object):
    """
    Graph summary stat similarity without measuring the expensive to compute clustering coefficient
    """
    def __init__(self, **kwargs):
        pass

    @staticmethod
    def graph_properties(g):
        """
        compute the properties of the graph in question (nx graph object)
        
        max degree, ks_dist, degree_assortativity,
        """
        degree = g.degree()
        
        degree_diff = [np.abs(degree[i] - degree[j]) for i, j in g.edges()]
        c_degree_diff = [1.0/len(degree_diff) *
                         len([k for k in degree_diff if k <= x]) for x in range(max(degree_diff) + 1)]

        return c_degree_diff

    @staticmethod
    def graph_fitness(degree_diff_a, degree_diff_b):
        """
        distance between cumulative degree difference distributions
        """

        if len(degree_diff_a) > len(degree_diff_b):
            degree_diff_b += [1.0] * (len(degree_diff_a) - len(degree_diff_b))
        elif len(degree_diff_a) < len(degree_diff_b):
            degree_diff_a += [1.0] * (len(degree_diff_b) - len(degree_diff_a))

        return np.mean(np.abs(np.array(degree_diff_a) - np.array(degree_diff_b))) 