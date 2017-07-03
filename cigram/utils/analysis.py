from __future__ import division
from community import modularity as _modularity
from community import best_partition
import numpy as np
from collections import defaultdict
import networkx as nx
from numpy import choice


def modularity(g, part=None):
    """
    Return the newman modularity of graph g
    Optional parameter part, is a partition of a graph
    If part is None, this will be maximised with the louvian approach
    """
    if part is None:
        part = best_partition(g)

    return _modularity(part, g)


def partion_convert(part):
    """
    Swap a partition dict format from community:[nodes] to node:community_i
    """
    n_k = {}
    for k, nodes in part.items():
        for node in nodes:
            n_k[node] = k
            
    return n_k


def partion_groups(part):
    """
    Swap a partition dict format from node:community_i to community:[nodes]
    """
    partition = defaultdict(list)
    for node, group in part.items():
        partition[group].append(node)
            
    return partition


def gini_coefficient(data):
    """ 
    The gini coefficient measures the level of inequality in a distribution;
    0 is uniform equality, 1 is perfect inequality
    """
    
    # Values must be ranked
    data = sorted(data.values())
    n = len(data)
    
    nume = 0.0
    denome = 0.0
    for i in range(1, len(data) + 1):
        nume += data[i - 1] * i
        denome += data[i - 1]
        
    return ((nume * 2.0) / (n * denome) ) - ((n+1.0)/n)


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
        cum_vals_g.append(len([i for i in g_degree if i <= k])/len(g_degree))
        cum_vals_m.append(len([i for i in m_degree if i <= k])/len(m_degree))

    cum_vals_m = np.array(cum_vals_m[d_min:])
    cum_vals_g = np.array(cum_vals_g[d_min:])
    
    dists = np.abs(cum_vals_m - cum_vals_g)

    # make the distance uniformly sensitive across the range, useful for fat tails
    if normed:
        valid_d = np.where(dists > 0)
        # No point greater than 0, dists are identical
        if len(dists[valid_d]) == 0:
            return 0.0
        
        denom = np.sqrt(cum_vals_m[valid_d] * (1- cum_vals_m[valid_d]))
        dists = dists[valid_d]
        valid_d = np.where(denom > 0)

        ds =  dists[valid_d]/denom[valid_d]
        return np.max(ds)
    
    return np.max(dists)


def sample_model_function(m, m_args, measure_funcs, samples=100, seed=None):
    """
    Repeatedly run a model with a set of parameters and measure functions
    """
    
    results = defaultdict(list)
    
    for i in range(samples):
        m_args["seed"] = seed
        model_graph = m(**m_args)
        
        for measure, func in measure_funcs.items():
            results[measure].append(func(model_graph))
        
        if seed is not None:
            seed +=1
        
    return results


def result_summary(result):
    """
    returns Mean, standard deviation, range for an iterable
    """
    return np.mean(result), np.std(result), [min(result), max(result)]


def sample_shortest_path_length(g, sample_size=1000, seed=None):
    if seed is not None:
        np.random.seed(seed)
    
    # select sample size of node_pairs from n choose 2 possible combinations
    avg = 0.0
    selected=set()
    failcount = 0
    while len(selected) < sample_size and failcount < sample_size/2:
        pair = choice(g.nodes(), 2)
        key = tuple(sorted(pair)) 
        if key not in selected:
            selected.add(key)
            try:
                avg += nx.shortest_path_length(g,source=pair[0],target=pair[1])
            except nx.NetworkXNoPath:
                failcount += 1
                selected.remove(key)
    return avg / sample_size


def joint_degree_matrix(g):
    """
    generate the joint degree distribution for each edge in a network
    :param g: networkx graph
    :return: numpy matrix of joint degrees for each edge in g
    """
    deg = g.degree()
    kmax = max(deg.values()) + 1
    g_d = np.zeros((kmax, kmax), dtype=np.int)
    
    for i,j in g.edges():
        g_d[deg[i],deg[j]] += 1
        g_d[deg[j],deg[i]] += 1
        
    return g_d


def cumulative_joint_degree_matrix(g):
    """
    Returns the cumulative joint degree distribution of a network

    :param g: networkx graph object
    :return: matrix of cumulative joint degree distribution
    """
    deg = g.degree()
    
    kmax = max(deg.values()) + 2
    jdegs = [(deg[i], deg[j]) for i,j in g.edges()]
    
    m = np.zeros((kmax, kmax))
    for k in range(kmax):
        for l in range(kmax):
            m[k,l] = len([ _ for di, dj in jdegs if (di >= k and dj >=l) or (di >= l and dj >=k)])

    return m/g.number_of_edges()


def jdeg_dist(g, gp):
    """
    L2 norm of cumulative joint degree matricies
    """
    deg_g = g.degree()
    deg_gp = gp.degree()
    
    kmax = max(deg_gp.values() + deg_g.values())
    
    gdegs = [(deg[i], deg[j]) for i,j in g.edges()]
    pdegs = [(deg[i], deg[j]) for i,j in gp.edges()]
    
    m = np.zeros((kmax, kmax))
    mp = np.zeros((kmax, kmax))
    for k in range(kmax):
        for l in range(kmax):
            m[k,l] = len([(di, dj) for di, dj in gdegs if (di >= k and dj >=l) or (di >= l and dj >=k)]) / g.number_of_edges()
            mp[k,l] = len([(di, dj) for di, dj in pdegs if (di >= k and dj >=l) or (di >= l and dj >=k)]) / g.number_of_edges()

    return np.max(m - mp)