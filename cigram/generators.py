"""
CiGRAM - Gaussian Random grAph Model

Copyright (C) 2017  James Gilbert

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
more details: http://www.gnu.org/licenses/gpl.html
"""
from __future__ import division
import time
import networkx

import cigram.cmodel
import cigram.lfr_model
from collections import defaultdict

import logging

UNLIKELY = -214741

FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT, datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('cigram.generators')
logger.setLevel(logging.DEBUG)


def cigram_graph(n, avg_deg, k,
                 density=None, ek_per=0.01, p_o=0, a=0.0,
                 sigma_nodes=1.0, sigma_edges=1.0,
                 community_sigma_r=None, community_sigma_f=None, connected=True,
                 seed=None, ret_com_pos=False, min_degree=1,
                 min_community_size=0):
    """
    Generate a graph with a community structure

    ek_per relates to $E_k$ in the paper; the fraction of edges between communities
    We also allow for overlapping nodes with p_o but the model is a bit shakey for edges between nodes
    in the overlapping group so needs improvements
    The overlapping model can easily end up in deadlocks where the probability of assigning valid edges is so low they
    don't become assigned

    The probability distributions pos_prob is a random variate, the result should map directly to a value from edge_prob
    that gives nodes their connection probability
    """
    if seed is None:
        seed = int(time.time())

    conn = 0
    if connected:
        conn = 1

    if community_sigma_r is None:
        community_sigma_r = sigma_nodes

    if community_sigma_f is None:
        community_sigma_f = sigma_edges

    # Density overrides avg_k - weird behaviour but designed for some legacy code
    if density is None:
        density = 2 * avg_deg * 1/(n-1)

    params = (
        int(n),
        int(k),
        float(density),
        float(sigma_nodes),
        float(sigma_edges),
        float(a),
        float(community_sigma_r),
        float(community_sigma_f),
        float(ek_per),
        float(p_o),
        int(conn),
        int(min_degree),
        int(min_community_size),
        seed
    )

    edges, communities, node_positions, community_positions = cigram.cmodel.generate_graph(*params)
    graph = networkx.Graph()
    graph.add_nodes_from(range(int(n)))
    graph.add_edges_from(edges)

    communities = dict(enumerate(communities))
    node_positions = dict(enumerate(node_positions))
    community_positions = dict(enumerate(community_positions))

    graph.name = "cigram_" + "-".join([str(p) for p in params])

    if ret_com_pos:
        return graph, node_positions, communities, community_positions

    return graph, node_positions, communities


def single_process_graph(n, avg_deg,
                         density=None, a=0, sigma_nodes=1, sigma_edges=1, connected=True, min_degree=1, seed=None):
    """
    This is the simplest graph model, it generates no communities.
    Generates graphs with a tunable level of assortativity by modfiying parameter a.
    """
    if seed is None:
        seed = int(time.time())

    if density is None:
        density = 2 * avg_deg * 1/(n-1)

    conn = 0
    if connected:
        conn = 1

    params = (
        int(n),
        1,
        float(density),
        int(sigma_nodes),
        int(sigma_edges),
        int(a),
        1.0,
        1.0,
        0.0,
        0.0,
        bool(conn),
        int(min_degree),
        0,
        int(seed)
    )

    edges, _, pos, _ = cigram.cmodel.generate_graph(*params)

    graph = networkx.Graph()
    graph.add_nodes_from(range(int(n)))
    graph.add_edges_from(edges)

    graph.name = "cigram_" + "-".join([str(p) for p in params])
    return graph, pos


def lfr_benchmark_graph(n, average_degree, max_degree, mu, tau=2.0, tau2=1.0, minc_size=3, maxc_size=None,
                        overlapping_nodes=0,
                        overlapping_memberships=0, clustering=UNLIKELY,
                        rand=False, sup=False, inf=False,
                        nodewise_membership=False, seed=None):
    """
    Calls external C library implementation of the LFR benchmark model graph

    :param n: number of nodes
    :param average_degree: average degree of nodes
    :param max_degree: maximum degree of nodes
    :param mu: community mixing parameter
    :param tau: degree power law exponent
    :param tau2: community size power law exponent
    :param minc_size: minimum community size
    :param maxc_size: maximum community size
    :param overlapping_nodes: number of overlapping nodes
    :param overlapping_membership:
    :param clustering: Fixed clustering coefficient NOTE - not properly test. Can hang
    :param rand: generate a random network TODO: implement this
    :param sup: #TODO
    :param inf: #TODO
    :return: graph, memberships - networkx.Graph and dict of community membership
    """
    if sup and inf:
        raise ValueError("sup and inf cannot both be true.")

    excess = False
    defect = False
    fixed_range = False

    if maxc_size is None:
        maxc_size = n

    if seed is None:
        seed = int(time.time())

    # Order and types really matter here
    params = (
        int(n),
        float(average_degree),
        float(tau),
        float(tau2),
        float(mu),
        float(clustering),
        int(max_degree),
        int(overlapping_nodes),
        int(overlapping_memberships),
        int(minc_size),
        int(maxc_size),
        int(seed),
        bool(excess),
        bool(defect),
        bool(fixed_range),
    )

    # Call to external C function
    logger.debug("Calling LFR Code")
    edges, communities = cigram.lfr_model.generate_graph(*params)
    logger.debug("Exit LFR Code")

    graph = networkx.Graph()
    graph.add_edges_from(edges)
    memberships = defaultdict(list)

    for v, c in communities:
        if nodewise_membership:
            memberships[v].append(c)
        else:
            memberships[c].append(v)

    graph.name = "lfr_binary_" + "-".join([str(p) for p in params])

    return graph, memberships