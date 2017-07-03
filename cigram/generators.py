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
import time
import networkx

from cigram.cmodel import generate_graph


def community_graph(n, density, k, ek_per=0.01, p_o=0, a=0.0, sigma_nodes=1.0, sigma_edges=1.0, community_sigma_r=None,
                    community_sigma_f=None, connected=True, seed=None, ret_com_pos=False, min_degree=1,
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

    edges, communities, node_positions, community_positions = generate_graph(int(n), int(k), float(density),
                                                                             float(sigma_nodes), float(sigma_edges),
                                                                             float(a), float(community_sigma_r),
                                                                             float(community_sigma_f), float(ek_per),
                                                                             float(p_o), int(conn), int(min_degree),
                                                                             int(min_community_size), seed)
    graph = networkx.Graph()
    graph.add_nodes_from(range(int(N)))
    graph.add_edges_from(edges)

    communities = dict(enumerate(communities))
    node_positions = dict(enumerate(node_positions))
    community_positions = dict(enumerate(community_positions))

    if ret_com_pos:
        return graph, node_positions, communities, community_positions

    return graph, node_positions, communities


def single_process_graph(n, density, a=0, sigma_nodes=1, sigma_edges=1, connected=True, min_degree=1, seed=None):
    """
    This is the simplest graph model, it generates no communities.
    Generates graphs with a tunable level of assortativity by modfiying parameter a.

    """
    if seed is None:
        seed = int(time.time())

    conn = 0
    if connected:
        conn = 1

    edges, _, pos, _ = generate_graph(n, 1, density, sigma_nodes, sigma_edges, a, 1.0, 1.0, 0.0, 0.0, conn, min_degree,
                                      seed)

    graph = networkx.Graph()
    graph.add_nodes_from(range(int(n)))
    graph.add_edges_from(edges)

    return graph, pos