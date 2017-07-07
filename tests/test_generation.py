from __future__ import division, print_function
from cigram import community_graph, single_process_graph
import pytest
import networkx as nx


def test_generation():
    """Test that the cmodel functions"""
    graph, positions, communities = community_graph(1000, 4.965, 3)

    assert graph.number_of_nodes() == 1000
    assert graph.number_of_edges() == 4965


def test_large_graph():
    """Test that larger graphs don't freeze. Should take a while to run"""
    graph, positions, communities = community_graph(20000, 4.965, 3, sigma_edges=0.8, sigma_nodes=0.8)

    assert graph.number_of_nodes() == 20000
    assert graph.number_of_edges() == 99300


def test_sp_generation():
    """Test for graphs without community structure"""
    graph, positions = single_process_graph(1000, 4.965)
    assert graph.number_of_nodes() == 1000
    assert graph.number_of_edges() == 4965


@pytest.mark.parametrize("params,test_funcs", [
    ({"a": -5}, {"assort": (lambda g, p, c: nx.degree_assortativity_coefficient(g), lambda r: r < -0.2)}),
    ({"a": 5}, {"assort": (lambda g, p, c: nx.degree_assortativity_coefficient(g), lambda r: r > 0.2)}),
    ({"sigma_nodes": 0.5, 'sigma_edges': 0.5}, {}),
    ({"sigma_nodes": 5.0, 'sigma_edges': 5.0}, {}),
    ({"density": 1.0}, {"desnsity": (lambda g, p, c: nx.density(g), lambda r: r == 1.0)}),
    ({"community_sigma_f": 0.5, "community_sigma_r": 0.5}, {}),
    ({"community_sigma_f": 5.0, "community_sigma_r": 5.0}, {})
])
def test_param_bound(params, test_funcs):
    """
    Test the boundaries of the generator to ensure functioning behaviour
    The pytest mark format expects the specified model parameters dictionary,
    then a dict of tuples (functions to run, assertion for the results of the functions)
    """
    unmodified_params = params.copy()
    params['seed'] = 1337

    if 'n' not in params:
        params['n'] = 2000

    if 'avg_deg' not in params:
        params['avg_deg'] = 5

    if 'k' not in params:
        params['k'] = 3

    graph, positions, communities = community_graph(**params)
    for k, (func, res_t) in test_funcs.items():
        # print statements are used for debug
        print(unmodified_params)
        print(k)
        fr = func(graph, positions, communities)
        print(fr)
        assert res_t(fr)
