import networkx as nx
from cigram.network_opt.cigram_optimiser import CigramOptimiser
from cigram.network_opt.model_parameters import optimise_model
from cigram.network_opt.fitness_measures import NetLSDFit


def test_fitter():
    """
    General test
    :return:
    """
    graph = nx.karate_club_graph()
    mdl = CigramOptimiser(graph, NetLSDFit, k_bounds=2)
    optimise_model(mdl)
