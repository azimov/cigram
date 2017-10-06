Generating networks
===================

Generating networks in CiGRAM is straightforward.

.. code-block:: python

    from cigram import generate_graph
    n = 1000
    avg_k = 4.95
    k = 1
    graph, positions, communities = generate_graph(n, avg_k, k)

Here n is the number of nodes, avg_k is the desired average degree and k is the number of communities.
Note that the number of communities here is fixed at 1, so the graph will not generate artificial clusters.

The resulting returned tuple of ``graph``, ``positions`` and ``communities`` is a networkx graph object, a dictionary
of points upon a unit circle for each node, and a dictionary for the community membership of each vertex.

More complex parameters allow you to generate different heterogenous degree distributions.
This is controlled by the parameters sigma_f and sigma_r which have an effect on the underlying probability space for
the connections between nodes.


.. code-block:: python

    sigma_r = 0.8
    sigma_f = 0.8
     graph, positions, communities = generate_graph(n, avg_k, k)


To generate networks with assortativity this can be specified with the parameter a (by default this is 0).