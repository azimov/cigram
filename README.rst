CiGRAM - Circular Gaussian Random gRAph Model
#############################################

CiGRAM is a method for generating large complex networks with heterogeneous properties observed in real-world networks such as heavy tailed node degree distributions,  assortative connectivity and community structure.
Currently the model only supports unweighted, undirected graphs.

For a demonstration and brief definition of some of the features of cigram at glance please see the web visualisation developed with help from Paweł Widera at from the School of Computing Science at Newcastle University:

http://cigram.ico2s.org

This work was developed as part of a PhD. If using this work in publications please cite the PhD thesis which also
includes a comprehensive description of the network, its features and some applications (Chapters 4-5).

Gilbert, J. P. "A probabilistic model for the evaluation of module extraction algorithms in complex biological networks." PhD disseration., University of Nottingham, 2015.

Bibtex entry:

.. code-block:: guess

    @phdthesis{gilbert2015probabilistic,
      title={A probabilistic model for the evaluation of module extraction algorithms in complex biological networks},
      author={Gilbert, JP},
      year={2015},
      school={University of Nottingham}
    }

This library also includes Lancichinetti–Fortunato–Radicchi (LFR) binary benchmark graphs.
Future versions plan to include all the LFR benchmarks (weighted, directed and hierarchical networks).

Installation
------------

Requires a C++ complier, libboost and python (versions 2.7+ and 3.5+ tested).

Clone repository and use pip:

.. code-block:: guess

    pip install -r requirements.txt
    pip install .

Alternatively, the package developed here can be installed by cloning the git repository and running the command:

.. code-block:: guess

    python setup.py install

Usage
-----

The current version of the model can be used to generate graphs with the wrapper to the C++ model

To generate a random model with 1000 nodes, average degree of 5 split in to 5 communities

.. code-block:: python

    import cigram
    graph, vertex_positions, community_memberships = cigram.cigram_graph(1000, 5, 5)

graph is a networkx model, vertex_positions is a dict of vertex positions by node id, and community memberships is a dict of vertex memberships of each cluster.
