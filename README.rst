CiGRAM - Circular Gaussian Random gRAph Model
#############################################

|Build Status| |Coverage Status| |PyPI| |docs|

CiGRAM is a method for generating large complex networks with heterogeneous properties observed in real-world networks such as heavy tailed node degree distributions,  assortative connectivity and community structure.
Currently the model only supports unweighted, undirected graphs.

For a demonstration and brief definition of some of the features of cigram at glance please see the web visualisation developed with help from Paweł Widera at from the School of Computing Science at Newcastle University:

http://cigram.ico2s.org

This work was developed as part of a PhD. If using this work in publications please cite the PhD thesis which also
includes a comprehensive description of the network, its features and some applications (Chapters 4-5).

Gilbert, J. P. "A probabilistic model for the evaluation of module extraction algorithms in complex biological networks." PhD disseration., University of Nottingham, 2015.

Bibtex entry:

.. code-block::

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

Requires a C++ complier, libboost and python (versions 2.7, 3.5, 3.6 tested).

Install libbosst from http://www.boost.org/

Should work with libboost version 1.58 to 1.66. Note that the libboost version in ubuntu 16.04 repositories (1.55)
does not work.

Clone repository and use pip (virtualenv is suggested):

.. code-block:: bash

    pip install cigram

Alternatively, the package developed here can be installed by cloning the git repository and running the command:

.. code-block:: bash

    python setup.py install

Usage
-----

The current version of the model can be used to generate graphs with the wrapper to the C++ model

To generate a random model with 1000 nodes, average degree of 5 split in to 5 communities

.. code-block:: python

    import cigram
    graph, vertex_positions, community_memberships = cigram.cigram_graph(1000, 5, 5)

graph is a networkx model, vertex_positions is a dict of vertex positions by node id, and community memberships is a dict of vertex memberships of each cluster.


.. |Build Status| image:: https://travis-ci.org/azimov/cigram.svg?branch=master
   :target: https://travis-ci.org/azimov/cigram
.. |Coverage Status| image:: https://codecov.io/github/azimov/cigram/coverage.svg?branch=master
   :target: https://codecov.io/github/azimov/cigram
.. |Build status2| image:: https://ci.appveyor.com/api/projects/status/
   :target: https://ci.appveyor.com/project/azimov/cigram/branch/master
.. |PyPI| image:: https://badge.fury.io/py/cigram.svg
   :target: https://pypi.python.org/pypi/cigram
.. |docs| image:: https://readthedocs.org/projects/cigram/badge/?style=flat
    :target: https://readthedocs.org/projects/cigram
    :alt: Documentation Status
