.. CiGRAM documentation master file, created by
   sphinx-quickstart on Fri Oct  6 11:39:48 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CiGRAM's documentation!
==================================

.. toctree::
    generation
    :maxdepth: 2
    :caption: Contents:

Introduction
------------
CiGRAM is a Circular Gaussian Random Graph Generator written in python and C++.
The objective of cigram is to generate large, complex networks with realistic properties such as high clustering,
heterogenous degree distributions, community structure and assortativity.

These documents intend to describe how CiGRAM works, how you can generate networks with it and how you can use the
network fitting package within this project to match the desired properties of real world networks.
This documentation intends only to give a high level overview of the project.
For more details please see the initial publication of CiGRAM:

Installation
------------

Install with pip from pypi
.. code-block:: guess

    pip install cigram

Alternatively, install from sources:

.. code-block:: guess

    pip install -e .

It is reccomended that installation from source is done inside a virtualenv.


Testing
-------

To confirm that cigram works on your system and the build is functioning use the included tests.
These are written in py.test. To run simply install py.test (a requirement of CiGRAM and run).

.. code-block:: guess

    pytest

This should run without error.

Acknowledgements
----------------
None of this work would have been possible without the help of numerous talented people.
This includes the supervision of Jamie Twycross, Andrew Wood, Andrezj Bargeila, Natalio Krasnogor and Michael Holdsworth.
Much of this work was also supported by Paweł Widera.

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
