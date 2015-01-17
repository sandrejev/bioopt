=============
Installation
=============

Currently BioOpt package requires Python 2.7 or later. Python 2.7 is available on many linux distributions and OS X by default.
If you are using Windows we recommend 32bit `Python distribution from ActiveState <http://www.activestate.com/activepython/downloads>`_.
At this moment Python 3 is not supported and we will look into the matter when python 3 will become more widespread.

Python installation
====================

Mac OS X
---------
0. [install pip] (http://pip.readthedocs.org/en/latest/installing.html).
1. In a terminal, run ``sudo pip install bioopt``

GNU/Linux
----------
0. [install pip] (http://pip.readthedocs.org/en/latest/installing.html).
1. You will need python C headers. On debian-based systems (including Ubuntu and Mint) this can be done by running
   ``sudo apt-get install python-dev``
2. Afterwards run ``sudo pip install cobra`` in terminal

On Microsoft Windows
---------------------
0. Download and install `Python <http://www.activestate.com/activepython/downloads>`_
1. Open PyPM package manager which comes with ActiveState Python distribution and run `pypm install bioopt`

Manual installation (system independent)
-----------------------------------------
1. Clone the git repository ``git clone https://github.com/sandrejev/bioopt.git``
2. Run **setup.py** script which ships with the package: ``python setup.py``

Installing optional dependencies
=================================
libSBML
--------
libSBML is a library for reading and writing `SBML <http://sbml.org/Main_Page>`_ files. BioOpt package provides
converters to libSBML objects. libSBML requirement is optional. If you don't need to convert anything then there is no
need to install libSBML.

To install libSBML o systems with python package manager run ``sudo pip install libsbml`` or ``pypm install libsbml``
(on Windows).

COBRApy
--------
`COBRAPy <https://github.com/opencobra/cobrapy>`_ is Python port of `COBRA toolbox <http://opencobra.github.io/>`_.
COBRApy provides an API to manipulate metabolic models and perform different sort of analysis on resulting models
(like FBA or OptKnock). BioOpt package uses COBRApy to perform some of the simulations. If you are performing any
simulations most likely you will need to install COBRApy

Instructions on how to install COBRApy can be found
`here <https://github.com/opencobra/cobrapy/blob/master/INSTALL.md>`_

Pandas
-------
`Pandas <http://pandas.pydata.org/>`_ is a data analysis library for Python. It is heavily inspired by R and provides
similar classes and methods. Pandas dependency is optional and only needed if you are using legacy module ``output_parser``.
Output parser support and requirement for Pandas library will be dropped in near future.

Solvers
========
By default BioOpt uses GLPK solver which comes with COBRApy package but other solvers will be supported as well.

CPLEX
------
CPLEX is a very fast optimization package from IBM and it comes for free for
`academic use <https://www.ibm.com/developerworks/university/academicinitiative/>`_. CPLEX support both linear and
quadratic optimizations. To download CPLEX for academia follow
`this link <https://www.ibm.com/ibm/university/academic/pub/jsps/assetredirector.jsp?asset_id=1070>`_

Testing your installation
==========================

1. Start python
2. Type the following into the Python shell

.. code:: python

    from bioopt.test import *
    test_all()


