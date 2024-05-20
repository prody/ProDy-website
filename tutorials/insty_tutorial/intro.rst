Introduction
===============================================================================

This tutorial shows how to predict and display the interactions within
protein structure using the coordinates (single PDB) or ensemble of
conformations (multi-model PDB or dcd file gerenarted by NAMD_ program).

This module can be successfully used to study the distribution of different
types of interactions (``hydrogen bonds``, ``salt bridges``, ``repulsive ionic binding``,
``pi-cation``, ``pi-stacking``, ``hydrophobic interactions``, and ``disulfide bonds``) within
protein structure. It may help to distinguish the difference between wild type
protein and its mutant, indentify regions or simply residues which are privileged
to create a larger number of potential interactions in single PDB, ensemble PDB
(NMR data) or during the molecular dynamics simulation.   


Required Programs
-------------------------------------------------------------------------------

Latest version of ProDy_ is required.


Recommended Programs
-------------------------------------------------------------------------------

Besides ProDy_, the Matplotlib_ library and VMD_ program are required for
some steps in the tutorial. IPython_ is highly recommended for interactive usage.

To take full advantage of InSty's capabilities in predicting hydrophobic
interactions, we should download an additional component written in C++ and
Fortran. The additional component (hpb.so file) is available in the ProDy
repository (prody/proteins/hpbmodule folder) and should be chosen and
selected according to the version of Python that is being used by user.
This file should be copied to the ProDy folder (prody/proteins/) or to the
local directory.

Moreover, in the case of the lack of hydrogen atoms in protein structure,
additional package such as Openbabel_ or PDBfixer_ are required for
predicting hydrogen bonds.

.. _Openbabel: https://github.com/openbabel
.. _PDBfixer: https://github.com/openmm/pdbfixer


Getting Started
-------------------------------------------------------------------------------

To follow this tutorial, you will need the following files:

.. files.txt will be automatically generated

.. literalinclude:: files.txt


We recommend that you will follow this tutorial by typing commands in an
IPython session, e.g.::

  $ ipython

or with pylab environment::

  $ ipython --pylab


First, we will make necessary imports from ProDy and Matplotlib_ packages.

.. ipython:: python

   from prody import *
   from pylab import *
   import matplotlib
   ion()

We have included these imports in every part of the tutorial, so that
code copied from the online pages is complete. You do not need to repeat
imports in the same Python session.
