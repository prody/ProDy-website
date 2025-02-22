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

Latest version of ProDy_ (2.5.0) is required. This should be installed from GitHub (see manual_)
or conda-forge.


Recommended Programs
-------------------------------------------------------------------------------

Besides ProDy_, the Matplotlib_ library and VMD_ program are required for some steps in the tutorial. 
IPython_ is highly recommended for interactive usage.

Moreover, in the case of the lack of hydrogen atoms in protein structure,
additional package such as Openbabel_ or PDBfixer_ are required for
predicting hydrogen bonds.

.. _Openbabel: https://github.com/openbabel
.. _PDBfixer: https://github.com/openmm/pdbfixer
.. _manual: http://www.bahargroup.org/prody/manual/devel/develop.html


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
