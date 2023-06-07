Introduction
===============================================================================

This tutorial shows how to predict and display the interactions within
protein structure, protein-protein, and protein-ligand structure, using the
coordinates or trajectories (dcd files).

This module can be successfully used to study the distribution of different
types of interactions (hydrogen bonds, salt bridges, repulsive ionic bindnig,
pi-cation, pi-stacking, and hydrophobic) within protein structure for a
single PDB file as well as for the trajectory. It may help to distinguish the
difference between protein wt and mutant, stability of protein interactions
in the course of the simulation, and to indentify regions or simply residues
which residues that are privileged to create a larger number of potential
interactions or energetically more stable.   


Required Programs
-------------------------------------------------------------------------------

Latest version of ProDy_ is required.


Recommended Programs
-------------------------------------------------------------------------------

Besides ProDy_, the Matplotlib_ library and VMD_ program are required for
some steps in the tutorial. IPython_ is highly recommended for interactive usage.
Moreover, in the case of predicting protein-ligand interactions or in the case of 
the lack of hydrogen atoms in protein structure, PLIP_ and Openbabel_ or 
PDBfixer_ are required, respectively.

.. _PLIP: https://github.com/pharmai/plip
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


First, we will make necessary imports from ProDy and Matplotlib
packages.

.. ipython:: python

   from prody import *
   from pylab import *
   import matplotlib
   ion()

We have included these imports in every part of the tutorial, so that
code copied from the online pages is complete. You do not need to repeat
imports in the same Python session.
