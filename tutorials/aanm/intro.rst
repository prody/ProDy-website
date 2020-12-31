Introduction
===============================================================================

This tutorial shows how to transition between two conformational states using 
adaptive ANM. The example used in this tutorial is the GroEL chaperonin in the 
R'' and T states (chain A of PDB structures 1GRU and 1GR5), which we previously 
studied using this method (see Figure 4 of [ZY09]_).


Required Programs
-------------------------------------------------------------------------------

The latest version of ProDy_ is required.

Getting Started
-------------------------------------------------------------------------------

We recommend that you will follow this tutorial by typing commands in an
IPython session, e.g.::

  $ ipython


First, we will make necessary imports from ProDy and other
packages.

.. ipython:: python

   from prody import *
   from pylab import *
   ion()

We have included these imports in every part of the tutorial, so that
code copied from the online pages is complete. You do not need to repeat
imports in the same Python session.

.. [ZY09] Yang Z, Májek P, Bahar I.
   Allosteric transitions of supramolecular systems explored by network models: 
   application to chaperonin GroEL. *PLoS Comput Biol.* **2009** 5(4):e1000360. 
   