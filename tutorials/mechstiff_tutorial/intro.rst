Introduction
===============================================================================

This tutorial shows how to calculate mechanical stiffness of protein structure
based on the coarse-grained elastic network models with the coordinates of the
native state as an *input*. Based on the topology of interresidue contacts 
in the native state this method can show a fundamental property that 
dominates in the observed behavior of PDB structure.

This methodology can be successfully used for interpretation of single molecule 
manipulation techniques results such as Single Molecule Frorce Spectroscopy 
(AFM) or optical tweezers, in prediction of applied external force location
without performance of expensive steered molecular dynamics simulations.
It can be worn as a support in defining the direction of pulling in Steered MD
simulations as well.

Theory and example of usage has been described in [EB08]_.

Required Programs
-------------------------------------------------------------------------------

Latest version of ProDy_ us required. Matplotlib_ library and VMD_ progran is 
required for some options. IPython_ is highly recommended for interactive usage.


Getting Started
-------------------------------------------------------------------------------

To follow this tutorial, you will not need any additional files.

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
   ion()

We have included these imports in every part of the tutorial, so that
code copied from the online pages is complete. You do not need to repeat
imports in the same Python session.


How to Cite
-------------------------------------------------------------------------------

If you benefited from Mechanical Stiffness Calculations in your research, 
please cite the following paper:

.. [EB08] Eyal E., Bahar I. Toward a Molecular Understanding of 
   the Anisotropic Response of Proteins to External Forces: Insights from 
   Elastic Network Models. *Biophys. J.* **2008** 94:3424-34355.
