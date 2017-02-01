Introduction
===============================================================================

This tutorial describes how to calculate the mechanical stiffness map for a 
given protein structure, using the coordinates in the native state as *input*. 
The map provides a measure of the effective resistance of all pairs of residues
to increasing their inter-residue separation, in the context of the entire network.

This methodology can be successfully used for interpreting data from single 
molecule manipulation techniques results such as Single Molecule Force Spectroscopy 
(AFM) or optical tweezers. It may also be used for predicting the response to 
tension without performing expensive steered molecular dynamics simulations. 

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
