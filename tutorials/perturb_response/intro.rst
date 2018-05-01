Introduction
===============================================================================

This tutorial describes how to perform perturbation response scanning analysis for a 
given protein structure, using the coordinates in the native state as *input*. 
This provides a measure of the effectiveness and sensitivity of residues
in propagation of allosteric signals through the protein.

The theory has been extensively described. See [AA09]_ for the original description 
and [GB14]_ for a description of our version.

Required Programs
-------------------------------------------------------------------------------

The latest version of ProDy_ (v 1.10) is required along with Numpy_ and Matplotlib_. 
IPython_ is highly recommended for interactive usage.


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

If you benefited from Perturbation Response Analysis in your research, 
please cite the following papers:

.. [AA09] Atilgan C, Atilgan AR. Perturbation-response scanning reveals 
   ligand entry-exit mechanisms of ferric binding protein. *PLoS Comput. Biol.* 
   **2009** 5:e1000544.
.. [GB14] General IJ, Liu Y, Blackburn ME, Mao W, Gierasch LM, Bahar I.
   ATPase subdomain IA is a mediator of interdomain allostery in Hsp70 molecular 
   chaperones. *PLoS Comput. Biol.* **2014** 10:e1003624.