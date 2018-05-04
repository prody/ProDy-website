Introduction
===============================================================================



This tutorial demonstrates how to use perturbation response scanning (PRS) to determine 
sensors and effectors, which are important for allosteric signal transduction. The PRS 
approach is derived from linear response theory where perturbation forces are applied 
via a covariance matrix, which can be derived from elastic network models or MD simulations.

The example used in this tutorial is the AMPA-type ionotropic glutamate receptor (AMPAR; PDB 3kg2), 
which we recently studied using this method (see Figure 6 of [DB15]_).

The theory was originally described in [AA09]_ and extended to include sensors in [GB14]_.

Required Programs
-------------------------------------------------------------------------------

The latest version of ProDy_ (v 1.10.1) is required along with Numpy_ and Matplotlib_. 
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

.. [AA09] Atilgan C, Atilgan AR. 
   Perturbation-response scanning reveals ligand entry-exit mechanisms of ferric 
   binding protein. *PLoS Comput. Biol.* **2009** 5:e1000544.
.. [GB14] General IJ, Liu Y, Blackburn ME, Mao W, Gierasch LM, Bahar I.
   ATPase subdomain IA is a mediator of interdomain allostery in Hsp70 molecular 
   chaperones. *PLoS Comput. Biol.* **2014** 10:e1003624.
.. [DB15] Dutta A, Krieger J, Lee JY, Garcia-Nafria J, Greger IH, Bahar I.
   Cooperative Dynamics of Intact AMPA and NMDA Glutamate Receptors: Similarities 
   and Subfamily-Specific Differences. *Structure* **2015** 23:1692-1704