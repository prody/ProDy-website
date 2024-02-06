Introduction
===============================================================================

This tutorial shows how to use ProDy within Scipion including experimental ensemble analysis,
NMA and ClustENMD simulations (see Figure 3 of [KJ23]_).


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

.. [KJ23] Krieger JM, Sorzano COS, Carazo JM.
   Scipion-EM-ProDy: A Graphical Interface for the ProDy Python Package within the Scipion Workflow Engine 
   Enabling Integration of Databases, Simulations and Cryo-Electron Microscopy Image Processing. *Int. J. Mol. Sci.* **2023** 24(18):14245.