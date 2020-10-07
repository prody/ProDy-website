Introduction
===============================================================================

This tutorial describes how to calculate signature dynamics for a family of
proteins with similar structures using Elastic Network Models. This method creates 
an ensemble of aligned structures and calculates statistics such as means and 
standard deviations on various dynamic properties including mode profiles, 
mean square fluctuations and cross-correlation matrices. It also includes tools 
for classifying family members based on their sequence, structure and dynamics.

The theory and usage of this toolkit will be described in [SZ18]_.

Required Programs
-------------------------------------------------------------------------------

The latest version of ProDy_ is recommended along with NumPy_ and Matplotlib_. 
IPython_ is highly recommended for interactive usage.


Getting Started
-------------------------------------------------------------------------------

This tutorial contains three parts. In the first part, we will have a quick 
walk-through on the SignDy_ calculations and functions using the example of type-I 
periplasmic binding protein (PBP-I) domains, in which case the data is convienient 
collected from the `Dali server`_ [LH10]_. The second part will be a more detailed 
+periplasmic binding protein (PBP-I) domains, in which case the data is convieniently 
+collected from the `Dali server`_ [LH10]_. The second part will be review how to use 
+the CATH database to build the ensemble. The third part will be a more detailed 
tutorial on building an ensemble 'manually' from scratch, and try to reproduce the 
figures presented in [SZ18]_.

We recommend that you will follow this tutorial by typing commands in an
IPython session, e.g.::

  $ ipython

or with pylab environment::

  $ ipython --pylab


First, we will make necessary imports from the ProDy_, NumPy_ and Matplotlib_
packages.

.. ipython:: python

    from prody import *
    from numpy import *
    from matplotlib.pyplot import *
    ion()

We have included these imports in every part of the tutorial, so that
code copied from the online pages is complete. You do not need to repeat
imports in the same Python session.


How to Cite
-------------------------------------------------------------------------------

If you benefited from SignDy Calculations in your research, 
please cite the following paper:

.. [SZ18] Zhang S, Li H, Krieger J, Bahar I. 
    Shared signature dynamics tempered by local fluctuations enables fold adaptability and specificity.
    *Mol. Biol. Evol.* **2019** (In Press)

.. [LH10] Holm L, Rosenström P.
    Dali server: conservation mapping in 3D.
    *Nucleic Acids Res.* **2010** 10(38):W545-9

.. _`Dali server`: http://ekhidna2.biocenter.helsinki.fi/dali/

.. _`SignDy`: http://prody.csb.pitt.edu/test_prody/tutorials/signdy_tutorial/
