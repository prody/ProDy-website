Introduction
===============================================================================

This tutorial describes how to calculate signature dynamics for a family of
proteins with similar structures using Elastic Network Models. This method creates 
an ensemble of aligned structures and calculates statistics such as means and 
standard deviations on various dynamic properties including mode profiles, 
mean square fluctuations and cross-correlation matrices. It also includes tools 
for classifying family members based on their sequence, structure and dynamics.

The theory and usage of this toolkit will be described in [ZB18]_.

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

If you benefited from SignDy Calculations in your research, 
please cite the following paper:

.. [ZB18] Zhang S., Krieger J., Li H., Bahar I. 
   SignDy: Discovering the Signature Dynamics of Protein Families with Elastic Network Model 
   Analysis. *(in preparation)* **2018**
