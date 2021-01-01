Introduction
===============================================================================

This tutorial shows how to use cryo-em electron density data to perform 
elastic network model analysis as performed in [YZ20]_. The electron density maps could be gathered 
from the electron density map database EMDDatabank_. 

We demonstrate this method for an AMPA-type iGluR rather than the chaperonin CCT 
to speed up calculations as it is smaller.


Required Programs
-------------------------------------------------------------------------------

Latest version of ProDy_ is required.

Recommended Programs
-------------------------------------------------------------------------------

We recommend a visualization program too, such as VMD_.

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
   ion()

We have included these imports in every part of the tutorial, so that
code copied from the online pages is complete. You do not need to repeat
imports in the same Python session.

.. _EMDDatabank: http://www.emdatabank.org

.. [YZ20] Yan Zhang, James Krieger, Karolina Mikulska-Ruminska, Burak Kaynak, 
    Carlos Oscar S. Sorzano, José-María Carazo, Jianhua Xing, Ivet Bahar 
    State-dependent sequential allostery exhibited by chaperonin TRiC/CCT 
    revealed by network analysis of Cryo-EM maps,
    *Prog. Biophys. Mol. Biol.* **In Press** 10.1016/j.pbiomolbio.2020.08.006.
