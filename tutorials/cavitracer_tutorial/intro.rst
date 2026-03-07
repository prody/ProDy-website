Introduction
===============================================================================

This tutorial shows how to detect intraprotein channels, tunnels, and
cavities. This method can be applied to single PDB file, PDB enseble, and
molecular dynamics (MD) trajectory.


Required Programs
-------------------------------------------------------------------------------

The latest version of ProDy_ is required.


Recommended Programs
-------------------------------------------------------------------------------

Besides ProDy_, the Matplotlib_ library and VMD_ program are required for
some steps in the tutorial. IPython_ is highly recommended for interactive usage.

Moreover, ProDy provides visualization of the predictions using Open3D_. 

.. _Open3D: https://pypi.org/project/open3d/

It can be installed using a conda or pip.

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


First, we will make necessary imports from ProDy and Matplotlib packages.

.. ipython:: python

   from prody import *
   from pylab import *
   import matplotlib

We have included these imports in every part of the tutorial so that the
code copied from the online pages is complete. You do not need to repeat
imports in the same Python session.


How to Cite
-------------------------------------------------------------------------------

If you benefited from CaviTracer in your research, 
please cite [ET25]_. Since the implementation is inspired by the methods described in the
publication [DS13]_, we recommend to cite it too.


.. [ET25] E. Trzcinski, J. M. Krieger, J. Duda, K. Mikulska-Ruminska,
BPS2025-CaviTracer: A ProDy tool for mapping protein tunnels and channels
with application in lipoxygenases. *Biophys J* 124 (3) **2025** 552a (abstract).

.. [DS13] D. Sehnal, et al.,
MOLE 2.0 advanced approach for analysis of biomacromolecular channels.
*J Chemoinform* 5 (39) **2013**.
