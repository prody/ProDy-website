Introduction
===============================================================================

This tutorial shows how to use cryo-EM electron density data to perform 
elastic network model analysis as in [YZ20]_. Electron density maps can 
be provided by the user or retrieved from the EMDB_. 

Cryo-electron microscopy (cryo-EM) is a type of transmission 
electron microscopy where the studied sample, e.g. protein 
solution, is rapidly frozen. In structural biology, single 
particle cryo-EM is gaining more and more popularity due to 
its advantage of allowing the observation of protein under 
near-native conditions at ever-improving resolution, and as 
a result, the structures of numerous mega-Dalton biomolecules 
are captured and stored in the form of cryo-EM density maps, 
with resolutions ranging from 2 to 100 Å. 

Fitting all-atom structures to a density map is time-consuming, 
and sometimes not feasible due to low resolution. In addition, 
for the purpose of studying the dynamics, MD simulations of 
supercomplexes with atomic details are extremely 
computationally expensive. It is thus desirable to apply 
coarse-grained methods, for example the Anisotropic Network Model 
(ANM), to such data to study the “big motions” of molecular 
machines beyond the atomic level. 

ProDy can be used for constructing a bead-and-spring model from cryo-EM data 
using a previously published algorithm ([TM94]_), which can be used for ANM.


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

.. _EMDB: https://www.ebi.ac.uk/pdbe/emdb/

.. [YZ20] Yan Zhang, James Krieger, Karolina Mikulska-Ruminska, Burak Kaynak, 
    Carlos Oscar S. Sorzano, José-María Carazo, Jianhua Xing, Ivet Bahar 
    State-dependent sequential allostery exhibited by chaperonin TRiC/CCT 
    revealed by network analysis of Cryo-EM maps,
    *Prog. Biophys. Mol. Biol.* **In Press** 10.1016/j.pbiomolbio.2020.08.006.

.. _[TM94]: Martinetz and Schulten, 1994
   Topology representing networks
   *Neural Networks* **1994** 7:507-552.
