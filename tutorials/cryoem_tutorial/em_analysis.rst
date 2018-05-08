.. em_analysis:

Parsing Cryo-EM Electron Density Map
==================================================================

Lets start with essential import statements:

.. ipython:: python

   from prody import *
   import numpy as np
   import matplotlib


Parse Density Map
-----------------------------------------------------------------

The first step is to parse a .map file, which contains information
about a density map as the electron density at points on a grid.
This file format is a binary format also known as CCP4 or MRC2015. 

.. ipython:: python

   emd = parseEMD('emd_1960.map', cutoff=1.2, n_nodes=8000, num_iter=30)
   emd

This function returns an atom group from the electron density
map. The cutoff parameter discards any electron density lower than
the given number, n_nodes is the parameter that describes the
total number of beads that initialize the system. The method to
build coordinates is also called the "Topology Representing
Network" algorithm, which is an iterative algorithm and num_iter 
is therefore a required parameter for this algorithm. 
The resultant structure is given in the following figure. 

Note: Depending on your hardware, this may take a while - around
30 minutes. To skip this step, you can load the structure file
directly in PDB format (emd_1960.pdb).

.. figure:: ../../_static/figures/8000nodes_with_map.png
   :scale: 80%

Elastic Network Model Analysis
-----------------------------------------------------------------

The 144x144x144 density grid is converted into an :class:`.AtomGroup`
class and elastic network model analysis can be applied to the 
constructed structure as usual. 

anm = ANM('TRiC EMDMAP ANM Analysis')
anm.buildHessian(emd)
anm.calcModes()
writeNMD('tric_anm_3_modes_8000nodes.nmd', anm[:3], emd)


