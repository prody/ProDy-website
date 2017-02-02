.. em_analysis:

Parsing Cryo-EM Electron Density Map
==================================================================

Lets start with essential import statements:

.. ipython:: python

   from prody import *
   import numpy as np
   import matplotlib

Tutorial files
-----------------------------------------------------------------

We need the following file for this tutorial: 

.. literalinclude:: files.txt

* :file:`emd_1960.map`
* :file:`emd_1960.pdb`


The first file is mandatory for inputs since it contains the electron
densities in each grid. The second file is a short cut in tutorial to
skip the coordinate building step. 

Parse Density Map
-----------------------------------------------------------------

Density map contains the electron density for 144x144x144 grids
given in the .map file. This file format is in a binary format 
also known as CCP4. 

.. ipython:: python

   emd = parseEMD('emd_1960.map', cutoff=1.2, n_nodes = 8000, num_iter=30)

This function returns an atom group from the electron density
map. Cutoff parameter discards the electron density lower than
the given number, n_nodes is the parameter that describes the
total number of beads that initialize the system. The method to
build coordinates is also called as "Topology Representing
Network" which is an iterative algorithm and num_iter is the 
parameter for this algorithm. The resultant structure is given
in the following figure. 

Note: Depending on your hardware, this may take a while. Around
30 minutes. To skip this step, you can load the structure file
directly in PDB format (emd_1960.pdb).    

.. figure:: _static/figures/mol_rep.png
   :scale: 80%

Elastic Network Model Analysis
-----------------------------------------------------------------

The 144x144x144 density grid converted into an :class:`.AtomGroup`
class and elastic network model analysis are applied to the 
constructed structure. 

anm = ANM('TRiC EMDMAP ANM Analysis')
anm.buildHessian(emd)
anm.calcModes()
writeNMD('tric_anm_3_modes_8000nodes.nmd', anm[:3], emd)


