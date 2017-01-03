.. _saxs:

Small Angle X-ray Scattering (SAXS) Analysis
===============================================================================
SAXS is an experimental technique which gives information about solution
structure of biological macromolecules. SAXS data is a 1D data of Magnitude of
Scattering Vector vs Scattering Intesity. In this tutorial, we will see how to
combine SAXS data with a known protein structure to obtain information about
solution structure of a protein. Basically, we are trying to obtain an unknown
protein conformation from a known pdb file and SAXS data by utilizing normal
modes obtained with anisotropic network model. We will use a widely known test
case, adenylate kinase. We assume that we have pdb file for open conformation
(4ake) and SAXS profile of closed conformation (1ake). We are trying to recover
the closed conformation by combining initial pdb file with SAXS data. 

Parse known structure
-------------------------------------------------------------------------------
We start with parsing a PDB file by passing an identifier.
Note that if a file is not found in the current working directory, it will be
downloaded.

.. ipython:: python

	        protein = parsePDB('4ake')

We want to use only Cα atoms in chain A, so we select them:

.. ipython:: python

	        calphas = protein.select('protein and name CA and chain A')
		calphas
		   

Parse SAXS data
-------------------------------------------------------------------------------
We need to parse experimental  SAXS data and save it into Numpy
arrays for further use. An experimental SAXS profile contains three columns:
    Magnitude of scattering vector [q]
    Scattering intensity [I(q)]
    Errors for each I(q) value [Sigma(q)]
If the data is simulated data, it generally contains only q and I(q) values.
Now, lets start parsing SAXS data:
.. ipython:: python

		from saxs import *
		Q_exp, I_q_exp, sigma_q=parseSaxsData(saxs_file, \
		simulated=False, isLogScale=True)

Sometimes, experimental intensities are saved in log scale. Therefore, we
have to specify this information if the SAXS data is in log scale or not. 
		
Compare SAXS profiles of known structure and experimental (simulated) data
-------------------------------------------------------------------------------
We have to produce theoretical SAXS profile of the open conformation. We will
use Fast-SAXS method as explained in [SY09]_.

.. [SY09] Yang S, Park S, Makowski L, and Roux B. A Rapid Coarse Residue-Based
   Computational Method for X-Ray Solution Scattering Characterization of
   Protein Folds and Multiple Conformational States of Large Protein Complexes.
   *Biophysical Journal*  **2009** 96:4449–4463.


.. image:: gnu.png

Build hessian and get normal modes
-------------------------------------------------------------------------------
We will try to reduce Chi value between SAXS profile of initial pdb and
experimental (or simulated) SAXS data by interpolating normal modes obtained
with anisotropic network model. First, lets calculate normal modes.
.. ipython:: python

		anm = ANM('ANM Analysis')
		anm.buildHessian(calphas, cutoff=15.0)
		numCalphas=calphas.numAtoms()
		modes=anm.calcModes(n_modes=3, zeros=False)
		

Iterate a single normal mode
-------------------------------------------------------------------------------



Perform iterations of over multiple low frequnce normal modes
-------------------------------------------------------------------------------

