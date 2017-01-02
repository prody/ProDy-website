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


Compare SAXS profiles of known structure and experimental (simulated) data
-------------------------------------------------------------------------------


Build hessian and get normal modes
-------------------------------------------------------------------------------


Iterate a normal mode
-------------------------------------------------------------------------------



