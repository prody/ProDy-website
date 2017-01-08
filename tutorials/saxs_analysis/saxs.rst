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

Let start with essential import statements:
.. ipython:: python
   import numpy as np
   import matplotlib

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

* Magnitude of scattering vector [q]
* Scattering intensity [I(q)]
* Errors for each I(q) value [Sigma(q)]

If the data is simulated data, it generally contains only q and I(q) values.
Now, lets start parsing SAXS data:

.. ipython:: python
   from saxs import *
   Q_exp, I_q_exp, sigma_q=parseSaxsData(saxs_file, simulated=False, isLogScale=True)

Sometimes, experimental intensities are saved in log scale. Therefore, we
have to specify this information if the SAXS data is in log scale or not. 
		
Compare SAXS profiles of known structure and experimental (simulated) data
-------------------------------------------------------------------------------
We have to produce theoretical SAXS profile of the open conformation. We will
use Fast-SAXS method for this purpose [SY09]_.

.. [SY09] Yang S, Park S, Makowski L, and Roux B. A Rapid Coarse Residue-Based
   Computational Method for X-Ray Solution Scattering Characterization of
   Protein Folds and Multiple Conformational States of Large Protein Complexes.
   *Biophysical Journal*  **2009** 96:4449–4463.

We will call a function to compute SAXS profile of our initial pdb file. 


.. ipython:: python
   I_model=np.zeros(len(Q_exp))
   calcSaxsPerModel(calphas, numCalphas, I_model, Q_exp)

We can write this SAXS profile to a file using Numpy.
It is time to quantify the agreement between the experimental profile and
profile of the initial pdb file. We use Chi value for this purpose. Our purpose
is to produce models by using anisotropic network model that are in better agreement
with the experimental SAXS profile. Of course, better agreement means a lower chi
value. 

.. ipython:: python
   maxChi=calcSaxsChi(Q_exp, I_q_exp, sigma_q, Q_exp, I_model):

  
.. ipython:: python
   showSaxsProfiles(exp_data_file, model_data_file)

.. image:: saxs_4ake_vs_1ake.png
   :scale: 200 %

Build hessian and get normal modes
-------------------------------------------------------------------------------
We will try to reduce Chi value between SAXS profile of initial pdb and
experimental (or simulated) SAXS data by interpolating normal modes obtained
with anisotropic network model. First, lets calculate normal modes.

.. ipython:: python

   anm = ANM('ANM Analysis')
   anm.buildHessian(calphas, cutoff=15.0)
   numCalphas=calphas.numAtoms()
   modes=anm.calcModes(n_modes=5, zeros=False)
		

Interpolate a single normal mode to get a model
-------------------------------------------------------------------------------
If we want to interpolate a single mode and see if it reduces chi values, we
can issue the following. Return value of this function will be two lists,
one for the chi values and the other one for the frame numbers. 

.. ipython:: python

   chi_list, frames_list=interpolateMode(calphas, modes[0], Q_exp, I_q_exp,\
	     sigma_q, maxChi)

Lets see if we have a model in mode 0 that reduces chi value. 
.. ipython:: python

   showChivsFrames(chi_list, frames_list, numFrames=20)

Frame 20 of mode 0 reduces chi value significantly. In this model, lid domain of
adenylate kinase gets closed. Generally, checking just a few normal modes can
be sufficient to observe a large scale conformational change.    

Perform iterations of over multiple low frequnce normal modes
-------------------------------------------------------------------------------
If we want to do iterations over multiple modes and get a model with the lowest
chi value, we can leave interactive ipython prompt and issue the following
command:

    prody saxs 4ake_chainA.pdb 1ake_chainA_saxs_w_yerrorbars.dat -n 5
   
Here, we specified the maximum number of modes to be used as with -n 5 parameter.
   
