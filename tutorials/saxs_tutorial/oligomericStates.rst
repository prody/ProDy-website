.. _oligomericStates:

Determining oligomeric states of a protein
===============================================================================
Proteins can exist in various oligomeric forms in solution. They can be found in
solutions as monomers, dimers, tetramers or even a combination of these
multimeric forms. SAXS can be used to determine oligomeric states of a protein
in solution. We will show how to perform this analysis by using Aldehyde
dehydrogenase 7A1 (ALDH7A1) as an example.

Lets start with essential import statements:

.. ipython:: python

   import numpy as np
   import matplotlib


Tutorial files
-------------------------------------------------------------------------------
We need three files for this tutorial:

.. literalinclude:: files.txt

* :file:`4zul.pdb`
* :file:`ALDH7A1SAXS.dat`

We should not that experimental data file :file:`ALDH7A1SAXS.dat` was obtained
from [JT15]_
   
.. [JT15] Tanner J. SAXS Fingerprints of Aldehyde Dehydrogenase Oligomers.
	  *Data in Brief* **2015** 5:745-751.

We wonder if solution structure of ALDH7A1 is dimeric or tetrameric. The essential
assumption here is that the protein given in PDB file do not undergo a significant
conformational transition in solution. 
   

Parse known structure
-------------------------------------------------------------------------------
We start with parsing a PDB file by passing an identifier.
Note that if a file is not found in the current working directory, it will be
downloaded.

.. ipython:: python

   protein = parsePDB('4zul.pdb')

We want to use only Cα atoms in chain A, so we select them:

.. ipython:: python

   dimer_calphas = protein.select('protein and name CA and (chain A or chain B)')

.. ipython:: python

   tetramer_calphas = protein.select('protein and name CA and (chain A or chain B or chain C or chain D)')


Parse SAXS data
-------------------------------------------------------------------------------
We need to parse experimental SAXS data and save it into Numpy
arrays for further use. An experimental SAXS profile contains three columns:

* q: Magnitude of scattering vector or scattering angle. We will use 1/Angstrom a
  as q unit. 
* I(q): Scattering intensity, which may be in logarithmic scale. 
* Sigma(q): Errors for each I(q) value 


Please note that the SAXS data file should only contain data. It should not
contain any comment or other text lines. Now, lets start parsing SAXS data:

.. ipython:: python

   from prody.dynamics.saxs import *
   Q_exp, I_q_exp, sigma_q=parseSaxsData('ALDH7A1SAXS.dat', isLogScale=False)


Sometimes, experimental intensities are saved in log scale as in our case.
Therefore, we have to specify this information if the SAXS data is in log scale
or not. 
		
Compare SAXS profiles of known structure and experimental data
-------------------------------------------------------------------------------
We have to produce theoretical SAXS profile of the open conformation. We will
use Fast-SAXS method for this purpose [SY09]_.

.. [SY09] Yang S, Park S, Makowski L, and Roux B. A Rapid Coarse Residue-Based
   Computational Method for X-Ray Solution Scattering Characterization of
   Protein Folds and Multiple Conformational States of Large Protein Complexes.
   *Biophysical Journal*  **2009** 96:4449–4463.

We will call :function:`calcSaxsPerModel` to compute SAXS profile of our initial
pdb file for each experimental q value. Therefore, we have to use Q_exp array
we have just read from the experimental data file. Theoretical
SAXS intensities will be saved to I_q_model array. We should note that the
theoretical SAXS intensities are also produced in log scale. Lets calculate SAXS
profile for dimeric form of ALDH7A1 and compare it to the experimental SAXS data
first.

.. ipython:: python

   I_q_model_dimer=np.zeros(len(Q_exp))
   calcSaxsPerModel(dimer_calphas, I_q_model_dimer, Q_exp)

We can write this SAXS profile to a file using :function:`writeSaxsProfile`.
Lets write this model to a file for further investigation. 
   
.. ipython:: python

   writeSaxsProfile(I_q_model_dimer, Q_exp, sigma_q, 'dimerSAXS.dat')

Now, lets calculate the SAXS profile for the tetrameric form of ALDH7A1:

.. ipython:: python

   I_q_model_tetramer=np.zeros(len(Q_exp))
   calcSaxsPerModel(tetramer_calphas, I_q_model_tetramer, Q_exp)


Lets write this model as well to a file for further investigation. 
   
.. ipython:: python

   writeSaxsProfile(I_q_model_tetramer, Q_exp, sigma_q, 'tetramerSAXS.dat')

Finally, we have to write experimental SAXS data to a file in log scale to
make a comparison.

.. ipython:: python

   writeSaxsProfile(I_q_exp, Q_exp, sigma_q, 'ALDH7A1SAXS_logscale.dat')

Initially, we can plot SAXS profile of the dimer vs the experimental SAXS
profile. 

.. ipython:: python

   showSaxsProfiles('ALDH7A1SAXS_logscale.dat', 'dimerSAXS.dat')
   
.. image:: images/dimerVsExperiment.png 
   :width: 4in
   
Obviously, the theoretical SAXS profile do not agree well with the
experiment. Lets plot the theoretical profile for the tetrameric form vs
experimental SAXS profile.  
   
.. ipython:: python

   showSaxsProfiles('ALDH7A1SAXS_logscale.dat', 'tetramerSAXS.dat')
   @savefig images/tetramerVsExperiment.png width=4in


We can see that the SAXS data suggests that the protein is in tetrameric
form in solution.    
