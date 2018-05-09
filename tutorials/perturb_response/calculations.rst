Calculations
===============================================================================

Here are the required imports again. You do not need to repeat them if you are
still in the same python session.

.. ipython:: python

   from prody import *
   from pylab import *
   ion()


Loading a structure and applying the anisotropic network model
-------------------------------------------------------------------------------
First, we parse a structure that we want to analyse with PRS. 
For this tutorial, we will fetch the near-intact full length AMPAR structure 3kg2 from the PDB. 
We import a subset containing the calpha atoms, which we will use in downstream steps.

.. ipython:: python

    ampar_ca = parsePDB('3kg2', subset='ca')


Next, create an ANM instance and calculate modes from which the covariance matrix can be calculated. 
We ask for all modes rather than the default subset of the first 20 global modes. We could alternatively 
apply the PRS to another model from which a covariance matrix could be derived such as PCA, GNM or an MD simulation.

.. ipython:: python

    anm_ampar = ANM('AMPAR 3kg2')
    anm_ampar.buildHessian(ampar_ca)
    anm_ampar.calcModes(n_modes='all')


You can also save the model in a .npz file for loading again and write an .nmd for visualizing in the NMWiz_.

.. ipython:: python

    saveModel(anm_ampar,'3kg2',matrices=True)
    writeNMD('anm_3kg2.nmd',anm_ampar,ampar_ca)


Showing the PRS matrix with effectiveness and sensitivity profiles
-------------------------------------------------------------------------------

The PRS matrix is then calculated from the covariance matrix from the ANM of the AMPAR. 
This can be automatically carried out together with plotting using the following code.
We supply atoms to allow us to delineate residue numbers and chains.

.. ipython:: python

    show = showPerturbResponse(model=anm_ampar, atoms=ampar_ca)
    @savefig 3kg2_prs.png width=4in
    plt.show()

The next page illustrates some more detailed analyses.

.. _NMWiz: http://prody.csb.pitt.edu/nmwiz/
