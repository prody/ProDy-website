Core Calculations
===============================================================================

In order to infer signature dynamics features we create mode ensembles from the 
PDB ensembles by calculating normal modes for each member of the :class:`.PDBEnsemble`. 

First, make necessary imports from ProDy and Matplotlib packages if you haven't already.

.. ipython:: python

   from prody import *
   from pylab import *
   ion()


Mode Ensemble
-------------------------------------------------------------------------------

For this analysis we'll use the :class:`.PDBEnsemble` object we built with Dali. 
There are also options to select the model (by default it uses the GNM) and 
the way of considering non-aligned residues (default is to use :func:`.reduceModel`, 
which treats them as environment). 

We first load the ensemble:

.. ipython:: python

   dali_ens = loadEnsemble('dali_ensemble.ens.npz')

Then we calculated GNM modes for each member of the ensemble:

.. ipython:: python

   ens_gnms = calcEnsembleENMs(dali_ens)
   ens_gnms


Signature dynamics
-------------------------------------------------------------------------------

Signatures are calculated as the mean and standard deviation of various properties 
such as mode shapes and mean square fluctations.

For example, we can show the average and standard deviation of the shape of the first 
mode (second index 0). The first index of the mode ensemble is over conformations.

 .. ipython:: python

   @savefig ens_gnms_signature_mode1.png width=4in
   show = showSignatureMode(ens_gnms[:,0])


We can also show such things for properties involving multiple modes such as the mean 
square fluctuations from the first 5 modes or the cross-correlations from the first 20.

 .. ipython:: python

   @savefig ens_gnms_signature_sqflucts_mode1-5.png width=4in
   show = showSignatureSqFlucts(ens_gnms[:,:5])


 .. ipython:: python

   @savefig ens_gnms_signature_cross-corr.png width=4in
   show = showSignatureCrossCorr(ens_gnms[:,:20])


We can also look at distributions over values across different members of the ensemble 
such as inverse eigenvalue. We can show a bar above this with individual members labelled 
like [KB15]_.

 .. ipython:: python

    highlights= ['3h5vA_ca', '3o21C_ca',
                 '3h6gA_ca', '3olzA_ca', 
                 '5kc8A_ca']
    protnames = {'3h5vA_ca': 'GluA2','3o21C_ca': 'GluA3',
                 '3h6gA_ca': 'GluK2', '3olzA_ca': 'GluK3', 
                 '5kc8A_ca': 'GluD2'}

    plt.figure();
    shape = (10, 1)
    gs = plt.GridSpec(ncols=1, nrows=2, height_ratios=[1, 10], hspace=0.15)

    plt.subplot(gs[0]);
    bar, annotations = showVarianceBar(ens_gnms[:, :5], fraction=True, highlights=highlights)
    for ann in annotations:
        text = ann.get_text()
        if text in protnames:
            ann.set_text(protnames[text])
    plt.xlabel('');

    plt.subplot(gs[1]);
    show = showSignatureVariances(ens_gnms[:, :5], fraction=True, bins=80, alpha=0.7)
    plt.xlabel('Fraction of inverse eigenvalue');

    @savefig ens_gnms_signature_variance_mode1-5.png width=4in
    plt.show()


Saving the ModeEnsemble
-------------------------------------------------------------------------------

Finally we save the mode ensemble for later processing:

.. ipython:: python

   saveModeEnsemble(ens_gnms, 'gnm_ensemble')

.. [KB15] Krieger J, Bahar I, Greger IH
    Structure, Dynamics, and Allosteric Potential of Ionotropic Glutamate Receptor N-Terminal Domains.
    *Biophys. J.* **2015** 109(6):1136-48
