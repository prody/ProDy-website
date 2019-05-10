.. _signdy-core:

Core Calculations
===============================================================================

In order to infer signature dynamics features we create mode ensembles from the 
PDB ensembles by calculating normal modes for each member of the 
:class:`.PDBEnsemble`. For details, please refer to [SZ18]_.

First, make necessary imports from ProDy and Matplotlib packages if you haven't 
already.

.. ipython:: python

   from prody import *
   from pylab import *
   ion()


Mode Ensemble
-------------------------------------------------------------------------------

For this analysis we'll use the provided :class:`.PDBEnsemble` object (link), 
which we built and refined in [SZ18]_, to build a :class:`.ModeEnsemble`. 
There are options to select the model (:class:`.GNM` by default) and the way of 
considering non-aligned residues (default is :func:`.reduceModel`, which treats 
them as environment). 

We first load the ensemble:

.. ipython:: python

   ens = loadEnsemble('LeuT.ens.npz')

Then we calculated first 20 GNM modes for each member of the ensemble:

.. ipython:: python

   gnms = calcEnsembleENMs(ens, model=GNM, trim='reduce', n_modes=20, match=True)
   gnms

In this way, we will obtain one :class:`.ModeSet` for each member and 85 in total. 
Finding a consistent order of modes across different mode sets is critical to the 
accuracy of SignDy_ calculations. Therefore, in the above code, *match* is set 
to *True* so that all other :class:`.ModeSet`s are sorted to match the order of 
the reference :class:`.ModeSet`, which is the first :class:`.ModeSet` in the 
:class:`.ModeEnsemble` by default.


Signature Dynamics
-------------------------------------------------------------------------------

:class:`.Signature`s are calculated as the mean and standard deviation of various 
properties such as mode shapes and mean square fluctations.

For example, we can show the average and standard deviation of the shape of the first 
mode (second index 0). The first index of the mode ensemble is over conformations.

 .. ipython:: python

   @savefig ens_gnms_signature_mode1.png width=4in
   showSignatureMode(gnms[:, 0]);

In the plot, the curve shows the mean values, the darker shade shows the standard 
deviations, and the lighter shade shows the range (minimum and maximum values).
We can also show such things for properties involving multiple modes such as the mean 
square fluctuations from the first 5 modes,

 .. ipython:: python

   @savefig ens_gnms_signature_sqflucts_mode1-5.png width=4in
   showSignatureSqFlucts(gnms[:, :5]);

or the cross-correlations from the first 20.

 .. ipython:: python

   @savefig ens_gnms_signature_cross-corr.png width=4in
   showSignatureCrossCorr(gnms[:, :20]);


We can also look at distributions over values across different members of the ensemble 
such as inverse eigenvalue. We can show a bar above this with individual members labelled 
like [SZ18]_.

 .. ipython:: python

    highlights = {'2A65A': 'LeuT:OF', '3TT1A': 'LeuT:OF', 
                  '3TT3A': 'LeuT:IF', '4US4A': 'MhsT:IF', 
                  '3NCYA': 'AdiC:OF', '2X79A': 'Mhp1:IF', 
                  '2WITA':'BetP', '4M48A':'DAT'}

    figure();
    gs = plt.GridSpec(ncols=1, nrows=2, height_ratios=[1, 10], hspace=0.15)

    subplot(gs[0]);
    showVarianceBar(gnms[:, :5], fraction=True, highlights=highlights);
    xlabel('');

    subplot(gs[1]);
    showSignatureVariances(gnms[:, :5], fraction=True, bins=80, alpha=0.7);
    xlabel('Mode weight');

    @savefig ens_gnms_signature_variance_mode1-5.png width=4in

Saving the ModeEnsemble
-------------------------------------------------------------------------------

Finally we save the mode ensemble for later processing:

.. ipython:: python

   saveModeEnsemble(gnms, 'LeuT')

