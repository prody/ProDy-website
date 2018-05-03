Core Calculations
===============================================================================

In order to infer signature dynamics features we create mode ensembles from the 
PDB ensembles by calculating normal modes for each member of the PDB ensemble. 

First, make necessary imports from ProDy and Matplotlib packages if you haven't already.

.. ipython:: python

   from prody import *
   from pylab import *
   ion()

PDB Ensemble
-------------------------------------------------------------------------------

The first step in signature dynamics analysis is to collect a set of related 
protein structures and build a :class:`PDBEnsemble`. This can be achieved by 
multiple routes. In this example, we start with a user-defined list of PDB IDs 
corresponding to transporters with the LeuT fold. These are aligned using the 
CE algorithm

.. ipython:: python

   PDB_IDs = ['2A65A', '2Q6HA', '2Q72A', '2QB4A',
 '2QEIA',
 '2QJUA',
 '3F3AA',
 '3F3CA',
 '3F3DA',
 '3F3EA',
 '3F48A',
 '3F4IA',
 '3F4JA',
 '3GJCA',
 '3GJDA',
 '3GWUA',
 '3GWVA',
 '3GWWA',
 '3MPNA',
 '3MPQA',
 '3QS4A',
 '3QS5A',
 '3QS6A',
 '3TT1A',
 '3TU0A',
 '3USGA',
 '3USIA',
 '3USJA',
 '3USKA',
 '3USLA',
 '3USMA',
 '3USOA',
 '3USPA',
 '4FXZA',
 '4FY0A',
 '4HMKA',
 '4HODA',
 '4MM4A',
 '4MM5A',
 '4MM6A',
 '4MM7A',
 '4MM8A',
 '4MM9A',
 '4MMAA',
 '4MMBA',
 '4MMCA',
 '4MMDA',
 '4MMEA',
 '4MMFA',
 '3TT3A',
 '4US4A',
 '4US3A',
 '4M48A',
 '4XNUA',
 '4XNXA',
 '4XP1A',
 '4XP4A',
 '4XP5A',
 '4XP6A',
 '4XP9C',
 '4XPAA',
 '4XPBA',
 '4XPFA',
 '4XPGA',
 '4XPHA',
 '4XPTA',
 '2JLNA',
 '2X79A',
 '4D1AA',
 '4D1BA',
 '4D1CA',
 '4D1DA',
 '2WITA',
 '2WITB',
 '2WITC',
 '3P03A',
 '3P03B',
 '3P03C',
 '4AINA',
 '4AINB',
 '4AINC',
 '4C7RA',
 '4C7RB',
 '4C7RC',
 '4DOJA',
 '4DOJB',
 '4DOJC',
 '4LLHA',
 '4LLHB',
 '4LLHC',
 '2XQ2A',
 '3L1LA',
 '3LRBA',
 '3LRCA',
 '3NCYA',
 '3OB6A',
 '5J4IA',
 '5J4NA',
 '4M8JA',
 '2WSXA',
 '2WSXB',
 '2WSXC',
 '2WSWA',
 '3HFXA']


Mode Ensemble
-------------------------------------------------------------------------------

For this analysis we'll use the PDB ensemble we built with Dali. There are also options 
to select the model (by default it uses the GNM) and way of considering non-aligned residues 
(default is to use :func:`reduceModel`, which treats them as environment). 

.. ipython:: python

   mode_ens = calcEnsembleENMs(dali_ens)
   mode_ens


Signature dynamics
-------------------------------------------------------------------------------

Signatures are calculated as the mean and standard deviation of various properties 
such as mode shapes and mean square fluctations.

For example, we can show the average and standard deviation of the shape of the first 
mode (second index 0). The first index of the mode ensemble is over conformations.

 .. ipython:: python

   show = showSignatureMode(ens_gnms[:,0])


We can also show such things for properties involving multiple modes such as the mean 
square fluctuations from the first 5 modes or the cross-correlations from the first 20.

 .. ipython:: python

   show = showSignatureSqFlucts(ens_gnms[:,:5]
   plt.figure()
   show = showSignatureCrossCorr(ens_gnms[:,:20]


We can also look at distributions over values across different members of the ensemble 
such as inverse eigenvalue. We can show a bar above this with individual members labelled 
like [KB15]_.

 .. ipython:: python

    highlights= ['3H5V_AB', '3O21_CD', '4MS3_AB', '1ISR_AA', '5CNM_AA', '1EWT_AB', '3QEM_CD']
    protnames = {'3H5V_AB': 'GluA2','3O21_CD': 'GluA3', 
                '4MS3_AB': 'GABAB:A', '1ISR_AA': 'mGlu1:A',
                '5CNM_AA': 'mGlu3:R', '1EWT_AB': 'mGlu1:R',
                '3QEM_CD': 'NMDAR'}

    plt.figure()
    shape = (10, 1)
    gs = plt.GridSpec(ncols=1, nrows=2, height_ratios=[1, 10], hspace=0.15)

    ### cumulative modes variance bar ###
    plt.subplot(gs[0])
    bar, annotations = showVarianceBar(ens_gnms[:, 0], fraction=True, highlights=highlights)
    for ann in annotations:
        text = ann.get_text()
        if text in protnames:
            ann.set_text(protnames[text])
    plt.xlabel('')

### mode variance distributions ###
plt.subplot(gs[1])
showSignatureVariances(ens_gnms[:, :5], fraction=True, bins=80, alpha=0.7)
plt.xlabel('Fraction of inverse eigenvalue')


.. [KB15] Krieger J, Bahar I, Greger IH
    Structure, Dynamics, and Allosteric Potential of Ionotropic Glutamate Receptor N-Terminal Domains.
    *Biophys. J.* **2015** 109(6):1136-48
