.. _signdy-overview:

Overview
===============================================================================

The first step in signature dynamics analysis is to collect a set of related 
protein structures and build a :class:`.PDBEnsemble`. This can be achieved by 
multiple routes: a query search of the PDB using :func:`.blastPDB` or 
:func:`.searchDali`, extraction of PDB IDs from the Pfam or CATH database, or 
input of a pre-defined list. 

We demonstrate the `Dali`_ method here in the first part of the tutorial. The usage of 
`Pfam`_ and `CATH`_ methods are described in the database tutorial (under construction) 
and the function :func:`.blastPDB` is described in the `Structure Analysis Tutorial`_.

We apply these methods to the PBP-I domains, a group of protein structures originally 
found in bacteria for transport of solutes across the periplasmic space and later 
seen in various eukaryotic receptors including ionotropic and metabotropic glutamate 
receptors. We use the N-terminal domain of AMPA receptor subunit GluA2 (gene name 
`GRIA2 <https://www.uniprot.org/uniprot/P42262>`_) as a query.

First, make necessary imports from ProDy and Matplotlib packages if you haven't already.

.. ipython:: python

    from prody import *
    from pylab import *
    ion()

Prepare ensemble (via Dali)
-------------------------------------------------------------------------------

First we use the function :func:`.searchDali` to search the PDB, which returns a 
:class:`.DaliRecord` that contains a list of PDB IDs and their corresponding 
mappings to the reference structure. 

.. ipython:: python

    dali_rec = searchDali('3H5V','A')
    dali_rec

Next, we get the lists of PDB IDs and mappings from *dali_rec*, parse the *pdb_ids* 
to get a list of :class:`.AtomGroup` instances:

.. ipython:: python

    pdb_ids = dali_rec.getPDBs()
    mappings = dali_rec.getMappings()

.. ipython:: python

    pdbs = parsePDB(*pdb_ids, subset='ca')
    len(pdbs)

Then we provide *pdbs* together with *mappings* to :func:`.buildPDBEnsemble`. We 
set the keyword argument ``seqid=20`` to account for the low sequence identity 
between some of the structures.

.. ipython:: python

    dali_ens = buildPDBEnsemble(pdbs, mapping=mappings, seqid=20)
    dali_ens

Finally we save the ensemble for later processing:

.. ipython:: python

   saveEnsemble(dali_ens, 'PBP-I')

Mode ensemble
-------------------------------------------------------------------------------

For this analysis we'll build a :class:`.ModeEnsemble` by calculating normal 
modes for each member of the :class:`.PDBEnsemble`. First, we load the ensemble:

.. ipython:: python

   dali_ens = loadEnsemble('PBP-I.ens.npz')

Then we calculated :class:`.GNM` modes for each member of the ensemble. There 
are options to select the *model* (:class:`.GNM` by default) and the way of 
considering non-aligned residues by setting the *trim* option (default is 
:func:`.reduceModel`, which treats them as environment). :

.. ipython:: python

   gnms = calcEnsembleENMs(dali_ens, model='GNM', trim='reduce')
   gnms


Signature dynamics
-------------------------------------------------------------------------------

Signatures are calculated as the mean and standard deviation of various properties 
such as mode shapes and mean square fluctations.

For example, we can show the average and standard deviation of the shape of the first 
mode (second index 0). The first index of the mode ensemble is over conformations.

 .. ipython:: python

   @savefig signdy_dali_mode1.png width=4in
   showSignatureMode(gnms[:, 0]);


We can also show such things for properties involving multiple modes such as the mean 
square fluctuations from the first 5 modes or the cross-correlations from the first 20.

 .. ipython:: python

   @savefig signdy_dali_mode1-5.png width=4in
   showSignatureSqFlucts(gnms[:, :5]);


 .. ipython:: python

   @savefig signdy_dali_cross-corr.png width=4in
   showSignatureCrossCorr(gnms[:,:20]);


We can also look at distributions over values across different members of the ensemble 
such as inverse eigenvalue. We can show a bar above this with individual members labelled 
like [KB15]_.

 .. ipython:: python

    highlights = {'3h5vA_ca': 'GluA2','3o21C_ca': 'GluA3',
                 '3h6gA_ca': 'GluK2', '3olzA_ca': 'GluK3', 
                 '5kc8A_ca': 'GluD2'}

    @savefig signdy_dali_variance_mode1-5.png width=4in
    figure();
    gs = GridSpec(ncols=1, nrows=2, height_ratios=[1, 10], hspace=0.15)

    subplot(gs[0]);
    showVarianceBar(gnms[:, :5], fraction=True, highlights=highlights);
    xlabel('');

    subplot(gs[1]);
    showSignatureVariances(gnms[:, :5], fraction=True, bins=80, alpha=0.7);
    xlabel('Fraction of inverse eigenvalue');

Finally we save the mode ensemble for later processing:

.. ipython:: python

   saveModeEnsemble(gnms, 'PBP-I')


Spectral overlap and distance
-------------------------------------------------------------------------------

Spectral overlap is also known as covariance overlap as defined in [BH02]_. 
Covariance overlap measures the distance between two covariance matrices, but in 
our case, we will generalize it to calculate the overlap of a subset of the modes 
(or spectrum). We first load the :class:`.ModeEnsemble`:

.. ipython:: python

   gnms = loadModeEnsemble('PBP-I.modeens.npz')

We calculate the spectral overlap matrix, calculate a tree from its arccosine 
(to convert the overlap to distance):

.. ipython:: python

    so_matrix = calcEnsembleSpectralOverlaps(gnms[:, :1])
    labels = gnms.getLabels()
    so_tree = calcTree(names=labels, 
                       distance_matrix=arccos(so_matrix), 
                       method='upgma')

We can reorder the spectral overlap matrix using the tree as follows: 

.. ipython:: python

    reordered_so, new_so_indices = reorderMatrix(so_matrix, 
                                                 so_tree, 
                                                 names=labels)

Both :class:`.PDBEnsemble` and :class:`.ModeEnsemble` objects can be reordered 
based on the new indices:

.. ipython:: python

    reordered_ens = dali_ens[new_so_indices]
    reordered_gnms = gnms[new_so_indices, :]


Compare with sequence and structural distances
-------------------------------------------------------------------------------

The sequence distance is given by the (normalized) Hamming distance, which is 
calculated by subtracting the percentage identity (fraction) from 1, and the 
structural distance is the RMSD. We can also calculate and show the matrices 
and trees for these from the PDB ensemble.

First we calculate the sequence distance matrix:

.. ipython:: python

    seqid_matrix = buildSeqidMatrix(ens.getMSA())
    seqd_matrix = 1. - seqid_matrix

We can visualize the matrix using :func:`.showMatrix`:

.. ipython:: python

    @savefig signdy_dali_seqd_matrix.png width=4in
    showMatrix(seqd_matrix);

We can also construct a tree based on the distance matrix:

.. ipython:: python

    seqd_tree = calcTree(names=labels, 
                         distance_matrix=seqd_matrix, 
                         method='upgma')

Similarily, once we obtain the RMSD matrix using :meth:`.PDBEnsemble.getRMSD`, we 
can calculate the structure-based tree:

.. ipython:: python

    rmsd_matrix = ens.getRMSDs(pairwise=True)
    @savefig signdy_dali_rmsd_matrix.png width=4in
    showMatrix(rmsd_matrix);

    rmsd_tree = calcTree(names=labels, 
                         distance_matrix=rmsd_matrix, 
                         method='upgma')

It could be of interest to put all three trees constructed based on different 
distance metrics side by side and compare them:

.. ipython:: python

    @savefig signdy_trees.png width=4in
    figure();
    subplot(1, 3, 1);
    showTree(seqd_tree, format='plt');
    title('Sequence');
    subplot(1, 3, 2);
    showTree(rmsd_tree, format='plt');
    title('Structure');
    subplot(1, 3, 3);
    showTree(so_tree, format='plt');
    title('Dynamics');

This analysis is quite sensitive to how many modes are used. As the number of modes approaches the full number, 
the dynamic distance order approaches the RMSD order. With smaller numbers, we see finer distinctions. This is 
particularly clear in the current case where we used just one mode.

.. [KB15] Krieger J, Bahar I, Greger IH.
    Structure, Dynamics, and Allosteric Potential of Ionotropic Glutamate Receptor N-Terminal Domains.
    *Biophys. J.* **2015** 109(6):1136-48

.. _`Structure Analysis Tutorial`: http://prody.csb.pitt.edu/tutorials/structure_analysis/blastpdb.html
.. _`list_comprehensions`: https://docs.python.org/2/tutorial/datastructures.html#list-comprehensions
.. _`Dali`: http://ekhidna2.biocenter.helsinki.fi/dali/
.. _`Pfam`: https://pfam.xfam.org/
.. _`CATH`: http://www.cathdb.info/
