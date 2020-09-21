.. _signdy-class:

Classification using sequence, structure and dynamics distances
===============================================================================

We can compare the dynamics of individual proteins using the spectral overlap, 
also known as covariance overlap. The arccosine of this value provides a distance 
metric. Calculating this for all pairs in a mode ensemble gives us the spectral distance 
matrix, which can be used to calculate a dynamics-based "phylogenetic" tree. This can be 
compared against matrices and trees calculated using sequence and structure distances.

Again here are the imports if you need them.

.. ipython:: python

    from prody import *
    from pylab import *
    ion()

Load PDBEnsemble and ModeEnsemble
-------------------------------------------------------------------------------

We first load the :class:`.PDBEnsemble`:

.. ipython:: python

    ens = loadEnsemble('LeuT.ens.npz')

Then we load the :class:`.ModeEnsemble`:

.. ipython:: python

    gnms = loadModeEnsemble('LeuT.modeens.npz')

Spectral overlap and distance
-------------------------------------------------------------------------------

We calculate the spectral overlap matrix, calculate a tree from its arccosine and 
reorder the spectral overlap matrix using the tree as follows: 

.. ipython:: python

    so_matrix = calcEnsembleSpectralOverlaps(gnms[:,:1])
    labels = gnms.getLabels()
    so_tree = calcTree(names=labels, 
                       distance_matrix=arccos(so_matrix), 
                       method='upgma')

    reordered_so, new_so_indices = reorderMatrix(names=labels,
                                                 matrix=so_matrix, 
                                                 tree=so_tree)


We can show the original and reordered spectral distance matrices and the tree as follows.
:func:`.showTree` has multiple *format* options. Here we show the output of using *plt*.
This layout allows us to directly compare against the output from :func:`.showMatrix`
using the option *origin='upper'*.

.. ipython:: python

    @savefig ens_gnms_so_matrix.png width=4in
    showMatrix(arccos(so_matrix), origin='upper')
	
    @savefig ens_gnms_so_tree.png width=4in
    showTree(so_tree, format='plt')
	
    @savefig ens_gnms_so_reordered_so_matrix.png width=4in
    showMatrix(arccos(reordered_so), origin='upper')

    plt.close('all')

Sequence and structural distances
-------------------------------------------------------------------------------

The sequence distance is given by the Hamming distance, which is calculated by 
subtracting the percentage identity (fraction) from 1, and the structural distance 
is the RMSD. We can also calculate and show the matrices and trees for these from 
the PDB ensemble.

.. ipython:: python

    seqid_matrix = buildSeqidMatrix(ens.getMSA())
    seqd_matrix = 1. - seqid_matrix
    @savefig ens_gnms_seqd_matrix.png width=4in
    showMatrix(seqd_matrix, origin='upper')

    plt.figure()
    seqd_tree = calcTree(names=labels, 
                         distance_matrix=seqd_matrix, 
                         method='upgma')
    @savefig ens_gnms_seqd_tree.png width=4in
    showTree(seqd_tree, format='plt')

    reordered_seqd, indices = reorderMatrix(labels, seqd_matrix, seqd_tree)
    plt.figure();
    @savefig ens_gnms_seqd_reordered_seqd_matrix.png width=4in
    showMatrix(reordered_seqd, origin='upper');

    plt.close('all')

.. ipython:: python

    rmsd_matrix = ens.getRMSDs(pairwise=True)
    @savefig ens_gnms_rmsd_matrix.png width=4in
    showMatrix(rmsd_matrix, origin='upper')

    plt.figure()
    rmsd_tree = calcTree(names=labels, 
                         distance_matrix=rmsd_matrix, 
                         method='upgma')
    @savefig ens_gnms_rmsd_tree.png width=4in
    showTree(rmsd_tree, format='plt')

    plt.figure()
    reordered_rmsd, indices = reorderMatrix(labels, rmsd_matrix, rmsd_tree)
    @savefig ens_gnms_rmsd_reordered_rmsd_matrix.png width=4in
    showMatrix(reordered_rmsd, origin='upper')

    plt.close('all')

Comparing sequence, structural and dynamic classifications
-------------------------------------------------------------------------------

We can reorder the seqd and sod matrices by the RMSD tree too to compare them:

.. ipython:: python

    reordered_seqd, indices = reorderMatrix(names=labels, matrix=seqd_matrix, tree=rmsd_tree)
    reordered_sod, indices = reorderMatrix(names=labels, matrix=so_matrix, tree=rmsd_tree)

.. ipython:: python

    @savefig ens_gnms_rmsd_reordered_seqd_matrix.png width=4in
    showMatrix(reordered_seqd, origin='upper')

    @savefig ens_gnms_rmsd_reordered_rmsd_matrix.png width=4in
    showMatrix(reordered_rmsd, origin='upper')

    @savefig ens_gnms_rmsd_reordered_sod_matrix.png width=4in
    showMatrix(arccos(reordered_sod), origin='upper')

    plt.close('all')

This analysis is quite sensitive to how many modes are used. As the number of modes approaches the full number, 
the dynamic distance order approaches the RMSD order. With smaller numbers, we see finer distinctions. This is 
particularly clear in the current case where we used just one mode.