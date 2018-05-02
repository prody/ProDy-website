Classification using matrices and trees of sequence, structure and dynamics distances
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


Spectral overlap and distance
-------------------------------------------------------------------------------

We calculate the spectral overlap matrix, calculate a tree from its arccosine and 
reorder the spectral overlap matrix using the tree as follows: 

.. ipython:: python

    so_matrix = calcEnsembleSpectralOverlaps(mode_ens[:,1)
    so_tree = calcTree(names=mode_ens.getLabels(), distance_matrix=arccos(so_matrix), method='upgma')
    reordered_so, new_so_indices = reorderMatrix(so_matrix, so_tree, names=new_ens.getLabels())


We can show the original and reordered spectral distance matrices and the tree as follows.

.. ipython:: python

    show = showMatrix(arccos(so_matrix))
    plt.figure()
    showTree(so_tree, format='networkx', label_colors=colors_dict)
    plt.figure()
    show = showMatrix(arccos(reordered_so))


Sequence and structural distances
-------------------------------------------------------------------------------

The sequence distance is given by the Hamming distance, which is calculated by 
subtracting the percentage identity (fraction) from 1, and the structural distance 
is the RMSD. We can also calculate and show the matrices and trees for these from 
the PDB ensemble.

.. ipython:: python

    rmsd_matrix = dali_ens.getRMSDs(pairwise=True)
    rmsd_tree = calcTree(names=dali_ens.getLabels(), distance_matrix=rmsd_matrix, method='upgma')

    seqid_matrix = buildSeqidMatrix(dali_ens.getMSA)
    seqd_matrix = 1. - seqid_matrix
    seqd_tree = calcTree(names=dali_ens.getLabels(), distance_matrix=seqd_matrix, method='upgma')

Comparing sequence, structural and dynamic classifications
-------------------------------------------------------------------------------

We can reorder all these matrices by the RMSD tree to compare them:

.. ipython:: python

    reordered_seqd, new_seqd_indices = reorderMatrix(seqd_matrix, rmsd_tree, names=new_ens.getLabels())
    reordered_rmsd, new_rmsd_indices = reorderMatrix(rmsd_matrix, rmsd_tree, names=new_ens.getLabels())
    reordered_sod, new_sod_indices = reorderMatrix(so_matrix, rmsd_tree, names=new_ens.getLabels())

    show = showMatrix(arccos(reordered_seqd))
    plt.figure()
    show = showMatrix(arccos(reordered_rmsd))
    plt.figure()
    show = showMatrix(arccos(reordered_sod))


This analysis is quite sensitive to how many modes are used. As the number of modes approaches the full number, 
the dynamic distance order approaches the RMSD order. With smaller numbers, we see finer distinctions. This is 
particularly clear in the current case where we used just one mode.