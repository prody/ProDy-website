.. _signdy-data:

Data Collection
===============================================================================

The first step in signature dynamics analysis is to collect a set of related 
protein structures and build a :class:`.PDBEnsemble`. This can be achieved by 
multiple routes: a query search of the PDB using :func:`.blastPDB` or :func:`.searchDali`, 
extraction of PDB IDs from the Pfam or CATH database, or input of a pre-defined list. 

We demonstrate the Dali method here. The Pfam and CATH methods are still under development 
and the function :func:`.blastPDB` is described in the `Structure Analysis Tutorial`_.

We apply these methods to the type-I periplasmic binding protein domains, 
a group of protein structures originally found in bacteria for transport of solutes 
across the periplasmic space and later seen in various eukaryotic receptors including 
ionotropic and metabotropic glutamate receptors. We use the N-terminal domain of AMPA
receptor subunit GluA2 (gene name GRIA2) as a query.

First, make necessary imports from ProDy and Matplotlib packages if you haven't already.

.. ipython:: python

    from prody import *
    from pylab import *
    ion()

Searching the PDB with Dali and building a PDBEnsemble
-------------------------------------------------------------------------------

First we use the function searchDali to search the PDB, which returns a :class:`.DaliRecord` 
that contains a list of PDB IDs and their corresponding mappings to the reference structure. 

.. ipython:: python

    dali_rec = searchDali('3H5V','A')
    dali_rec

Next, we get the lists of PDB IDs and mappings from *dali_rec*, parse the *pdb_ids* to get 
a list of :class:`.AtomGroup` instances, and provide them together with mappings to 
:func:`.buildPDBEnsemble`. We provide the keyword argument ``seqid=10`` to account for the 
low sequence identity between some of the structures.

.. ipython:: python

    pdb_ids = dali_rec.getPDBs()
    mappings = dali_rec.getMappings()

.. ipython:: python

    pdbs = parsePDB(*pdb_ids, subset='ca')
    len(pdbs)

.. ipython:: python

    dali_ens = buildPDBEnsemble(pdbs, mapping=mappings, seqid=20)
    dali_ens

Finally we save the ensemble for later processing:

.. ipython:: python

   saveEnsemble(dali_ens, 'dali_ensemble')



.. _`Structure Analysis Tutorial`: http://prody.csb.pitt.edu/tutorials/structure_analysis/blastpdb.html
.. _`list_comprehensions`: https://docs.python.org/2/tutorial/datastructures.html#list-comprehensions