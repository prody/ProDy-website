.. _signdy-data:

Data Collection
===============================================================================

The first step in signature dynamics analysis is to collect a set of related 
protein structures. This can be achieved by multiple routes, belonging to three 
main types: a query search of the PDB using BLAST or Dali, extraction of PDB IDs 
from the Pfam or CATH database, or input of a pre-defined list. We demonstrate the 
Dali, Pfam and CATH methods here. The function :func:`blastPDB` is described in 
the structure_analysis_ tutorial.

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

Dali
-------------------------------------------------------------------------------

First we use the function searchDali to search the PDB, which returns a DaliRecord 
that contains a list of PDBs and their corresponding alignments to the reference structure. 
This can be used to build a PDB ensemble directly using the method associated to this class.

Dali 

.. ipython:: python

    dali_rec = searchDali('3H5V','A')
    dali_ens = dali_rec.buildDaliEnsemble()


Pfam
-------------------------------------------------------------------------------

Alternatively, we can Pfam functions to parse PDBs belonging to a particular family 
and use the associated MSA to build a PDB ensemble using them. The first step is to 
find the appropriate Pfam family by searching with the Uniprot code:

.. ipython:: python
    
    searchPfam('GRIA2_RAT')


We see three domains of which the first one in the sequence (PF01094) is the one of interest 
and we parse the corresponding MSA and PDBs. 

.. ipython:: python

    fetchPfamMSA('PF01094')
    pfam_msa = parseMSA('PF01094_full.sth')
    pfam_pdbs = parsePfamPDBs(query='GRIA2_RAT', subset='ca')


We then extract the individual chains and build a PDB ensemble from them.

.. ipython:: python

    pfam_chains = []
    for pdb in pfam_pdbs:
        for chain in pdb.getHierView():
            pfam_chains.append(chain)

    pfam_ens = buildPDBEnsemble(pfam_chains, alignments=pfam_msa)


CATH
-------------------------------------------------------------------------------




.. _structure_analysis: http://prody.csb.pitt.edu/tutorials/structure_analysis/blastpdb.html
