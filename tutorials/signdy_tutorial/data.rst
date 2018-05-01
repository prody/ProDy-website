.. _signdy-data:

Data Collection
===============================================================================

The first step in signature dynamics analysis is to collect a set of related 
protein structures. This can be achieved by multiple routes, belonging to three 
main types: a query search of the PDB using BLAST or Dali, extraction of PDB IDs 
from the Pfam or CATH database, or input of a pre-defined list. We demonstrate the 
Dali, Pfam and CATH methods here. The function :func:`blastPDB` is described in 
the structure_analysis_ tutorial.

First, we make necessary imports from ProDy and Matplotlib packages.

.. ipython:: python

   from prody import *
   from pylab import *
   ion()

Dali
-------------------------------------------------------------------------------

We can use Dali functions to search for type-I periplasmic binding protein domains, 
a group of protein structures originally found in bacteria for transport of solutes 
across the periplasmic space and later seen in various eukaryotic receptors including 
ionotropic and metabotropic glutamate receptors. We use the N-terminal domain of AMPA
receptor subunit GluA2 as a query.

.. ipython:: python

   dali_rec = searchDali('3H5V','A')
   dali_ens = dali_rec.buildDaliEnsemble()


Pfam
-------------------------------------------------------------------------------

Alternatively, we can Pfam functions to parse PDBs belonging to a particular family 
and build a PDB ensemble as usual. 

.. ipython:: python

   pdbs = parsePfamPDBs(query='GRIA2_RAT', end=380, subset='ca')
   

.. _structure_analysis: http://prody.csb.pitt.edu/tutorials/structure_analysis/blastpdb.html
