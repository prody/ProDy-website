.. _signdy-data:

Data Collection
===============================================================================

The first step in signature dynamics analysis is to collect a set of related 
protein structures. This can be achieved by multiple routes, belonging to three 
main types: a query search of the PDB using BLAST or Dali, extraction of PDB IDs 
from the Pfam or CATH database, or input of a pre-defined list. We demonstrate a 
Dali search and a Pfam data extraction here.

First, we make necessary imports from ProDy and Matplotlib packages.

.. ipython:: python

   from prody import *
   from pylab import *
   ion()

Dali
-------------------------------------------------------------------------------

Next, we use Dali functions to search for type-I periplasmic binding protein domains, 
a group of protein structures originally found in bacteria for transport of solutes 
across the periplasmic space and later seen in various eukaryotic receptors including 
ionotropic and metabotropic glutamate receptors. We use the N-terminal domain of AMPA
receptor subunit GluA2 as a query.

.. ipython:: python

   dali_rec = searchDali('3H5V','A')
   ens = dali_rec.buildDaliEnsemble()

   