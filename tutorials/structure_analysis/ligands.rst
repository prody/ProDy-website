.. _extract-ligands:

Ligand Extraction
=============================================================================

This example shows how to align structures of the same protein and extract
bound ligands from these structures.

:func:`.matchAlign` function can be used for aligning protein structures.
This example shows how to use it to extract ligands from multiple PDB
structures after superposing the structures onto a reference.
Output will be PDB files that contain ligands superposed onto the reference
structure.


Parse reference and blast search
-------------------------------------------------------------------------------

We start by importing everything from the ProDy package:

.. ipython:: python

   from prody import *
   from pylab import *
   ion()

First, we parse the reference structure and blast search PDB for similar
structure:

.. ipython:: python

   p38 = parsePDB('5uoj')
   seq = p38['A'].getSequence()

.. ipython::
   :verbatim:

   In [1]: blast_record = blastPDB(seq)

It is a good practice to save this record on disk, as NCBI may not respond to
repeated searches for the same sequence. We can do this using Python standard
library :mod:`pickle` as follows:

.. ipython:: python

   import pickle

Record is save using :func:`~pickle.dump` function into an open file:

.. ipython::
   :verbatim:

   In [10]: pickle.dump(blast_record, open('p38_blast_record.pkl', 'wb'))


Then, it can be loaded using :func:`~pickle.load` function:

.. ipython:: python

   blast_record = pickle.load(open('p38_blast_record.pkl', 'rb'))


Align structures and extract ligands
-------------------------------------------------------------------------------

Then, we parse the hits one-by-one, superpose them onto the reference
structure, and extract ligands:

.. ipython:: python

   for pdb_id in blast_record.getHits(90, 70):
       # blast search may return PDB identifiers of deprecated structures,
       # so we parse structures within a try statement
       try:
           pdb = parsePDB(pdb_id)
           pdb = matchAlign(pdb, p38)[0]
       except:
           continue
       else:
           ligand = pdb.select('not protein and not water')
           repr(ligand)
           if ligand:
               writePDB(pdb_id + '_ligand.pdb', ligand)

   !ls *_ligand.pdb

Ligands bound to p38 are outputted. Note that output PDB files may contain
multiple ligands.

The output can be loaded into a molecular visualization tool for analysis.

