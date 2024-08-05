.. _watfinder_prep:

Protein Preparation
===============================================================================

Since PDB structures often lack hydrogen atoms in water molecules, we need to
add them. We can use the :func:`.addMissingAtoms` function. This function utilizes
either Openbabel_ or PDBFixer_. Those are external packages; therefore, they should be
installed independently (see Recommended Programs). 

If we prefer not to install additional tools, we can alternatively provide
a PDB structure with hydrogens already added by other software.

Here, we will fetch a structure of LMW-PTP (PDB: **5KQM**) from 
Protein Data Bank (PDB) in uncompressed form using ``compressed=False`` 
and add missing atoms using :func:`.addMissingAtoms`:

.. _Openbabel: https://github.com/openbabel
.. _PDBfixer: https://github.com/openmm/pdbfixer


.. ipython:: python
   :verbatim:

   pdb = '5kqm'
   filename = fetchPDB(pdb, compressed=False)

.. parsed-literal::

   @> Connecting wwPDB FTP server RCSB PDB (USA).
   @> 5kqm downloaded (5kqm.pdb)
   @> PDB download via FTP completed (1 downloaded, 0 failed).

.. ipython:: python
   :verbatim:

   filename2 = addMissingAtoms(filename)

.. parsed-literal::

   @> Hydrogens were added to the structure. 
   Structure addH_5kqm.pdb is saved in the local directry.

A new file containing hydrogens is now generated, prefixed with ``'addH_'``. For protein structures lacking
hydrogens, results will also be computed without applying angle criteria.

