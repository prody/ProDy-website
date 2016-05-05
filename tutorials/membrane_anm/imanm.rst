Implicit Membrane ANM
===============================================================================

Here we will make use of ProDy's implicit membrane ANM (imANM) capabilities to investigate the motions of a neurotransmitter transporter in the presence of the plasma membrane.  The procedure is based on the methods devescribed in [TL12]_ and relies on the RTB (Rotations and Translations of Blocks) method [FT00]_ of reducing complexity within ENMs. You will need following files:

  * `Membrane-aligned outward-facing structure file <2NWL-opm.pdb>`_
  * `Membrane-aligned inward-facing structure file <3KBC-opm.pdb>`_
  * `Outward-facing block definition file <2nwl_blocks.txt>`_
  * `Inward-facing block definition file <3kbc_blocks.txt>`_

Files in the following archives can be used to follow this tutorial:

  * `membrane ANM Tutorial Files (TGZ) <membanm_tutorial_files.tgz>`_
  * `membrane ANM Tutorial Files (ZIP) <membanm_tutorial_files.zip>`_

The first file contains the outward-facing structure of the glutamate transporter after insertion into the plasma membrane.  It is obtained from the `Orientations of Proteins in Membranes <http://opm.phar.umich.edu/>`_ database.

.. [TL12] Lezon TR, Bahar I. Constraints Imposed by the Membrane Selectively Guide the Alternating Access Dynamics of the Glutamate Transporter GltPh. *Biophys J* **2012** 102 1331-1340.

.. [FT00] Tama F, Gadea FJ, Marques O, Sanejouand YH. Building-block approach for determining low-frequency normal modes of macromolecules. *Proteins* **2000** 41 1-7.


Preparing the structures
-------------------------------------------------------------------------------
Begin by firing up ProDy in the usual manner:

.. ipython:: python

   from prody import *
   from pylab import *
   ion()

The imANM assumes that the membrane is normal to the z-axis, so it is important to use a structure that is properly aligned.  The structure from the `OPM <http://opm.phar.umich.edu/>`_ database will work.

.. ipython:: python

   of_all = parsePDB('2NWL-opm.pdb')  # Outward-facing structure
   if_all = parsePDB('3KBC-opm.pdb')  # Inward-facing structure


There will be warnings saying that ProDy wants to read beta factors, but the coordinates should be read properly.  In addition to atoms, the OPM file contains points to indicate the boundaries of the membrane. 

We now make two selections, one from each structure.  The selections are chosen so that the final structures are homotrimers with an equal number (398) of atoms in each subunit.  We also want to remove the three aspartate ligands, which are indicated as chain D in the outward-facing structure and have resid 500, in the inward-facing structure.

.. ipython:: python

   of_ca = of_all.select('protein and name CA and not (chain A and resid 119 to 122) and not (chain C and resid 119 to 123) and not chain D')
   if_ca = if_all.select('protein and name CA and not (resid 6 to 9) and not (resid 119 to 127) and resid < 500')

As a last step in preparation, we can align the structures so that we can calculate a deformation vector and compare the modes to it, as shown in the `ENM Tutorial`.

.. ipython:: python

   superpose(if_ca, of_ca)



Assigning Blocks
-------------------------------------------------------------------------------
ProDy's RTB method can be used for any system, whether or not a membrane is involved.  It just happens to be attached to the imANM part of ProDy.  The RTB method allows us to decompose the protein into pre-defined rigid blocks.  Atoms within a block do not move relative to each other (hence the descriptor "rigid"), but blocks can move relative to other blocks.  There are two main benefits of using blocks:  First, the Hessian for a good blocking scheme is smaller than the Hessian for an all-residue representation, so the modes can be calculated more quickly.  This is particularly useful when one is considering very large systems (currently, those containing thousands of residues).  Second, the use of rigid blocks reduces unphysical distortions of the structure, such as stretching of backbone bonds, that may result from the harmonic approximation.  These benefits come with the price of accuracy.  Imposing rigidity reduces the amount of dynamical detail that can be uncovered from the model.

In ProDy, a rigid block is defined as a set of atoms that co-move (i.e., the distances between them are fixed).  Typically the constituent atoms of a rigid block are spatially adjacent (i.e., they all belong to the same domain or secondary structure), but users are free to define blocks however they wish.  The only restriction is that a block cannot contain exactly two particles.  This restriction is in place because it is mathematically inconvenient to deal with two-particle blocks.  

We can either define blocks within our python session, or define them externally in a separate file and write a little bit of code to handle the tasks of reading the file and assigning residues to blocks.  This latter approach can be useful when exploring and comparing many different blocking schemes.  We have developed one such format for a `block file`, examples of which can be found in ``2nwl_blocks.txt`` and ``3kbc_blocks.txt``.  The first ten lines of ``2nwl_blocks.txt`` are::

    1 TYR A     10  VAL A     12
    4 LEU A     13  LYS A     15
    5 ILE A     16  TYR A     33
    6 GLY A     34  ALA A     36
    7 HIS A     37  VAL A     43
    8 LYS A     44  ALA A     70
    9 ALA A     71  ALA A     71
    10 SER A     72  SER A     72
    11 ILE A     73  ILE A     73
    12 SER A     74  LEU A     78


The columns, separated by whitespace, are formatted as follows:

      * *1.* Integer identifier of the block.
      * *2.* Three-letter code for first residue in block.
      * *3.* Chain ID of first residue in block.
      * *4.* Sequential number of first residue in block.
      * *5.* Three-letter code for last residue in block.
      * *6.* Chain ID of last residue in block.
      * *7.* Sequential number of last residue in block.

This is just one way of storing information on how the protein is deconstructed into blocks.  You can think of others.  We can read blocks from ``2nwl_blocks.txt`` into the array ``blocks`` as follows:

.. ipython:: python

   blk='2nwl_blocks.txt'
   with open(blk) as inp:
        for line in inp:
             b, n1, c1, r1, n2, c2, r2 = line.split()
             sel = of_ca.select('chain {} and resnum {} to {}'
                              .format(c1, r1, r2))
             if sel != None:
                sel.setBetas(b)


   of_blocks = of_ca.getBetas()

We will do the same for the blocks of the inward-facing structure.  The block definitions are based on secondary structures, which vary slightly between the structures.  We therefore have two separate blocking schemes.

.. ipython:: python

   blk='3kbc_blocks.txt'
   with open(blk) as inp:
        for line in inp:
             b, n1, c1, r1, n2, c2, r2 = line.split()
             sel = if_ca.select('chain {} and resnum {} to {}'
                              .format(c1, r1, r2))
             if sel != None:
                sel.setBetas(b)


   if_blocks = if_ca.getBetas()




Calculating the Modes
-------------------------------------------------------------------------------
To use the blocks in an RTB ANM calculation, we instantiate an RTB object for each structure:

.. ipython:: python

   of_rtb = RTB('2nwl')
   if_rtb = RTB('3kbc')

and we build a couple of Hessians using the coordinates of the crystal structures

.. ipython:: python

   of_coords = of_ca.getCoords()
   if_coords = if_ca.getCoords()
   of_rtb.buildHessian(of_coords, of_blocks, cutoff=11.0, scale=16., membrane_low=-1000.0, membrane_high=1000.0)
   if_rtb.buildHessian(if_coords, if_blocks, cutoff=11.0, scale=16., membrane_low=-1000.0, membrane_high=1000.0)

The scaling factor of 16 in this example means that the restoring force for any displacement in the x- or y-direction is 16 times greater than the force associated with a displacement in the z-direction.  The constraint on motions parallel to the membrane surface implicitly incorporates the membrane's effects into ANM.  To use RTB with no membrane effects, set ``scale=1.0`` (which is also the default value).  We have here set the boundaries of the membrane to extend well beyond the protein, effectively applying the implicit membrane scaling to the entire protein.

Now we calculate the modes and write them to a pair of .nmd files for viewing.

.. ipython:: python

   of_rtb.calcModes()
   if_rtb.calcModes()
   writeNMD('2nwl_im.nmd',of_rtb,of_ca.select('protein and name CA'))
   writeNMD('3kbc_im.nmd',if_rtb,if_ca.select('protein and name CA'))


.. figure:: _static/figures/membrane_anm-imanm_of3.png
   :scale: 100%

The third mode of the outward-facing structure moves all three transport domains simultaneously through the membrane in a 'lift-like' motion.

.. figure:: _static/figures/membrane_anm-imanm_if6.png
   :scale: 100%

A similar motion is shown in mode 6 of the inward-facing structure.
