.. _cavitracer_single:

Detection of channels across heterogeneous structures
===============================================================================


Now, we will illustrate how to detect intraprotein channels across various PDB
structures. As an example, we will use Cytochrome P450. 

First, we will provide a list of PDB codes with different Cytochrome P450
structures. Such a list can be also provided by BLAST, Dali, or Foldseek
(see :func:`.runBLAST`, :func:`.runDali`, :func:`.runFoldseek` in the
`InSty tutorial`_).

.. ipython:: python
   :verbatim:

   pdb_files = ["1TQN", "1W0E", "4I3Q", "5A1P", "5VCC", "6BD6", "6BD8", "6BDI", 
            "6BDM", "6DA8", "6DAJ", "6DAL", "6MA6", "6MA7", "6OOA", "6UNE", "6UNG"]


Protein preparation
-------------------------------------------------------------------------------

Before performing the analysis, we will align all the structures onto the first
structure from our list (``target``) and save it under a new name with ``align__``
prefix. Such an approach was also shown in other ProDy tutorials and explained
in detail (see Structure Composition of the `Structure Analysis tutorial`_).

.. ipython:: python
   :verbatim:

   structures = parsePDB(pdb_files)
   target = structures[0]
   rmsds = []

   for mobile in structures[1:]:
       try:
           i = mobile.getTitle()
           print (i)
           matches = matchChains(mobile.protein, target.protein, subset='bb')
           m = matches[0]
           m0_alg, T = superpose(m[0], m[1], weights=m[0].getFlags("mapped"))
           rmsds.append(calcRMSD(m[0], m[1], weights=m[0].getFlags("mapped")))
           writePDB('align__'+i+'.pdb', mobile)
       except: pass

.. parsed-literal::

   @> 17 PDBs were parsed in 39.35s.
   @> Checking AtomGroup 1W0E: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 1W0E (len=450) and Chain A from 1TQN (len=468):
   @> 	Match: 447 residues match with 99% sequence identity and 96% overlap.
   
   1W0E
   4I3Q
   
   @> Checking AtomGroup 4I3Q: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 4I3Q (len=465) and Chain A from 1TQN (len=468):
   @> 	Match: 465 residues match with 100% sequence identity and 99% overlap.
   @> Checking AtomGroup 5A1P: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 5A1P (len=465) and Chain A from 1TQN (len=468):
   @> 	Match: 463 residues match with 100% sequence identity and 99% overlap.
   
   5A1P
   5VCC
   
   @> Checking AtomGroup 5VCC: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 5VCC (len=457) and Chain A from 1TQN (len=468):
   @> 	Match: 457 residues match with 100% sequence identity and 98% overlap.
   @> Checking AtomGroup 6BD6: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6BD6 (len=457) and Chain A from 1TQN (len=468):
   @> 	Match: 457 residues match with 100% sequence identity and 98% overlap.

   6BD6
   6BD8

   @> Checking AtomGroup 6BD8: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6BD8 (len=457) and Chain A from 1TQN (len=468):
   @> 	Match: 457 residues match with 100% sequence identity and 98% overlap.
   @> Checking AtomGroup 6BDI: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6BDI (len=462) and Chain A from 1TQN (len=468):
   @> 	Match: 461 residues match with 100% sequence identity and 99% overlap.

   6BDI
   6BDM

   @> Checking AtomGroup 6BDM: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6BDM (len=456) and Chain A from 1TQN (len=468):
   @> 	Match: 456 residues match with 100% sequence identity and 97% overlap.
   @> Checking AtomGroup 6DA8: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6DA8 (len=451) and Chain A from 1TQN (len=468):
   @> 	Match: 450 residues match with 100% sequence identity and 96% overlap.
   @> Checking AtomGroup 6DAJ: 1 chains are identified

   6DA8
   6DAJ

   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6DAJ (len=451) and Chain A from 1TQN (len=468):
   @> 	Match: 450 residues match with 100% sequence identity and 96% overlap.
   @> Checking AtomGroup 6DAL: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6DAL (len=454) and Chain A from 1TQN (len=468):
   @> 	Match: 453 residues match with 100% sequence identity and 97% overlap.
   @> Checking AtomGroup 6MA6: 1 chains are identified

   6DAL
   6MA6

   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6MA6 (len=466) and Chain A from 1TQN (len=468):
   @> 	Match: 464 residues match with 100% sequence identity and 99% overlap.
   @> Checking AtomGroup 6MA7: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6MA7 (len=463) and Chain A from 1TQN (len=468):
   @> 	Match: 461 residues match with 99% sequence identity and 99% overlap.
   @> Checking AtomGroup 6OOA: 1 chains are identified

   6MA7
   6OOA

   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6OOA (len=452) and Chain A from 1TQN (len=468):
   @> 	Match: 452 residues match with 100% sequence identity and 97% overlap.
   @> Checking AtomGroup 6UNE: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6UNE (len=452) and Chain A from 1TQN (len=468):
   @> 	Match: 452 residues match with 100% sequence identity and 97% overlap.
   @> Checking AtomGroup 6UNG: 1 chains are identified

   6UNE
   6UNG
   ..
   ..


The names of the structurally aligned PDBs are further uploaded into ``pdb_files_new``.
This operation is required if we want to compare the localization of channels in
various PDB models.

.. ipython:: python
   :verbatim:

   from pathlib import Path
   pdb_files_new = [p.name for p in Path(".").glob("align__*.pdb")]
   pdb_files_new

.. parsed-literal::

   ['align__6OOA.pdb',
    'align__6BD8.pdb',
    'align__6BD6.pdb',
    'align__6BDI.pdb',
    'align__6DAJ.pdb',
    'align__5A1P.pdb',
    'align__6BDM.pdb',
    'align__6DA8.pdb',
    'align__4I3Q.pdb',
    'align__6UNG.pdb',
    'align__6UNE.pdb',
    'align__6DAL.pdb',
    'align__6MA6.pdb',
    'align__1W0E.pdb',
    'align__6MA7.pdb',
    'align__5VCC.pdb']


Channel prediction for individual PDB files
-------------------------------------------------------------------------------

Now, all PDB structures can be analyzed to detect channels, tunnels, or
pores in the protein structure using :func:`.calcChannels`. In this example,
the results will be saved in a single file (all detected channels in ``.pqr``
file with the name of the input file) and in multiple files (``separate``
must be set to ``True``) to save each detected channel independently.
Additionally, we will use :func:`.getChannelParameters` and
:func:`.getChannelResidueNames` to obtain information about channel
parameters and residues involved in its formation. To save this information
in the local directory, we provide ``param_file_name`` and
``residues_file_name``. The function will create files with the same name as
provided PDB with ``'_Parameters_All_channels.txt'`` and
``'_Residues_All_channels.txt'`` as suffixes.

.. ipython:: python
   :verbatim:

   for i in pdb_files_new:
       base_name = pdb_files_new[0]
       atoms = parsePDB(i).select('protein')
       channels2, surface2 = calcChannels(atoms, output_path=i[:-4], separate=True)
       getChannelParameters(channels2, param_file_name=i[:-4])
       getChannelResidueNames(atoms, channels2, residues_file_name=i[:-4])

.. parsed-literal::

   @> 3727 atoms and 1 coordinate set(s) were parsed in 0.12s.
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> WARNING structure has no hydrogens and r2=0.90 is below 1.2 A: the space left by the missing H is then wide enough for the probe to pass, and channels will be found through interstices that do not exist in the real protein (their number can rise several-fold). Either add hydrogens, or raise r2 to 1.2 A or more, where protonated and unprotonated structures give the same channels.
   @> Substituted 3630 atoms with 22830 homogeneous balls of radius 1.52 A in 0.20s.
   @> Delaunay tessellation of 22830 points constructed in 0.87s.
   @> Surface and inner simplices filtered in 1.85s.
   @> 10 surface cavities detected and filtered in 0.34s.
   @> Channel pathfinding (graph Dijkstra) over 10 cavities completed in 1.83s.
   @> Detected 35 channels.
   @> Saving multiple results to directory ..
   @> Channel calculation completed in 5.32s.
   @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> channel 0: 	45.24 		6.35 		0.92
   @> channel 1: 	45.17 		6.41 		0.91
   @> channel 2: 	91.45 		10.27 		0.96
   @> channel 3: 	27.86 		6.22 		0.93
   @> channel 4: 	30.59 		6.27 		0.84
   @> channel 5: 	36.62 		7.4 		0.84
   @> channel 6: 	68.78 		10.84 		0.96
   @> channel 7: 	96.88 		11.03 		0.94
   @> channel 8: 	31.08 		7.2 		0.88
   @> channel 9: 	82.38 		12.93 		0.96
   @> channel 10: 	83.38 		12.99 		0.95
   @> channel 11: 	28.95 		7.93 		0.9
   @> channel 12: 	55.61 		12.24 		0.91
   @> channel 13: 	45.36 		11.92 		0.87
   @> channel 14: 	153.25 		23.13 		0.86
   @> channel 15: 	158.78 		24.83 		0.86
   @> channel 16: 	92.84 		21.61 		0.9
   @> channel 17: 	130.74 		28.38 		0.92
   @> channel 18: 	121.29 		28.43 		0.92
   @> channel 19: 	193.41 		34.83 		0.9
   @> channel 20: 	163.44 		34.56 		0.89
   @> channel 21: 	310.39 		44.99 		0.81
   @> channel 22: 	416.28 		48.01 		0.81
   @> channel 23: 	385.01 		47.7 		0.81
   @> channel 24: 	254.88 		47.38 		0.83
   @> channel 25: 	318.51 		49.41 		0.81
   @> channel 26: 	250.6 		50.89 		0.57
   @> channel 27: 	261.86 		51.71 		0.57
   @> channel 28: 	354.45 		61.17 		0.81
   @> channel 29: 	377.43 		63.37 		0.81
   @> channel 30: 	350.68 		63.92 		0.81
   @> channel 31: 	370.19 		66.39 		0.81
   @> channel 32: 	359.68 		66.19 		0.81
   @> channel 33: 	385.8 		71.46 		0.81
   @> channel 34: 	388.4 		72.31 		0.81
   @> 3695 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 3690 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 3710 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3700 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3700 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3695 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 3710 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3705 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 3690 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 3725 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 3730 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 3725 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 3735 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 3755 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 3815 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 3835 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 3855 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 3905 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 3910 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 3905 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 3905 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 3975 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 4010 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 3990 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4020 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4040 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4070 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 4065 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 4110 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4125 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4125 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4145 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4135 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4205 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4220 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> Channel residues were saved to: align__6OOA_Residues_All_channels.txt
   @> 3756 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> 3756 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> WARNING structure has no hydrogens and r2=0.90 is below 1.2 A: the space left by the missing H is then wide enough for the probe to pass, and channels will be found through interstices that do not exist in the real protein (their number can rise several-fold). Either add hydrogens, or raise r2 to 1.2 A or more, where protonated and unprotonated structures give the same channels.
   @> Substituted 3675 atoms with 23131 homogeneous balls of radius 1.52 A in 0.19s.
   @> Delaunay tessellation of 23131 points constructed in 0.85s.
   @> Surface and inner simplices filtered in 1.70s.
   @> 8 surface cavities detected and filtered in 0.32s.
   @> Channel pathfinding (graph Dijkstra) over 8 cavities completed in 4.16s.
   @> Detected 70 channels.
   @> Saving multiple results to directory ..
   @> Channel calculation completed in 7.38s.
   @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> channel 0: 	53.33 		5.65 		0.85
   @> channel 1: 	36.47 		7.03 		0.94
   @> channel 2: 	42.3 		6.7 		0.84
   @> channel 3: 	38.71 		7.09 		0.91
   @> channel 4: 	31.3 		8.82 		0.86
   ..
   ..


Selection of channels in a certain protein area
-------------------------------------------------------------------------------

CaviTracer will generate multiple channels, but we might be interested only
in the channels that are localized in a certain region. For that reason, we
created :func:`.selectChannelBySelection` which can extract channels that
are localized near the region of our interest. Below, we first select which
``pqr_files`` files we want to analyze (to exclude all channels in one file
that is also saved by default). Next, we use the ``residue_sele`` option to
apply any selection that is understandable by ProDy select. In our case, we
select all ``FIL`` atoms (channel prediction artificial atoms) that are
within 5 Angstroms from the residue with the number 442. These values can be
changed using ``distA`` parameter of the function.

.. ipython:: python
   :verbatim:

   from pathlib import Path
   pdb_files_new_channels = [i.name for i in Path(".").glob("align__*chl*.pqr") if i.is_file()]
   pdb_files_new_channels

.. parsed-literal::

   [['align__6UNE_chl18.pqr',
     'align__6BDM_chl53.pqr',
     'align__6BDM_chl65.pqr',
     'align__6DAJ_chl16.pqr',
     'align__6BDI_chl14.pqr',
     'align__5A1P_chl29.pqr',
     'align__4I3Q_chl27.pqr',
     'align__6MA7_chl39.pqr',
     'align__6BDM_chl31.pqr',
     'align__6BDI_chl40.pqr',
     'align__5VCC_chl2.pqr',
     'align__5A1P_chl38.pqr',
     'align__6DAL_chl8.pqr',
     'align__1W0E_chl62.pqr',
     'align__6MA6_chl15.pqr',
     'align__6MA7_chl30.pqr',
     'align__6UNE_chl43.pqr',
     'align__6BD8_chl19.pqr',
     'align__6OOA_chl33.pqr',
     'align__5A1P_chl34.pqr',
     'align__6BD8_chl24.pqr',
     'align__6UNG_chl54.pqr',
     'align__6MA6_chl23.pqr',
     'align__6MA7_chl57.pqr',
     'align__1W0E_chl4.pqr',
     ..
     ..
     'align__6MA6_chl48.pqr',
     'align__1W0E_chl44.pqr',
     'align__6MA7_chl14.pqr',
     'align__5A1P_chl47.pqr',
     'align__6BDI_chl34.pqr',
     'align__1W0E_chl49.pqr',
     'align__6BDI_chl33.pqr',
     'align__6MA6_chl17.pqr']


.. ipython:: python
   :verbatim:

   atoms = parsePDB(pdb_files_new_channels[0].split('_ch')[0]+'.pdb')
   selectChannelBySelection(atoms, pqr_files=pdb_files_new_channels, residue_sele='resid 442')

.. parsed-literal::

   @> 3718 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 335 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 640 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 875 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 370 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 420 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: selected_files
   @> 345 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 515 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: selected_files
   @> 525 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: selected_files
   @> 685 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: selected_files
   @> 475 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 95 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 390 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 280 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 450 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: selected_files
   @> 345 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 535 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: selected_files
   @> 435 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 350 atoms and 1 coordinate sets were parsed in 0.00s.
   ..
   ..
   @> Filtered files are now in: selected_files
   @> 475 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: selected_files
   @> 470 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: selected_files
   @> 515 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: selected_files
   @> 405 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: selected_files
   @> 405 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 355 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Selected files: 
   @> align__6BDI_chl14.pqr align__4I3Q_chl27.pqr align__6MA7_chl39.pqr 
      align__6BDM_chl31.pqr align__1W0E_chl62.pqr align__6MA7_chl30.pqr 
      align__6BD8_chl19.pqr align__5A1P_chl34.pqr align__6MA7_chl52.pqr 
      align__6BDM_chl30.pqr align__6MA6_chl31.pqr align__6BD8_chl26.pqr 
      align__1W0E_chl24.pqr align__6BDI_chl22.pqr align__6BD6_chl48.pqr 
      align__6MA7_chl15.pqr align__6BDI_chl37.pqr align__6UNE_chl6.pqr 
      align__6BD6_chl54.pqr align__4I3Q_chl35.pqr align__6BD6_chl71.pqr 
      align__6UNE_chl42.pqr align__6BDM_chl28.pqr align__5A1P_chl23.pqr 
      align__1W0E_chl64.pqr align__5VCC_chl17.pqr align__1W0E_chl39.pqr 
      align__6BD8_chl43.pqr align__6BDM_chl29.pqr align__6UNE_chl7.pqr 
      align__6UNG_chl16.pqr align__6UNE_chl55.pqr align__1W0E_chl22.pqr 
      align__6UNE_chl21.pqr align__6BDI_chl57.pqr align__6BDI_chl47.pqr 
      align__6UNE_chl41.pqr align__1W0E_chl38.pqr align__5A1P_chl17.pqr 
      align__6BDI_chl43.pqr align__1W0E_chl57.pqr align__6MA7_chl34.pqr 
      align__6UNE_chl29.pqr align__6UNE_chl49.pqr align__6BD6_chl53.pqr 
      align__6UNG_chl44.pqr align__6BDI_chl69.pqr align__6BDI_chl17.pqr 
      align__4I3Q_chl58.pqr align__6UNE_chl48.pqr align__6BD6_chl42.pqr 
      align__6BDM_chl34.pqr align__6BD6_chl24.pqr align__6BDI_chl20.pqr 
      align__6BD8_chl42.pqr align__6UNE_chl39.pqr align__6BD6_chl21.pqr 
      align__6UNE_chl72.pqr align__6BD6_chl63.pqr align__1W0E_chl9.pqr 
      align__6UNE_chl56.pqr align__4I3Q_chl30.pqr align__5A1P_chl46.pqr 
      align__4I3Q_chl26.pqr align__5A1P_chl18.pqr align__1W0E_chl68.pqr 
      align__6BD6_chl12.pqr align__6BD8_chl50.pqr align__1W0E_chl51.pqr 
      align__5VCC_chl35.pqr align__4I3Q_chl69.pqr align__6BDI_chl68.pqr 
      align__6MA6_chl27.pqr align__5A1P_chl55.pqr align__6BDI_chl62.pqr 
      align__6BDM_chl47.pqr align__6BDM_chl45.pqr align__6MA6_chl49.pqr 
      align__6MA7_chl61.pqr align__6BDM_chl57.pqr align__6UNG_chl31.pqr 
      align__5A1P_chl40.pqr align__5VCC_chl32.pqr align__5VCC_chl36.pqr 
      align__4I3Q_chl61.pqr align__6BD6_chl58.pqr align__6MA7_chl43.pqr 
      align__6MA7_chl41.pqr align__6UNG_chl24.pqr align__6MA6_chl25.pqr 
      align__1W0E_chl59.pqr align__6UNG_chl62.pqr align__1W0E_chl61.pqr 
      align__5A1P_chl22.pqr align__6UNE_chl71.pqr align__6MA6_chl29.pqr 
      ..
      ..
      align__5A1P_chl47.pqr align__6BDI_chl34.pqr align__1W0E_chl49.pqr
   @> If newly created files are empty please check whether the parameter 
      names are: PDB_id+_Parameters_All_channels.txt


If we do not provide ``folder_name``, the results will be copied into the folder
``selected_file``, which is created automatically. That name can be changed using
``folder_name``, as shown below. Additionally, if we generated in the previous
step files with parameters and residues for the channels, we can also use
options ``residues_file`` and ``param_file`` set to ``True``. Then two new
files will be created **Selected_channel_residues.txt** and
**Selected_channel_parameters.txt**. In those files, we will find data for
selected channels only. 

.. ipython:: python
   :verbatim:

   atoms = parsePDB(pdb_files_new_channels[0].split('_ch')[0]+'.pdb')
   selectChannelBySelection(atoms, pqr_files=pdb_files_new_channels, residue_sele='resid 442',
   folder_name='res442', residues_file=True, param_file=True)

.. parsed-literal::

   @> 3718 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 335 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 640 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 875 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 370 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 420 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: res442
   @> 345 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 515 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: res442
   @> 525 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: res442
   @> 685 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: res442
   @> 475 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 95 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 390 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 280 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 450 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: res442
   @> 345 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 535 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: res442
   @> 435 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 350 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: res442
   @> 575 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 445 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: res442
   @> 220 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 455 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 490 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 450 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 75 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 660 atoms and 1 coordinate sets were parsed in 0.00s.
   ..
   ..
   @> Selected files: 
   @> align__6BDI_chl14.pqr align__4I3Q_chl27.pqr align__6MA7_chl39.pqr 
      align__6BDM_chl31.pqr align__1W0E_chl62.pqr align__6MA7_chl30.pqr 
      align__6BD8_chl19.pqr align__5A1P_chl34.pqr align__6MA7_chl52.pqr 
      align__6BDM_chl30.pqr align__6MA6_chl31.pqr align__6BD8_chl26.pqr 
      align__1W0E_chl24.pqr align__6BDI_chl22.pqr align__6BD6_chl48.pqr 
      align__6MA7_chl15.pqr align__6BDI_chl37.pqr align__6UNE_chl6.pqr 
      align__6BD6_chl54.pqr align__4I3Q_chl35.pqr align__6BD6_chl71.pqr 
      ..
      ..
      align__5A1P_chl47.pqr align__6BDI_chl34.pqr align__1W0E_chl49.pqr


Calculating overlapping channel regions across PDB structures
-------------------------------------------------------------------------------

Now, we will select files for analysis. In our case, we will select files with
all predicted channels in one file that starts with ``'align__'``. 

.. ipython:: python
   :verbatim:

   from pathlib import Path
   pdb_files_new_channels_ALL = [i.name for i in Path(".").glob("align__????.pqr") if i.is_file()]
   pdb_files_new_channels_ALL

.. parsed-literal::

   ['align__4I3Q.pqr',
    'align__6BDM.pqr',
    'align__6UNG.pqr',
    'align__6OOA.pqr',
    'align__6MA7.pqr',
    'align__6DAJ.pqr',
    'align__5A1P.pqr',
    'align__1W0E.pqr',
    'align__6BD6.pqr',
    'align__6UNE.pqr',
    'align__5VCC.pqr',
    'align__6BD8.pqr',
    'align__6BDI.pqr',
    'align__6DAL.pqr',
    'align__6MA6.pqr',
    'align__6DA8.pqr']   


To compute the so-called overlapping surface of the detected channels, we
should use :func:`.calcChannelSurfaceOverlaps`. To save the results, we
need to specify ``output_file_name``.

.. ipython:: python
   :verbatim:

   calcChannelSurfaceOverlaps(pqr_files=pdb_files_new_channels_ALL, output_file_name='overlapping_surf.pdb')

.. parsed-literal::

   [@> Number of PQR files: 16
    @> Resolution: 0.5
    @> max_proc: 2
    @> Calculating overlaps using 2 processes.
    @> 23365 atoms and 1 coordinate sets were parsed in 0.20s.
    @> 27720 atoms and 1 coordinate sets were parsed in 0.22s.
    @> 41040 atoms and 1 coordinate sets were parsed in 0.20s.
    @> 9280 atoms and 1 coordinate sets were parsed in 0.04s.
    @> 26615 atoms and 1 coordinate sets were parsed in 0.15s.
    @> 8530 atoms and 1 coordinate sets were parsed in 0.04s.
    @> 19370 atoms and 1 coordinate sets were parsed in 0.11s.
    @> 35700 atoms and 1 coordinate sets were parsed in 0.19s.
    @> 26080 atoms and 1 coordinate sets were parsed in 0.12s.
    @> 32005 atoms and 1 coordinate sets were parsed in 0.15s.
    @> 18525 atoms and 1 coordinate sets were parsed in 0.10s.
    @> 25105 atoms and 1 coordinate sets were parsed in 0.12s.
   [@> 31010 atoms and 1 coordinate sets were parsed in 0.16s.
    @> 11550 atoms and 1 coordinate sets were parsed in 0.06s.
    @> 35110 atoms and 1 coordinate sets were parsed in 0.19s.
    @> 380 atoms and 1 coordinate sets were parsed in 0.00s.
    @> Overlap written to: overlapping_surf.pdb
    @> Number of occupied overlap voxels: 143228

    'overlapping_surf.pdb'


The results can be displayed in any graphical visualization program. We
will use VMD_. Generated file ``overlapping_surf.pdb`` will contain a grid
with points and their corresponding values that will describe the occupation
of the channel in a specific space point across various structures that
were provided by ``pqr_files``.

Below, we can see a visualization of the grid that was created by
:func:`.calcChannelSurfaceOverlaps`. Values in the ``Occupancy`` column
correspond to the occupation of each point in ``pqr_files`` files. Values
range from 0 to 1, where ``1`` means that all ``pqr_files`` files have a
channel at this particular point. 

Below, we display the outcome of the prediction.


.. figure:: images/cavitracer_figure9.jpg
   :scale: 50 %


To display the most significant information about channel occupancy, it is
good to display values that are higher than 0.8 in the ``Occupancy``
column. In such a case, we will see only minimal channels that occur in at
least 80% of analyzed files from ``pqr_files``.

.. figure:: images/cavitracer_figure10.jpg
   :scale: 50 %


.. _InSty tutorial: http://www.bahargroup.org/prody/tutorials/insty_tutorial/
.. _Structure Analysis tutorial: http://www.bahargroup.org/prody/tutorials/structure_analysis/