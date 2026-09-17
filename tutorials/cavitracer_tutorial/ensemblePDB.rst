.. _cavitracer_single:

I. Detection of channels across heterogeneous structures
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

   @> 17 PDBs were parsed in 14.24s.
   1W0E
   @> Checking AtomGroup 1W0E: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 1W0E (len=450) and Chain A from 1TQN (len=468):
   @> 	Match: 447 residues match with 99% sequence identity and 96% overlap.
   4I3Q
   @> Checking AtomGroup 4I3Q: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 4I3Q (len=465) and Chain A from 1TQN (len=468):
   @> 	Match: 465 residues match with 100% sequence identity and 99% overlap.
   5A1P
   @> Checking AtomGroup 5A1P: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 5A1P (len=465) and Chain A from 1TQN (len=468):
   @> 	Match: 463 residues match with 100% sequence identity and 99% overlap.
   5VCC
   @> Checking AtomGroup 5VCC: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 5VCC (len=457) and Chain A from 1TQN (len=468):
   @> 	Match: 457 residues match with 100% sequence identity and 98% overlap.
   6BD6
   @> Checking AtomGroup 6BD6: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6BD6 (len=457) and Chain A from 1TQN (len=468):
   @> 	Match: 457 residues match with 100% sequence identity and 98% overlap.
   6BD8
   @> Checking AtomGroup 6BD8: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6BD8 (len=457) and Chain A from 1TQN (len=468):
   @> 	Match: 457 residues match with 100% sequence identity and 98% overlap.
   6BDI
   @> Checking AtomGroup 6BDI: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6BDI (len=462) and Chain A from 1TQN (len=468):
   @> 	Match: 461 residues match with 100% sequence identity and 99% overlap.
   6BDM
   @> Checking AtomGroup 6BDM: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6BDM (len=456) and Chain A from 1TQN (len=468):
   @> 	Match: 456 residues match with 100% sequence identity and 97% overlap.
   6DA8
   @> Checking AtomGroup 6DA8: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6DA8 (len=451) and Chain A from 1TQN (len=468):
   @> 	Match: 450 residues match with 100% sequence identity and 96% overlap.
   6DAJ
   @> Checking AtomGroup 6DAJ: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6DAJ (len=451) and Chain A from 1TQN (len=468):
   @> 	Match: 450 residues match with 100% sequence identity and 96% overlap.
   6DAL
   @> Checking AtomGroup 6DAL: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6DAL (len=454) and Chain A from 1TQN (len=468):
   @> 	Match: 453 residues match with 100% sequence identity and 97% overlap.
   6MA6
   @> Checking AtomGroup 6MA6: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6MA6 (len=466) and Chain A from 1TQN (len=468):
   @> 	Match: 464 residues match with 100% sequence identity and 99% overlap.
   6MA7
   @> Checking AtomGroup 6MA7: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6MA7 (len=463) and Chain A from 1TQN (len=468):
   @> 	Match: 461 residues match with 99% sequence identity and 99% overlap.
   6OOA
   @> Checking AtomGroup 6OOA: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6OOA (len=452) and Chain A from 1TQN (len=468):
   @> 	Match: 452 residues match with 100% sequence identity and 97% overlap.
   6UNE
   @> Checking AtomGroup 6UNE: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6UNE (len=452) and Chain A from 1TQN (len=468):
   @> 	Match: 452 residues match with 100% sequence identity and 97% overlap.
   6UNG
   @> Checking AtomGroup 6UNG: 1 chains are identified
   @> Checking AtomGroup 1TQN: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from 6UNG (len=453) and Chain A from 1TQN (len=468):
   @> 	Match: 453 residues match with 100% sequence identity and 97% overlap.


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

    @> 3727 atoms and 1 coordinate set(s) were parsed in 0.04s.
    @> The atoms supplied to calcChannels contain protein atoms only.
    @> Substituted 3630 atoms with 22830 homogeneous balls of radius 1.52 Å in 0.17s.
    @> Delaunay tessellation of 22830 points constructed in 0.78s.
    @> Surface and inner simplices filtered in 1.93s.
    @> Cavities: 145 found, 9 deeper than min_depth=5.0 Å and searched for channels, in 0.31s.
    @> Chambers (probe 1.40 Å): 6 of the 9 searched cavities have them; the other 4 are searched whole.
    @>     cavity 0: 26 chambers, 5 of them seeded.
    @>     cavity 1: 3 chambers, 1 of them seeded.
    @>     cavity 2: 1 chamber, seeded.
    @>     cavity 3: 3 chambers, 1 of them seeded.
    @>     cavity 4: 2 chambers, all seeded.
    @>     cavity 7: 1 chamber, none of them deep and large enough to seed; searched whole.
    @> 14 search sites (sp) in 0.06s: one per seeded chamber, one per cavity searched whole.
    @> Channel search (Dijkstra) over 14 search sites in 9 cavities completed in 0.22s.
    @> Found 16 channels and 2 links (a link joins a deep chamber to a shallower one and never reaches the surface).
    @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
    @>     site  start_point [Å]              void                   volume [Å³]  depth [Å]  channels  links
    @>     sp0   [-22.739, -17.762, -10.690]  cavity 0, chamber 1/5         4060        8.2         6      -
    @>     sp1   [-27.340, -26.783, -1.430]   cavity 0, chamber 2/5          262        8.0         1      -
    @>     sp2   [-13.431, -37.424, -4.075]   cavity 5, whole                221        5.5         1      -
    @>     sp3   [-25.317, -31.316, -6.063]   cavity 0, chamber 3/5          206       16.1         -      -  sealed
    @>     sp4   [-24.901, -41.642, -11.119]  cavity 2, chamber 1/1          189        8.9         -      -  sealed
    @>     sp5   [-28.738, -46.393, -21.582]  cavity 6, whole                184        5.1         1      -
    @>     sp6   [-17.407, 3.345, -18.881]    cavity 7, whole                173        5.3         1      -
    @>     sp7   [-11.248, -24.855, 4.687]    cavity 8, whole                159        5.4         1      -
    @>     sp8   [-22.709, -14.923, -23.709]  cavity 0, chamber 4/5          139       14.5         -      1  -> sp0
    @>     sp9   [-19.542, -43.656, -16.249]  cavity 4, chamber 1/2           85       11.5         -      1  -> sp10
    @>     sp10  [-15.393, -42.639, -12.731]  cavity 4, chamber 2/2           80        5.5         1      -
    @>     sp11  [-9.092, -31.126, -5.657]    cavity 1, chamber 1/1           67        5.4         2      -
    @>     sp12  [-27.968, -18.011, -23.364]  cavity 0, chamber 5/5           59        5.3         2      -
    @>     sp13  [-27.143, -36.831, -31.674]  cavity 3, chamber 1/1           54        9.1         -      -  sealed
    @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
    @> The 3 sites marked sealed above report neither a channel nor a link: no route out of them survived - either narrower than bottleneck=1.20 Å, or dropped as a duplicate of a shallower site's, or the void is its own mouth and has nowhere to path to. Lower bottleneck to see how the first kind connect.
    @> Saving 16 channels and 2 links to directory ., one file per object named sp<site>_chl<n> and sp<site>_lnk<n>.
    @> Channel calculation completed in 3.43s.
    @> Channel ID:      Volume [Å³]     Length [Å]      Bottleneck [Å]
    @> channel 0:       775.73          8.2             3.22
    @> channel 1:       967.38          12.29           2.8
    @> channel 2:       97.89           5.69            1.87
    @> channel 3:       637.99          11.02           1.6
    @> channel 4:       591.75          10.78           1.68
    @> channel 5:       72.55           5.45            1.42
    @> channel 6:       1011.37                 16.64           1.65
    @> channel 7:       54.3            5.33            1.23
    @> channel 8:       1032.07                 18.43           1.94
    @> channel 9:       134.88          8.02            1.45
    @> channel 10:      52.55           5.2             1.25
    @> channel 11:      50.39           5.51            1.26
    @> channel 12:      57.23           5.5             1.24
    @> channel 13:      57.52           5.4             1.21
    @> channel 14:      46.83           5.5             1.2
    @> channel 15:      136.16          10.37           1.28
    @> Channel residues were saved to: align__6OOA_Residues_All_channels.txt
    ..
    ..
    @> 3882 atoms and 1 coordinate set(s) were parsed in 0.04s.
    @> The atoms supplied to calcChannels contain protein atoms only.
    @> Substituted 3720 atoms with 23368 homogeneous balls of radius 1.52 Å in 0.16s.
    @> Delaunay tessellation of 23368 points constructed in 0.73s.
    @> Surface and inner simplices filtered in 1.85s.
    @> Cavities: 145 found, 10 deeper than min_depth=5.0 Å and searched for channels, in 0.29s.
    @> Chambers (probe 1.40 Å): 8 of the 10 searched cavities have them; the other 6 are searched whole.
    @>     cavity 0: 14 chambers, 3 of them seeded.
    @>     cavity 2: 4 chambers, 1 of them seeded.
    @>     cavity 3: 1 chamber, seeded.
    @>     cavity 4: 2 chambers, 1 of them seeded.
    @>     cavity 5: 1 chamber, none of them deep and large enough to seed; searched whole.
    @>     cavity 6: 2 chambers, none of them deep and large enough to seed; searched whole.
    @>     cavity 7: 1 chamber, none of them deep and large enough to seed; searched whole.
    @>     cavity 8: 1 chamber, none of them deep and large enough to seed; searched whole.
    @> 12 search sites (sp) in 0.06s: one per seeded chamber, one per cavity searched whole.
    @> Channel search (Dijkstra) over 12 search sites in 10 cavities completed in 0.22s.
    @> Found 13 channels and 1 link (a link joins a deep chamber to a shallower one and never reaches the surface).
    @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
    @>     site  start_point [Å]              void                   volume [Å³]  depth [Å]  channels  links
    @>     sp0   [-20.317, -17.687, -11.117]  cavity 0, chamber 1/3         4379       13.6         3      -
    @>     sp1   [-10.376, -21.751, -13.490]  cavity 1, whole               1083        5.2         -      -  sealed
    @>     sp2   [-28.046, -16.350, -21.344]  cavity 5, whole                328        5.6         1      -
    @>     sp3   [-30.637, -30.877, -28.063]  cavity 6, whole                212        5.0         1      -
    @>     sp4   [-25.189, -42.244, -11.167]  cavity 3, chamber 1/1          206        8.0         1      -
    @>     sp5   [-16.336, -41.507, -13.054]  cavity 7, whole                169        5.3         1      -
    @>     sp6   [-28.885, -46.170, -21.266]  cavity 8, whole                153        5.1         1      -
    @>     sp7   [-9.029, -11.855, -25.151]   cavity 9, whole                128        7.3         -      -  sealed
    @>     sp8   [-9.562, -30.917, -5.403]    cavity 4, chamber 1/1          106        5.4         1      -
    @>     sp9   [-12.520, -24.020, 7.203]    cavity 0, chamber 2/3           59        5.6         2      -
    @>     sp10  [-24.705, -25.684, -11.042]  cavity 0, chamber 3/3           57       23.6         -      1  -> sp0
    @>     sp11  [-28.256, -30.724, -21.745]  cavity 2, chamber 1/1           54        5.3         2      -
    @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
    @> The 2 sites marked sealed above report neither a channel nor a link: no route out of them survived - either narrower than bottleneck=1.20 Å, or dropped as a duplicate of a shallower site's, or the void is its own mouth and has nowhere to path to. Lower bottleneck to see how the first kind connect.
    @> Saving 13 channels and 1 links to directory ., one file per object named sp<site>_chl<n> and sp<site>_lnk<n>.
    @> Channel calculation completed in 3.27s.
    @> Channel ID:      Volume [Å³]     Length [Å]      Bottleneck [Å]
    @> channel 0:       128.0           5.35            1.72
    @> channel 1:       91.93           5.61            1.4
    @> channel 2:       66.55           5.4             1.46
    @> channel 3:       58.31           5.11            1.48
    @> channel 4:       76.47           5.48            1.24
    @> channel 5:       85.11           6.45            1.52
    @> channel 6:       733.63          20.27           2.43
    @> channel 7:       54.88           5.16            1.22
    @> channel 8:       51.96           5.6             1.32
    @> channel 9:       456.21          13.71           1.64
    @> channel 10:      119.14          8.21            1.27
    @> channel 11:      87.29           8.26            1.36
    @> channel 12:      593.35          22.09           1.21
    @> Channel residues were saved to: align__6MA7_Residues_All_channels.txt
    @> 4124 atoms and 1 coordinate set(s) were parsed in 0.05s.
    @> The atoms supplied to calcChannels contain protein atoms only.
    @> Substituted 3675 atoms with 23131 homogeneous balls of radius 1.52 Å in 0.16s.
    @> Delaunay tessellation of 23131 points constructed in 0.71s.
    @> Surface and inner simplices filtered in 1.70s.
    @> Cavities: 150 found, 9 deeper than min_depth=5.0 Å and searched for channels, in 0.29s.
    @> Chambers (probe 1.40 Å): 8 of the 9 searched cavities have them; the other 6 are searched whole.
    @>     cavity 0: 11 chambers, 3 of them seeded.
    @>     cavity 1: 1 chamber, seeded.
    @>     cavity 2: 1 chamber, seeded.
    @>     cavity 3: 1 chamber, none of them deep and large enough to seed; searched whole.
    @>     cavity 4: 2 chambers, none of them deep and large enough to seed; searched whole.
    @>     cavity 5: 2 chambers, none of them deep and large enough to seed; searched whole.
    @>     cavity 6: 3 chambers, none of them deep and large enough to seed; searched whole.
    @>     cavity 7: 1 chamber, none of them deep and large enough to seed; searched whole.
    @> 11 search sites (sp) in 0.06s: one per seeded chamber, one per cavity searched whole.
    @> Channel search (Dijkstra) over 11 search sites in 9 cavities completed in 0.20s.
    @> Found 12 channels and 1 link (a link joins a deep chamber to a shallower one and never reaches the surface).
    @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
    @>     site  start_point [Å]              void                   volume [Å³]  depth [Å]  channels  links
    @>     sp0   [-17.602, -19.728, -11.814]  cavity 0, chamber 1/3         4491       15.1         3      -
    @>     sp1   [-27.975, -17.395, -22.773]  cavity 3, whole                439        5.1         1      -
    @>     sp2   [-15.364, -42.346, -12.635]  cavity 4, whole                314        5.1         1      -
    @>     sp3   [-23.632, -30.164, -27.285]  cavity 5, whole                284        6.8         1      -
    @>     sp4   [-30.186, -31.099, -28.502]  cavity 6, whole                224        5.0         1      -
    @>     sp5   [-25.527, -42.082, -11.325]  cavity 2, chamber 1/1          224        9.3         -      -  sealed
    @>     sp6   [-9.478, -31.244, -5.484]    cavity 1, chamber 1/1          174        8.4         2      -
    @>     sp7   [-17.420, 3.313, -18.466]    cavity 7, whole                159        5.3         1      -
    @>     sp8   [-30.850, -19.290, 3.062]    cavity 8, whole                116        5.1         1      -
    @>     sp9   [-29.245, -32.774, -21.098]  cavity 0, chamber 2/3           75        5.4         1      -
    @>     sp10  [-26.068, -29.878, -22.304]  cavity 0, chamber 3/3           66        6.3         -      1  -> sp9
    @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
    @> The 1 site marked sealed above report neither a channel nor a link: no route out of them survived - either narrower than bottleneck=1.20 Å, or dropped as a duplicate of a shallower site's, or the void is its own mouth and has nowhere to path to. Lower bottleneck to see how the first kind connect.
    @> Saving 12 channels and 1 links to directory ., one file per object named sp<site>_chl<n> and sp<site>_lnk<n>.
    @> Channel calculation completed in 3.07s.
    @> Channel ID:      Volume [Å³]     Length [Å]      Bottleneck [Å]
    @> channel 0:       132.35          5.49            1.84
    @> channel 1:       98.53           5.09            1.83
    @> channel 2:       69.55           5.08            1.63
    @> channel 3:       66.51           5.5             1.46
    @> channel 4:       53.65           5.4             1.29
    @> channel 5:       747.27          22.16           2.48
    @> channel 6:       432.11          15.22           1.63
    @> channel 7:       39.53           5.1             1.2
    @> channel 8:       121.43          6.8             1.28
    @> channel 9:       137.39          9.32            1.5
    @> channel 10:      96.81           8.68            1.22
    @> channel 11:      453.96          18.81           1.3
    @> Channel residues were saved to: align__5VCC_Residues_All_channels.txt


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

    ['align__6OOA_sp10_chl12.pqr',
     'align__1W0E_sp1_chl2.pqr',
     'align__4I3Q_sp1_chl5.pqr',
     'align__4I3Q_sp3_chl2.pqr',
     'align__6MA6_sp1_chl1.pqr',
     'align__6BD8_sp5_chl9.pqr',
     'align__6BD6_sp1_chl1.pqr',
     'align__1W0E_sp13_chl8.pqr',
     'align__6MA7_sp9_chl8.pqr',
     'align__6DAL_sp2_chl4.pqr',
     'align__6DAL_sp5_chl11.pqr',
     'align__5VCC_sp7_chl4.pqr',
     'align__6UNE_sp6_chl11.pqr',
     'align__6UNG_sp12_chl11.pqr',
     'align__1W0E_sp11_chl7.pqr',
     'align__6UNG_sp0_chl1.pqr',
     'align__6BDM_sp10_chl11.pqr',
     'align__5A1P_sp0_chl8.pqr',
     'align__1W0E_sp2_chl11.pqr',
     'align__6UNE_sp10_chl2.pqr',
     'align__6BDM_sp11_chl7.pqr',
     'align__6BD8_sp7_chl7.pqr',
     'align__6UNG_sp10_chl5.pqr',
     'align__6MA6_sp8_chl2.pqr',
     'align__6MA7_sp2_chl1.pqr',
     'align__6MA7_sp0_chl12.pqr',
     'align__6BD8_sp0_chl2.pqr',
     'align__6BD6_sp10_chl0.pqr',
     'align__6BDM_sp0_chl8.pqr',
     'align__6DAJ_sp6_chl16.pqr',
     'align__6MA6_sp5_chl5.pqr',
     'align__6UNG_sp0_chl4.pqr',
     'align__5VCC_sp0_chl5.pqr',
     'align__6UNG_sp8_chl7.pqr',
     'align__6BD6_sp13_chl2.pqr',
     'align__4I3Q_sp0_chl6.pqr',
     'align__5A1P_sp0_chl10.pqr',
     ..
     ..
     'align__1W0E_sp7_chl0.pqr',
     'align__6BD6_sp0_chl6.pqr',
     'align__6DAL_sp0_chl0.pqr',
     'align__5VCC_sp0_chl6.pqr',
     'align__6UNG_sp6_chl6.pqr',
     'align__6MA6_sp0_chl4.pqr',
     'align__6MA7_sp8_chl4.pqr',
     'align__1W0E_sp0_chl3.pqr',
     'align__6BD6_sp4_chl7.pqr']


.. ipython:: python
   :verbatim:

   atoms = parsePDB(pdb_files_new_channels[0].split('_sp')[0]+'.pdb')
   selectChannelBySelection(atoms, pqr_files=pdb_files_new_channels, 
					residue_sele='resid 442')

.. parsed-literal::

    @> 3727 atoms and 1 coordinate set(s) were parsed in 0.05s.
    @> 55 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 35 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 40 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 40 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 60 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 50 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 35 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 55 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 60 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 60 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 50 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 45 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 80 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 90 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 75 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 55 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 115 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 100 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 100 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 50 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 70 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 60 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 60 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 85 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 45 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 190 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 90 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 50 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 125 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 60 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 45 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 95 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 135 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 40 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 75 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 115 atoms and 1 coordinate sets were parsed in 0.00s.
    ..
    ..
    @> 150 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 100 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 55 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 90 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 60 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 45 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 40 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 40 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 40 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 60 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 125 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 55 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 100 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 40 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 95 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 55 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 80 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 50 atoms and 1 coordinate sets were parsed in 0.00s.
    @> Filtered files are now in: selected_files
    @> Selected files: 
    @> align__5A1P_sp0_chl8.pqr align__5A1P_sp0_chl10.pqr align__4I3Q_sp7_chl10.pqr align__5A1P_sp0_chl5.pqr align__6BD8_sp3_chl10.pqr align__5A1P_sp5_chl7.pqr align__5VCC_sp6_chl9.pqr


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

   atoms = parsePDB(pdb_files_new_channels[0].split('_sp')[0]+'.pdb')
   selectChannelBySelection(atoms, pqr_files=pdb_files_new_channels, residue_sele='resid 442',
   folder_name='res442', residues_file=True, param_file=True)

.. parsed-literal::

    @> 3727 atoms and 1 coordinate set(s) were parsed in 0.05s.
    @> 55 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 35 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 40 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 40 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 60 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 50 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 35 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 55 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 60 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 60 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 50 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 45 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 80 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 90 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 75 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 55 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 115 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 100 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 100 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 50 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 70 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 60 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 60 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 85 atoms and 1 coordinate sets were parsed in 0.00s.
    ..
    ..
    @> 125 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 55 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 100 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 40 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 95 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 55 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 80 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 50 atoms and 1 coordinate sets were parsed in 0.00s.
    @> Filtered files are now in: res442
    @> 7 residue row(s) saved to: Selected_channel_residues.txt
    @> 7 parameter row(s) saved to: Selected_channel_parameters.txt
    @> Selected files: 
    @> align__5A1P_sp0_chl8.pqr align__5A1P_sp0_chl10.pqr align__4I3Q_sp7_chl10.pqr align__5A1P_sp0_chl5.pqr align__6BD8_sp3_chl10.pqr align__5A1P_sp5_chl7.pqr align__5VCC_sp6_chl9.pqr


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

   calcChannelSurfaceOverlaps(pqr_files=pdb_files_new_channels_ALL, 
				output_file_name='overlapping_surf.pdb')

.. parsed-literal::

    @> Number of PQR files: 16
    @> Resolution: 0.5
    @> max_proc: 2
    @> Calculating overlaps using 2 processes.
    @> 1055 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1155 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1115 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1125 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 940 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1565 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 835 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 810 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 1015 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1035 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 875 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 810 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1050 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 960 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1020 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1030 atoms and 1 coordinate sets were parsed in 0.01s.
    @> Overlap written to: overlapping_surf.pdb
    @> Number of occupied overlap voxels: 63057

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