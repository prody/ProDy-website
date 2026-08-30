.. _cavitracer_single:

I. Detection of channels in multi-model PDBs
===============================================================================

In this example, we will use the NMR structure of proteorhodopsin with the
PDB code ``2L6X``, which contains 20 models in the PDB file to show how to
predict channels in multi-model PDBs. Multi-model PDB files can also contain
frames from molecular dynamics simulations (MD). Therefore, this analysis
can also be used to analyze MD trajectory in case we use some other the MD
format than ``DCD`` that is analyzed by ProDy tools. 

.. ipython:: python
   :verbatim:

   p3 = parsePDB('2L6X')

.. parsed-literal::

   @> Connecting wwPDB FTP server RCSB PDB (USA).
   @> Downloading PDB files via FTP failed, trying HTTP.
   @> 2l6x downloaded (2l6x.pdb.gz)
   @> PDB download via HTTP completed (1 downloaded, 0 failed).
   @> 3669 atoms and 20 coordinate set(s) were parsed in 0.23s.

To detect channels in multi-model PDB files or in MD trajectories, we need
to use :func:`.scalcChannelsMultipleFrames`. To speed up the calculations,
``max_proc`` can be used. Below, we are using four processors.

.. ipython:: python
   :verbatim:

   channels3, surfaces3 = calcChannelsMultipleFrames(p3, max_proc=4)

.. parsed-literal::

   @> Frame/model: 0
   @> Frame/model: 2
   @> Frame/model: 4
   @> Frame/model: 6
   @> WARNING The atoms supplied to calcChannels() contain non-protein components: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein structure, provide an appropriate selection, for example atoms.select('protein').
   @> WARNING The atoms supplied to calcChannels() contain non-protein components: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein structure, provide an appropriate selection, for example atoms.select('protein').
   @> WARNING The atoms supplied to calcChannels() contain non-protein components: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein structure, provide an appropriate selection, for example atoms.select('protein').
   @> WARNING The atoms supplied to calcChannels() contain non-protein components: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein structure, provide an appropriate selection, for example atoms.select('protein').
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 Å in 0.19s.
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 Å in 0.19s.
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 Å in 0.19s.
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 Å in 0.19s.
   @> Delaunay tessellation of 48645 points constructed in 2.38s.
   @> Delaunay tessellation of 48645 points constructed in 2.41s.
   @> Delaunay tessellation of 48645 points constructed in 2.40s.
   @> Delaunay tessellation of 48645 points constructed in 2.47s.
   @> Surface and inner simplices filtered in 3.01s.
   @> Surface and inner simplices filtered in 3.12s.
   @> Surface and inner simplices filtered in 3.07s.
   @> Surface and inner simplices filtered in 3.17s.
   @> Surface cavities: 246 found, 2 deeper than min_depth=5.0 Å and searched for channels, in 0.42s.
   @> Surface cavities: 244 found, 3 deeper than min_depth=5.0 Å and searched for channels, in 0.38s.
   @> Surface cavities: 252 found, 4 deeper than min_depth=5.0 Å and searched for channels, in 0.39s.
   @> Surface cavities: 265 found, 2 deeper than min_depth=5.0 Å and searched for channels, in 0.43s.
   @> Chambers (probe 1.40 Å): 1 of the 2 searched cavities have them; the other 1 is searched whole.
   @>     cavity 0: 28 chambers, 10 of them seeded.
   @> 11 search sites (sp) in 0.04s: one per seeded chamber, one per cavity searched whole.
   @> Chambers (probe 1.40 Å): 2 of the 3 searched cavities have them; the other 1 is searched whole.
   @>     cavity 0: 11 chambers, 7 of them seeded.
   @>     cavity 1: 5 chambers, 4 of them seeded.
   @> 12 search sites (sp) in 0.06s: one per seeded chamber, one per cavity searched whole.
   @> Chambers (probe 1.40 Å): 4 of the 4 searched cavities have them.
   @>     cavity 0: 20 chambers, 9 of them seeded.
   @>     cavity 1: 2 chambers, all seeded.
   @>     cavity 2: 1 chamber, seeded.
   @>     cavity 3: 1 chamber, seeded.
   @> 13 search sites (sp) in 0.05s: one per seeded chamber, one per cavity searched whole.
   @> Chambers (probe 1.40 Å): 1 of the 2 searched cavities have them; the other 1 is searched whole.
   @>     cavity 0: 32 chambers, 9 of them seeded.
   @> 10 search sites (sp) in 0.05s: one per seeded chamber, one per cavity searched whole.
   @> Channel search (Dijkstra) over 12 search sites in 3 cavities completed in 1.18s.
   @> Found 35 channels and 9 links (a link joins a deep chamber to a shallower one and never reaches the surface).
   @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
   @>     site  void                   volume [Å³]  depth [Å]  channels  links
   @>     sp0   cavity 0, chamber 1/7          522       10.2         6      1  -> sp3
   @>     sp1   cavity 0, chamber 2/7          482        9.5         6      1  -> sp0
   @>     sp2   cavity 1, chamber 1/4          184        5.3         5      -
   @>     sp3   cavity 0, chamber 3/7          182        7.6         3      -
   @>     sp4   cavity 2, whole                171        5.7         -      -  sealed
   @>     sp5   cavity 0, chamber 4/7          135        9.3         2      2  -> sp0, sp1
   @>     sp6   cavity 0, chamber 5/7          118       20.0         -      -  sealed
   @>     sp7   cavity 1, chamber 2/4          102        8.5         2      1  -> sp8
   @>     sp8   cavity 1, chamber 3/4           86        5.3         4      -
   @>     sp9   cavity 1, chamber 4/4           77        6.3         1      1  -> sp8
   @>     sp10  cavity 0, chamber 6/7           64        5.4         6      1  -> sp0
   @>     sp11  cavity 0, chamber 7/7           33       18.1         -      2  -> sp1, sp5
   @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
   @> The 1 site marked sealed above report neither a channel nor a link: every route out of them is narrower than bottleneck=0.90 Å. Lower it to see how they connect.
   @> No output path given.
   @> Channel calculation completed in 7.61s.
   ..
   ..
   @> Frame/model: 17
   @> WARNING The atoms supplied to calcChannels() contain non-protein components: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein structure, provide an appropriate selection, for example atoms.select('protein').
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 Å in 0.16s.
   @> Channel search (Dijkstra) over 15 search sites in 3 cavities completed in 1.94s.
   @> Found 58 channels and 11 links (a link joins a deep chamber to a shallower one and never reaches the surface).
   @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
   @>     site  void                     volume [Å³]  depth [Å]  channels  links
   @>     sp0   cavity 0, chamber 1/13           751       11.7         3      2  -> sp7, sp1
   @>     sp1   cavity 0, chamber 2/13           399        5.7         8      -
   @>     sp2   cavity 0, chamber 3/13           133        5.1         7      -
   @>     sp3   cavity 0, chamber 4/13           126        5.3        10      -
   @>     sp4   cavity 1, whole                  109        5.0         2      -
   @>     sp5   cavity 0, chamber 5/13            82       12.5         8      2  -> sp1, sp2
   @>     sp6   cavity 0, chamber 6/13            69        5.0         6      -
   @>     sp7   cavity 0, chamber 7/13            64        5.3         1      -
   @>     sp8   cavity 0, chamber 8/13            63       22.7         1      2  -> sp5, sp0
   @>     sp9   cavity 2, whole                   53        5.1         -      -  sealed
   @>     sp10  cavity 0, chamber 9/13            48        7.2         5      -
   @>     sp11  cavity 0, chamber 10/13           46        7.3         3      1  -> sp12
   @>     sp12  cavity 0, chamber 11/13           43        6.3         1      1  -> sp6
   @>     sp13  cavity 0, chamber 12/13           36        7.1         2      2  -> sp0, sp1
   @>     sp14  cavity 0, chamber 13/13           30        6.1         1      1  -> sp3
   @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
   @> No output path given.
   @> Channel calculation completed in 7.22s.
   @> Delaunay tessellation of 48645 points constructed in 1.91s.
   @> Surface and inner simplices filtered in 1.97s.
   @> Surface cavities: 230 found, 2 deeper than min_depth=5.0 Å and searched for channels, in 0.33s.
   @> Chambers (probe 1.40 Å): 2 of the 2 searched cavities have them.
   @>     cavity 0: 20 chambers, 6 of them seeded.
   @>     cavity 1: 4 chambers, 3 of them seeded.
   @> 9 search sites (sp) in 0.05s: one per seeded chamber, one per cavity searched whole.
   @> Channel search (Dijkstra) over 9 search sites in 2 cavities completed in 2.65s.
   @> Found 41 channels and 6 links (a link joins a deep chamber to a shallower one and never reaches the surface).
   @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
   @>     site  void                   volume [Å³]  depth [Å]  channels  links
   @>     sp0   cavity 0, chamber 1/6         1186       11.8        19      2  -> sp3, sp2
   @>     sp1   cavity 0, chamber 2/6          189        9.4         -      1  -> sp2
   @>     sp2   cavity 0, chamber 3/6          148        5.0        11      -
   @>     sp3   cavity 0, chamber 4/6          107        5.1         4      -
   @>     sp4   cavity 1, chamber 1/3           90        8.7         2      -
   @>     sp5   cavity 0, chamber 5/6           89        5.1         3      -
   @>     sp6   cavity 1, chamber 2/3           69        5.2         2      1  -> sp4
   @>     sp7   cavity 0, chamber 6/6           49       10.8         -      1  -> sp0
   @>     sp8   cavity 1, chamber 3/3           41       13.4         -      1  -> sp4
   @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
   @> No output path given.
   @> Channel calculation completed in 7.40s.


Channels are stored in a list:

.. ipython:: python
   :verbatim:

   channels3

.. parsed-literal::

    [[<prody.proteins.channels.Channel at 0x798986bc3100>,
    <prody.proteins.channels.Channel at 0x798986bc1ed0>,
    <prody.proteins.channels.Channel at 0x798986bc0250>,
    <prody.proteins.channels.Channel at 0x798986bc04f0>,
    <prody.proteins.channels.Channel at 0x798986bc01c0>,
    <prody.proteins.channels.Channel at 0x798987b3de10>,
    <prody.proteins.channels.Channel at 0x798987b3d750>,
    <prody.proteins.channels.Channel at 0x798987b3cb50>,
   [<prody.proteins.channels.Channel at 0x798a05e22c20>,
   ..]]

To have acess to a particular frame/model, we should treat it as a list of
elements, where elements are predicted channels. To display channels for
model #0, use :func:`.showChannels`:

.. ipython:: python
   :verbatim:

   showChannels(channels3[0])

.. figure:: images/cavitracer_figure11.jpg
   :scale: 50 %


We can also visualize one particular channel from model #2:

.. ipython:: python
   :verbatim:

   showChannels(channels3[2][1])

.. figure:: images/cavitracer_figure12.jpg
   :scale: 50 %


Visualization with protein required building a 3D model of the protein as a
TriangleMesh using :func:`.getVmdModel`. Below we will generate two models.
One for model #0 and second for model #2.

.. ipython:: python
   :verbatim:

   p3.setACSIndex(0)
   p3

.. parsed-literal::

   <AtomGroup: 2L6X (3669 atoms; active #0 of 20 coordsets)>

.. ipython:: python
   :verbatim:

   vmd_path = '/usr/local/bin/vmd'
   model3_0 = getVmdModel(vmd_path, p3)

.. parsed-literal::

   @> Model created successfully.

.. ipython:: python
   :verbatim:

   p3.setACSIndex(2)
   p3

.. parsed-literal::

   <AtomGroup: 2L6X (3669 atoms; active #2 of 20 coordsets)>

.. ipython:: python
   :verbatim:

   model3_2 = getVmdModel(vmd_path, p3)

.. parsed-literal::

   @> Model created successfully.

We generated two models, for model #0 and model #2, which contains 28210
points and 56400 triangles.

.. ipython:: python
   :verbatim:

   model3_0

.. parsed-literal::

   TriangleMesh with 28210 points and 56400 triangles.


.. ipython:: python
   :verbatim:

   showChannels(channels3[0], model=model3_0)

.. figure:: images/cavitracer_figure13.jpg
   :scale: 50 %


.. ipython:: python
   :verbatim:

   model3_2

.. parsed-literal::

   TriangleMesh with 28210 points and 56400 triangles.

.. ipython:: python
   :verbatim:

   showChannels(channels3[2][1], model=model3_2)

.. figure:: images/cavitracer_figure14.jpg
   :scale: 50 %


Access to the parameters of the channels is provided by
:func:`.getChannelParameters`:

.. ipython:: python
   :verbatim:

   getChannelParameters(channels3)

.. parsed-literal::

   @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame 0
   @> channel 0: 	87.62 		2.6 		2.06
   @> channel 1: 	107.63 		3.88 		2.16
   @> channel 2: 	148.44 		5.41 		2.16
   @> channel 3: 	173.05 		6.39 		2.13
   @> channel 4: 	187.64 		8.33 		1.97
   @> channel 5: 	93.46 		6.13 		1.01
   @> channel 6: 	79.86 		6.21 		1.32
   @> channel 7: 	139.58 		9.03 		1.83
   @> channel 8: 	149.08 		9.73 		1.83
   @> channel 9: 	106.1 		8.05 		1.55
   @> channel 10: 	121.32 		8.35 		1.55
   @> channel 11: 	97.23 		7.52 		1.46
   @> channel 12: 	107.48 		8.07 		1.47
   @> channel 13: 	133.14 		8.92 		1.37
   @> channel 14: 	102.46 		8.53 		1.38
   @> channel 15: 	90.95 		8.32 		1.16
   @> channel 16: 	98.57 		8.46 		1.27
   @> channel 17: 	113.81 		9.86 		1.39
   @> channel 18: 	186.39 		11.74 		1.39
   @> channel 19: 	115.66 		10.36 		1.33
   @> channel 20: 	73.21 		8.31 		1.06
   @> channel 21: 	113.18 		10.62 		1.22
   @> channel 22: 	106.94 		11.44 		1.25
   @> channel 23: 	88.67 		9.98 		1.13
   @> channel 24: 	92.99 		9.4 		1.14
   @> channel 25: 	49.62 		7.47 		0.97
   @> channel 26: 	116.08 		12.29 		0.91
   @> channel 27: 	125.73 		11.86 		0.94
   @> channel 28: 	98.8 		11.24 		0.94
   @> channel 29: 	187.48 		18.43 		1.34
   @> channel 30: 	143.3 		15.89 		1.14
   @> channel 31: 	140.87 		15.93 		1.04
   @> channel 32: 	127.16 		13.66 		0.94
   @> channel 33: 	122.64 		15.33 		0.97
   @> channel 34: 	116.73 		13.05 		0.95
   @> channel 35: 	130.76 		16.42 		0.96
   @> channel 36: 	172.19 		19.01 		1.14
   @> channel 37: 	150.78 		18.11 		0.97
   @> channel 38: 	165.42 		21.49 		1.08
   @> channel 39: 	167.51 		21.94 		1.13
   @> channel 40: 	162.91 		22.37 		1.08
   @> channel 41: 	178.45 		22.62 		1.13
   @> channel 42: 	203.18 		25.13 		0.95
   @> channel 43: 	195.4 		25.13 		1.09
   @> channel 44: 	195.18 		25.67 		1.05
   @> channel 45: 	214.16 		26.81 		1.09
   @> channel 46: 	175.04 		27.29 		0.94
   @> channel 47: 	174.68 		28.86 		0.94
   @> Frame 1
   @> channel 0: 	111.22 		5.44 		1.38
   @> channel 1: 	103.82 		6.92 		1.4
   @> channel 2: 	93.82 		6.28 		1.6
   @> channel 3: 	71.48 		6.19 		1.54
   @> channel 4: 	78.46 		6.39 		1.55
   @> channel 5: 	62.05 		5.92 		1.09
   @> channel 6: 	150.05 		8.81 		1.57
   @> channel 7: 	86.4 		7.77 		1.35
   @> channel 8: 	110.98 		8.41 		1.31
   @> channel 9: 	86.73 		8.76 		1.31
   @> channel 10: 	102.02 		7.3 		1.06
   @> channel 11: 	101.08 		9.75 		1.35
   @> channel 12: 	113.74 		10.47 		1.35
   ..
   ..
   @> Frame 18
   @> channel 0: 	117.63 		3.72 		2.05
   @> channel 1: 	111.04 		3.88 		2.03
   @> channel 2: 	159.0 		6.11 		2.06
   @> channel 3: 	178.08 		5.71 		2.27
   @> channel 4: 	110.01 		5.03 		1.86
   @> channel 5: 	142.49 		5.8 		1.92
   @> channel 6: 	178.68 		7.55 		1.7
   @> channel 7: 	91.09 		4.45 		1.48
   @> channel 8: 	210.35 		9.46 		1.85
   @> channel 9: 	128.41 		6.36 		1.52
   @> channel 10: 	161.35 		8.57 		1.74
   @> channel 11: 	105.98 		5.78 		1.31
   @> channel 12: 	203.89 		10.05 		1.66
   @> channel 13: 	71.13 		5.42 		1.37
   @> channel 14: 	90.79 		6.31 		1.59
   @> channel 15: 	61.05 		5.48 		1.31
   @> channel 16: 	209.48 		11.27 		1.8
   @> channel 17: 	140.91 		8.68 		1.32
   @> channel 18: 	117.62 		8.65 		1.35
   @> channel 19: 	217.04 		12.03 		1.59
   @> channel 20: 	87.53 		7.83 		1.18
   @> channel 21: 	71.3 		7.71 		1.1
   @> channel 22: 	131.67 		8.34 		1.01
   @> channel 23: 	112.78 		10.64 		1.35
   @> channel 24: 	83.38 		7.52 		1.06
   @> channel 25: 	73.77 		7.73 		0.94
   @> channel 26: 	24.11 		5.41 		0.94
   @> channel 27: 	65.45 		8.63 		1.14
   @> channel 28: 	129.16 		11.32 		1.14
   @> channel 29: 	69.33 		7.99 		1.01
   @> channel 30: 	85.1 		9.22 		1.06
   @> channel 31: 	51.38 		8.08 		1.1
   @> channel 32: 	118.08 		10.97 		1.01
   @> channel 33: 	51.89 		8.3 		0.99
   @> channel 34: 	82.51 		9.78 		0.96
   @> channel 35: 	225.55 		18.01 		1.35
   @> channel 36: 	83.76 		10.71 		1.01
   @> channel 37: 	251.52 		20.82 		1.38
   @> channel 38: 	37.75 		7.39 		0.94
   @> channel 39: 	48.1 		7.99 		0.93
   @> channel 40: 	65.08 		8.82 		0.94
   @> channel 41: 	155.6 		15.24 		0.91
   @> channel 42: 	129.39 		14.32 		0.91
   @> channel 43: 	242.72 		21.66 		1.38
   @> channel 44: 	256.59 		18.63 		1.14
   @> channel 45: 	90.24 		12.92 		0.91
   @> channel 46: 	254.75 		19.15 		0.95
   @> channel 47: 	37.66 		9.28 		0.91
   @> channel 48: 	250.39 		19.87 		0.95
   @> channel 49: 	216.89 		19.7 		0.96
   @> channel 50: 	254.69 		21.35 		0.96
   @> channel 51: 	199.4 		19.7 		0.96
   @> channel 52: 	92.98 		15.61 		0.97
   @> channel 53: 	146.14 		22.84 		0.91
   @> channel 54: 	132.1 		23.77 		0.91
   @> channel 55: 	143.94 		24.91 		0.91
   @> channel 56: 	133.77 		24.52 		0.91
   @> channel 57: 	107.29 		22.58 		0.9

   [([2.5986838910881516,
      3.876062698520582,
      5.4075068782393005,
      6.388078872187722,
      8.325246507785625,
      6.126309629480059,
      6.212556155969711,
      9.031218685795304,
      9.73465552906654,
      8.046363460710767,
      8.351422527140901,
      7.519061930609062,
      8.073020161229532,
      8.917532629970726,
      8.531961913191681,
      8.320946140304821,
      8.457584891836165,
      9.861993651120457,
      11.738731428780198,
      10.356479422519385,
      8.313864852983162,
      10.620597553894724,
      11.443141116821732,
      9.980806388554237,
      9.400418187229269,
      7.472384571818239,
      12.288866270330228,
      11.860717211900424,
      11.236259672259251,
      18.428052918278688,
      15.89225351104158,
      15.929484594257707,
      13.660741541929603,
      15.32834564051029,
      13.049979756695446,
      16.420477819689857,
      19.011331474497506,
      ..
      112.78137149103382,
      83.37867845402317,
      73.77245468867955,
      24.113779399290355,
      65.45170431480145,
      129.15978509016915,
      69.32608280672586,
      85.09598841324143,
      51.38223042596009,
      118.08387538617845,
      51.88960266096877,
      82.51212276104394,
      225.55115627400454,
      83.75513962413046,
      251.51981936034022,
      37.75104103212885,
      48.102225488741425,
      65.08388396250602,
      155.60442581153487,
      129.39332639465985,
      242.72352467964407,
      256.59437533958345,
      90.24296342294183,
      254.7472417569916,
      37.66293882758854,
      250.3912357501308,
      216.88535195444786,
      254.69390491156682,
      199.39525430000518,
      92.9751162058638,
      146.13903931740356,
      132.10008335569137,
      143.93957960720397,
      133.77351470361398,
      107.29493255596638])]


Access to the residues that are forming the channels is provided by
:func:`.getChannelResidueNames` function. Below, the example on how to
obtain information for frame #0.

.. ipython:: python
   :verbatim:

   p3.setACSIndex(0)
   getChannelResidueNames(p3, channels3[0])

.. parsed-literal::

   ['channel0: SER43:A, THR44:A, PHE47:A, THR69:A, GLY70:A, PHE73:A, ILE232:A, GLY235:A, LEU236:A',
    'channel1: PHE137:A, GLY141:A, GLY144:A, ALA147:A, ALA148:A, ALA151:A, TYR200:A, TYR204:A, TYR223:A, RET301:A',
    'channel2: PHE137:A, MET140:A, GLY141:A, GLY144:A, ALA147:A, ALA148:A, ALA151:A, TYR200:A, TYR204:A, LEU219:A, TYR223:A, RET301:A',
    'channel3: SER43:A, THR44:A, PHE47:A, THR69:A, GLY70:A, PHE73:A, LYS231:A, ILE232:A, GLY235:A, LEU236:A',
    'channel4: ARG80:A, TRP83:A, THR91:A, ARG94:A, PHE137:A, MET140:A, GLY141:A, GLY144:A, ALA147:A, ALA148:A, ALA151:A, TYR200:A, TYR204:A, TYR223:A, RET301:A',
    'channel5: SER43:A, PHE47:A, GLY66:A, THR69:A, GLY70:A, LYS231:A, ILE232:A, GLY235:A, LEU236:A',
    'channel6: PHE47:A, LEU62:A, SER65:A, GLY66:A, THR69:A, GLY70:A, PHE234:A, GLY235:A',
    'channel7: TYR110:A, LEU113:A, ALA114:A, ALA116:A, THR117:A, ALA120:A, GLY171:A, ALA174:A, ALA178:A, SER179:A, VAL182:A, TYR186:A',
    'channel8: TYR110:A, LEU113:A, ALA114:A, ALA116:A, THR117:A, ALA120:A, GLY171:A, ALA174:A, CYS175:A, ALA178:A, SER179:A, VAL182:A, TYR186:A',
    'channel9: TYR110:A, LEU113:A, ALA114:A, ALA120:A, GLY171:A, LYS172:A, VAL182:A, GLN183:A, TYR186:A',
    'channel10: TYR110:A, LEU113:A, ALA114:A, ALA120:A, LEU166:A, TRP167:A, ALA168:A, GLY169:A, GLY171:A, LYS172:A, VAL182:A, TYR186:A',
    'channel11: TYR110:A, LEU113:A, ALA114:A, THR117:A, VAL119:A, ALA120:A, GLY121:A, SER122:A, GLU170:A, GLY171:A, ALA174:A, VAL182:A, TYR186:A',
    'channel12: TYR110:A, LEU113:A, ALA114:A, THR117:A, VAL119:A, ALA120:A, GLU170:A, GLY171:A, ALA174:A, VAL182:A, TYR186:A',
    'channel13: VAL133:A, PHE152:A, GLY155:A, CYS156:A, ALA158:A, TRP159:A, TRP197:A, ALA198:A, PRO201:A, RET301:A',
    'channel14: TYR110:A, LEU113:A, ALA114:A, ALA120:A, GLY121:A, SER122:A, LYS125:A, GLY169:A, GLU170:A, GLY171:A, VAL182:A, TYR186:A',
    'channel15: SER43:A, PHE46:A, PHE47:A, THR69:A, ILE232:A, GLY235:A, LEU236:A',
    'channel16: LEU132:A, VAL133:A, VAL136:A, ILE154:A, GLY155:A, ALA158:A, TRP159:A, TRP197:A, RET301:A',
    'channel17: TYR110:A, LEU113:A, ALA114:A, VAL119:A, ALA120:A, GLY121:A, SER122:A, GLU170:A, GLY171:A, VAL182:A, TYR186:A',
    'channel18: TYR76:A, ARG80:A, TRP83:A, ARG94:A, TRP98:A, PHE137:A, GLY141:A, TYR200:A, TYR223:A, ALA226:A, ASP227:A, ASN230:A, RET301:A',
    'channel19: SER43:A, PHE46:A, PHE47:A, GLU50:A, THR69:A, ILE232:A, GLY235:A, LEU236:A, TRP239:A',
    'channel20: PHE47:A, LEU62:A, SER65:A, GLY66:A, THR69:A, GLY70:A, PHE234:A, GLY235:A',
    'channel21: TYR110:A, LEU113:A, ALA114:A, ALA120:A, GLY171:A, LYS172:A, VAL182:A, GLN183:A, TYR186:A',
    'channel22: LEU105:A, ILE106:A, MET189:A, ILE192:A, ILE193:A, PHE195:A, GLY196:A, TRP197:A, ILE199:A, VAL229:A, ASN230:A, LEU233:A, PHE234:A, RET301:A',
    'channel23: VAL133:A, GLY155:A, CYS156:A, ALA158:A, TRP159:A, ILE194:A, TRP197:A, ALA198:A, PRO201:A, RET301:A',
    'channel24: VAL136:A, PHE137:A, MET140:A, MET146:A, ALA147:A, ALA148:A, PRO150:A, ALA151:A, ILE154:A, TYR200:A, TYR204:A, TYR223:A, RET301:A',
    'channel25: ILE106:A, CYS107:A, TYR110:A, LEU111:A, ALA120:A, SER122:A, LEU123:A, LYS126:A',
    'channel26: TYR110:A, LEU113:A, ALA114:A, ALA120:A, GLY171:A, LYS172:A, THR177:A, ALA178:A, VAL182:A, GLN183:A, TYR186:A',
    'channel27: LYS57:A, TRP58:A, SER61:A, LEU62:A, SER65:A, PHE109:A, ILE112:A, LEU113:A, ALA116:A, VAL182:A, ILE238:A',
    'channel28: PHE137:A, ALA147:A, ALA148:A, ALA151:A, PHE152:A, TYR200:A, PRO201:A, TYR204:A, PHE205:A, TYR223:A, RET301:A',
    'channel29: TRP98:A, LEU99:A, VAL102:A, PRO103:A, LEU105:A, ILE106:A, LYS126:A, LEU127:A, GLY130:A, SER131:A, VAL133:A, MET134:A, MET189:A, ILE193:A, TRP197:A, ASN230:A, PHE234:A, RET301:A',
    'channel30: VAL102:A, PRO103:A, LEU105:A, ILE106:A, LEU123:A, LYS126:A, LEU127:A, GLY130:A, MET162:A, MET189:A, ILE193:A, TRP197:A, ASN230:A, PHE234:A, RET301:A',
    'channel31: LEU105:A, ILE106:A, MET189:A, ILE192:A, ILE193:A, GLY196:A, TRP197:A, ILE199:A, LEU225:A, ALA226:A, VAL229:A, ASN230:A, LEU233:A, PHE234:A, RET301:A',
    'channel32: PHE137:A, ALA147:A, ALA148:A, TRP149:A, ALA151:A, PHE152:A, TYR200:A, PRO201:A, TYR204:A, PHE205:A, TYR208:A, TYR223:A, RET301:A',
    'channel33: VAL133:A, ALA148:A, ALA151:A, PHE152:A, GLY155:A, CYS156:A, ALA158:A, TRP159:A, TRP197:A, ALA198:A, PRO201:A, TYR204:A, PHE205:A, RET301:A',
    'channel34: LEU132:A, VAL133:A, VAL136:A, PHE137:A, MET140:A, ALA147:A, ALA148:A, ALA151:A, ILE154:A, ALA158:A, TYR200:A, TYR204:A, TYR223:A, RET301:A',
    'channel35: SER43:A, PHE46:A, PHE47:A, PHE48:A, GLU50:A, ARG51:A, VAL54:A, LEU62:A, THR69:A, ILE232:A, GLY235:A, LEU236:A, TRP239:A',
    'channel36: SER43:A, PHE46:A, PHE47:A, GLU50:A, ARG51:A, VAL54:A, LEU62:A, THR63:A, GLY66:A, LEU67:A, THR69:A, ILE232:A, GLY235:A, LEU236:A, TRP239:A',
    'channel37: VAL133:A, ALA148:A, TRP149:A, ALA151:A, PHE152:A, GLY155:A, CYS156:A, ALA158:A, TRP159:A, TRP197:A, ALA198:A, PRO201:A, TYR204:A, PHE205:A, TYR208:A, RET301:A',
    'channel38: LEU105:A, ILE106:A, LYS126:A, CYS156:A, TRP159:A, MET162:A, ILE163:A, MET189:A, ILE193:A, ILE194:A, TRP197:A, ALA198:A, ASN230:A, PHE234:A, RET301:A',
    'channel39: SER43:A, PHE46:A, PHE47:A, GLU50:A, ARG51:A, VAL54:A, TRP58:A, LYS59:A, LEU62:A, THR63:A, THR69:A, ILE232:A, GLY235:A, LEU236:A, TRP239:A',
    'channel40: LEU105:A, ILE106:A, LYS126:A, CYS156:A, TRP159:A, MET162:A, ILE163:A, MET189:A, ILE193:A, ILE194:A, TRP197:A, ALA198:A, ASN230:A, PHE234:A, RET301:A',
    'channel41: SER43:A, PHE46:A, PHE47:A, GLU50:A, ARG51:A, ASP52:A, VAL54:A, TRP58:A, LYS59:A, LEU62:A, THR63:A, THR69:A, ILE232:A, GLY235:A, LEU236:A, TRP239:A',
    'channel42: LEU105:A, ILE106:A, PHE109:A, TYR110:A, LYS125:A, LYS126:A, VAL129:A, TYR161:A, MET162:A, ILE163:A, GLU165:A, LEU166:A, MET189:A, MET190:A, ILE193:A, TRP197:A, ASN230:A, PHE234:A, RET301:A',
    'channel43: TYR95:A, TRP98:A, LEU99:A, VAL102:A, LEU105:A, ILE106:A, LYS126:A, GLY130:A, VAL133:A, MET134:A, PHE137:A, GLY138:A, MET189:A, ILE193:A, TRP197:A, ASN230:A, PHE234:A, RET301:A',
    'channel44: TYR95:A, TRP98:A, LEU99:A, VAL102:A, LEU105:A, ILE106:A, LYS126:A, GLY130:A, VAL133:A, MET134:A, PHE137:A, GLY138:A, MET189:A, ILE193:A, TRP197:A, ASN230:A, PHE234:A, RET301:A',
    'channel45: THR91:A, TYR95:A, TRP98:A, LEU99:A, VAL102:A, LEU105:A, ILE106:A, LYS126:A, GLY130:A, VAL133:A, MET134:A, PHE137:A, GLY138:A, MET189:A, ILE193:A, TRP197:A, ASN230:A, PHE234:A, RET301:A',
    'channel46: TRP58:A, LEU105:A, ILE106:A, PHE109:A, TYR110:A, LEU113:A, LYS126:A, LEU166:A, VAL182:A, ALA185:A, TYR186:A, MET189:A, MET190:A, ILE193:A, TRP197:A, ASN230:A, PHE234:A, RET301:A',
    'channel47: TRP58:A, LEU105:A, ILE106:A, PHE109:A, TYR110:A, LEU113:A, LYS126:A, LEU166:A, ALA181:A, VAL182:A, ALA185:A, TYR186:A, MET189:A, MET190:A, ILE193:A, TRP197:A, ASN230:A, PHE234:A, RET301:A']


Second frame is accessible by setting :meth:`.setACSIndex` to ``1``.
To analyze the results for the second model, we need to select second set of
data in ``channels3`` prediction. Additionally, we will display residues
using one letter code and save the results to file by using
``residues_file_name``.

.. ipython:: python
   :verbatim:

   p3.setACSIndex(1)
   getChannelResidueNames(p3, channels3[1], one_letter_aa=True, residues_file_name='results_frame1')

.. parsed-literal::

   ['channel0: Y110:A, L113:A, A114:A, A120:A, G121:A, S122:A, G169:A, E170:A, G171:A, A174:A, V182:A, Y186:A',
    'channel1: F152:A, G155:A, C156:A, W159:A, W197:A, A198:A, P201:A, V202:A, RET301:A',
    'channel2: F137:A, M140:A, G141:A, G144:A, M146:A, A147:A, P150:A, I154:A, Y204:A, Y223:A, RET301:A',
    'channel3: F137:A, A147:A, A151:A, Y200:A, Y204:A, Y223:A, RET301:A',
    'channel4: F137:A, M140:A, M146:A, A147:A, P150:A, I154:A, Y200:A, Y204:A, Y223:A, RET301:A',
    'channel5: V136:A, F137:A, M140:A, M146:A, A147:A, P150:A, I154:A, RET301:A',
    'channel6: L113:A, A114:A, T117:A, A120:A, E170:A, G171:A, A174:A, T177:A, A178:A, V182:A, Y186:A',
    'channel7: K57:A, W58:A, S61:A, F109:A, I112:A, L113:A, A115:A, I238:A',
    'channel8: F137:A, M140:A, G141:A, G144:A, M146:A, A147:A, P150:A, I154:A, Y204:A, Y223:A, RET301:A',
    'channel9: F137:A, G141:A, G144:A, A147:A, A151:A, Y200:A, Y204:A, Y223:A, RET301:A',
    'channel10: L113:A, A114:A, T117:A, V119:A, A120:A, G121:A, E170:A, G171:A, A174:A, V182:A, Y186:A',
    'channel11: K57:A, W58:A, S61:A, F109:A, I112:A, L113:A, A116:A, A181:A, V182:A, I238:A',
    'channel12: K57:A, W58:A, S61:A, F109:A, I112:A, L113:A, A116:A, A181:A, V182:A, I238:A',
    'channel13: W98:A, L99:A, V102:A, P103:A, L127:A, G130:A, S131:A, V133:A, M134:A, W197:A, RET301:A',
    'channel14: S43:A, T44:A, F47:A, L62:A, S65:A, G66:A, T69:A, G70:A, F73:A, F109:A, K231:A, I232:A, F234:A, G235:A, I238:A, W239:A',
    'channel15: L113:A, A114:A, A116:A, T117:A, A120:A, E170:A, G171:A, A174:A, A178:A, S179:A, V182:A, Y186:A',
    'channel16: W98:A, L105:A, M189:A, I192:A, I193:A, F195:A, G196:A, I199:A, V229:A, N230:A, K231:A, L233:A, F234:A, RET301:A',
    'channel17: S43:A, F47:A, L62:A, S65:A, G66:A, T69:A, F109:A, K231:A, I232:A, F234:A, G235:A, L236:A, I238:A, W239:A',
    'channel18: W98:A, L105:A, M189:A, I192:A, I193:A, V229:A, N230:A, K231:A, L233:A, F234:A, I237:A, RET301:A',
    'channel19: W98:A, L105:A, M189:A, I192:A, I193:A, F195:A, G196:A, V229:A, N230:A, K231:A, L233:A, F234:A, I237:A, RET301:A',
    'channel20: F33:A, V36:A, F73:A, M77:A, Y78:A, R80:A, G81:A, N224:A, D227:A, F228:A',
    'channel21: Y110:A, L113:A, A120:A, G121:A, S122:A, K126:A, M162:A, L166:A, G169:A, E170:A, G171:A, Y186:A, M189:A, M190:A, I193:A',
    'channel22: Y76:A, T91:A, V92:A, R94:A, Y95:A, W98:A, G130:A, V133:A, M134:A, F137:A, G138:A, E142:A, W197:A, RET301:A',
    'channel23: S43:A, F46:A, F47:A, E50:A, L62:A, S65:A, G66:A, T69:A, F109:A, K231:A, I232:A, F234:A, G235:A, L236:A, I238:A, W239:A',
    'channel24: Y76:A, W83:A, T91:A, R94:A, Y95:A, W98:A, G130:A, V133:A, M134:A, F137:A, G138:A, G141:A, E142:A, W197:A, RET301:A',
    'channel25: V102:A, P103:A, I106:A, K126:A, L127:A, V129:A, G130:A, S131:A, M162:A, W197:A',
    'channel26: Y110:A, A120:A, G121:A, S122:A, K125:A, K126:A, M162:A, L166:A, G169:A, E170:A, G171:A, Y186:A, M189:A, M190:A, I193:A',
    'channel27: Y110:A, S122:A, K125:A, K126:A, M162:A, E165:A, L166:A, A168:A, G169:A, E170:A, Y186:A, M189:A, M190:A, I193:A',
    'channel28: Y76:A, W83:A, R94:A, W98:A, G130:A, V133:A, M134:A, F137:A, G141:A, W197:A, Y200:A, Y204:A, L219:A, Y223:A, RET301:A',
    'channel29: Y76:A, R80:A, W83:A, R94:A, W98:A, G130:A, V133:A, M134:A, F137:A, W197:A, Y200:A, L219:A, Y223:A, RET301:A',
    'channel30: Y76:A, W98:A, G130:A, V133:A, M134:A, F137:A, G196:A, W197:A, I199:A, Y200:A, Y223:A, L225:A, A226:A, D227:A, V229:A, N230:A, RET301:A',
    'channel31: Y110:A, L113:A, A114:A, A120:A, S122:A, K125:A, K126:A, L128:A, V129:A, Y161:A, M162:A, E165:A, L166:A, G169:A, E170:A, G171:A, A174:A, V182:A, Y186:A',
    'channel32: Y76:A, W98:A, G130:A, V133:A, M134:A, F137:A, G196:A, W197:A, I199:A, Y200:A, Y223:A, A226:A, D227:A, V229:A, N230:A, RET301:A',
    'channel33: Y76:A, R80:A, W83:A, I84:A, R94:A, W98:A, G130:A, V133:A, M134:A, F137:A, W197:A, Y200:A, A216:A, L219:A, N220:A, Y223:A, RET301:A',
    'channel34: Y76:A, W98:A, G130:A, V133:A, M134:A, F137:A, G196:A, W197:A, I199:A, Y200:A, G203:A, I222:A, Y223:A, L225:A, A226:A, D227:A, V229:A, N230:A, RET301:A',
    'channel35: Y76:A, W98:A, G130:A, V133:A, M134:A, F137:A, I192:A, G196:A, W197:A, I199:A, Y200:A, Y223:A, A226:A, D227:A, V229:A, N230:A, L233:A, F234:A, I237:A, RET301:A',
    'channel36: Y110:A, S122:A, K125:A, K126:A, L128:A, V129:A, Y161:A, M162:A, E165:A, L166:A, G169:A, E170:A, Y186:A, M189:A, M190:A, I193:A']



II. Detection of surface cavities in multi-model PDBs
===============================================================================


In this tutorial, we demonstrate how to identify and characterize surface
cavities  in the substrate-bound structure of *S. aureus* Sortase A using
the NMR structure with PDB ID ``2KID``. Sortase A is a membrane-associated
transpeptidase essential for bacterial virulence, and its substrate-binding
region provides a useful example of a shallow surface cavity located near
the catalytic site.

We first load the protein structure and select only ``chain A``, which
corresponds to the Sortase A protein. This step removes the bound substrate
peptide from the analysis, allowing the cavity detection procedure to
identify the surface groove that accommodates the substrate. 


.. ipython:: python
   :verbatim:

   PDB_ID = '2KID'
   atoms = parsePDB(PDB_ID).select('protein and chain A')


.. parsed-literal::

   @> Connecting wwPDB FTP server RCSB PDB (USA).
   @> Downloading PDB files via FTP failed, trying HTTP.
   @> 2kid downloaded (2kid.pdb.gz)
   @> PDB download via HTTP completed (1 downloaded, 0 failed).
   @> 2437 atoms and 20 coordinate set(s) were parsed in 0.23s.


Surface cavities are calculated for all available NMR models using
:func:`.calcSurfaceCavitiesMultipleFrames`. We used ``inner_radius=1.5`` controlling
the detection of accessible surface cavities. 
The results are also saved as separate PQR files for individual cavities
when ``separate`` parameter is set.


.. ipython:: python
   :verbatim:

   cavities, surface = calcSurfaceCavitiesMultipleFrames(atoms,
			inner_radius=1.5, 
			output_path=PDB_ID+'_CAV_', separate=True)   

.. parsed-literal::

   @> Frame/model: 0
   @> Frame/model: 3
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 2360 atoms with 30586 homogeneous balls of radius 1.20 Å in 0.17s.
   @> Substituted 2360 atoms with 30586 homogeneous balls of radius 1.20 Å in 0.15s.
   @> Delaunay tessellation of 30586 points constructed in 1.36s.
   @> Delaunay tessellation of 30586 points constructed in 1.42s.
   @> Surface and inner simplices filtered in 0.45s.
   @> Surface and inner simplices filtered in 0.42s.
   @> Surface cavities: 244 found, 26 deeper than min_depth=1.5 Å and kept, in 0.31s.
   @> Surface cavities: 238 found, 26 deeper than min_depth=1.5 Å and kept, in 0.31s.
   @> Returning surface cavities
   @> Saving multiple surface cavities to directory ..
   @> Surface cavity calculation completed in 2.43s.
   @> Frame/model: 1
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Returning surface cavities
   @> Saving multiple surface cavities to directory ..
   @> Surface cavity calculation completed in 2.46s.
   @> Frame/model: 4
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 2360 atoms with 30586 homogeneous balls of radius 1.20 Å in 0.10s.
   @> Substituted 2360 atoms with 30586 homogeneous balls of radius 1.20 Å in 0.11s.
   @> Delaunay tessellation of 30586 points constructed in 1.28s.
   @> Delaunay tessellation of 30586 points constructed in 1.36s.
   @> Surface and inner simplices filtered in 0.43s.
   @> Surface and inner simplices filtered in 0.40s.
   @> Surface cavities: 232 found, 32 deeper than min_depth=1.5 Å and kept, in 0.30s.
   @> Surface cavities: 220 found, 28 deeper than min_depth=1.5 Å and kept, in 0.29s.
   @> Returning surface cavities
   @> Saving multiple surface cavities to directory ..
   @> Surface cavity calculation completed in 2.26s.
   @> Frame/model: 2
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Returning surface cavities
   @> Saving multiple surface cavities to directory ..
   @> Surface cavity calculation completed in 2.31s.
   ..
   @> Frame/model: 17
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 2360 atoms with 30586 homogeneous balls of radius 1.20 Å in 0.11s.
   @> Delaunay tessellation of 30586 points constructed in 1.21s.
   @> Surface and inner simplices filtered in 0.43s.
   @> Delaunay tessellation of 30586 points constructed in 1.23s.
   @> Surface cavities: 221 found, 23 deeper than min_depth=1.5 Å and kept, in 0.30s.
   @> Returning surface cavities
   @> Saving multiple surface cavities to directory ..
   @> Surface cavity calculation completed in 2.19s.
   @> Surface and inner simplices filtered in 0.40s.
   @> Frame/model: 18
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 2360 atoms with 30586 homogeneous balls of radius 1.20 Å in 0.10s.
   @> Surface cavities: 255 found, 25 deeper than min_depth=1.5 Å and kept, in 0.33s.
   @> Returning surface cavities
   @> Saving multiple surface cavities to directory ..
   @> Surface cavity calculation completed in 2.22s.
   @> Delaunay tessellation of 30586 points constructed in 1.22s.
   @> Surface and inner simplices filtered in 0.38s.
   @> Surface cavities: 216 found, 25 deeper than min_depth=1.5 Å and kept, in 0.28s.
   @> Returning surface cavities
   @> Saving multiple surface cavities to directory ..
   @> Surface cavity calculation completed in 2.13s.
   @> Frame/model: 19
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 2360 atoms with 30586 homogeneous balls of radius 1.20 Å in 0.10s.
   @> Delaunay tessellation of 30586 points constructed in 1.17s.
   @> Surface and inner simplices filtered in 0.39s.
   @> Surface cavities: 230 found, 25 deeper than min_depth=1.5 Å and kept, in 0.27s.
   @> Returning surface cavities
   @> Saving multiple surface cavities to directory ..
   @> Surface cavity calculation completed in 2.07s.


Next, we generate a VMD molecular model of the protein, which can be used
together with :func:`.showSurfaceCavities` to visualize the detected
cavities and the protein surface. 

.. ipython:: python
   :verbatim:

   vmd_path = '/usr/local/bin/vmd'
   model = getVmdModel(vmd_path, atoms)

We then extract quantitative parameters for each cavity, such as its 
size and geometric descriptors, using
:func:`.getSurfaceCavityParametersMultipleFrames`. 

.. parsed-literal::

   @> Model created successfully.

In addition, we identify the amino acid residues surrounding each cavity
with :func:`.getSurfaceCavityResidueNamesMultipleFrames`, which allows the
detected cavities to be related to the substrate-binding and catalytic
regions of Sortase A.

.. ipython:: python
   :verbatim:

   parameters = getSurfaceCavityParametersMultipleFrames(cavities, 
			param_file_name=PDB_ID+'_param')

.. parsed-literal::

   @> Model/frame: 0
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	599.86 		6.67 		360
   @> cavity 1: 	521.74 		4.69 		308
   @> cavity 2: 	484.91 		6.86 		156
   @> cavity 3: 	389.91 		4.73 		297
   @> cavity 4: 	369.42 		8.6 		181
   @> cavity 5: 	343.3 		6.07 		160
   @> cavity 6: 	132.38 		3.98 		127
   @> cavity 7: 	96.49 		3.23 		70
   @> cavity 8: 	81.99 		3.14 		52
   @> cavity 9: 	80.0 		2.24 		95
   @> cavity 10: 	77.69 		1.86 		39
   @> cavity 11: 	67.58 		3.14 		66
   @> cavity 12: 	65.45 		3.93 		58
   @> cavity 13: 	54.74 		1.9 		43
   @> Model/frame: 1
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	1584.91 		7.96 		980
   @> cavity 1: 	303.16 		4.19 		215
   @> cavity 2: 	211.44 		6.01 		110
   @> cavity 3: 	180.99 		4.45 		105
   @> cavity 4: 	174.84 		6.11 		71
   @> cavity 5: 	164.29 		2.28 		157
   @> cavity 6: 	123.68 		3.5 		59
   @> cavity 7: 	118.85 		3.8 		63
   @> cavity 8: 	85.18 		5.62 		45
   @> cavity 9: 	78.4 		2.06 		75
   @> cavity 10: 	66.94 		2.9 		48
   @> cavity 11: 	64.48 		4.45 		78
   @> cavity 12: 	60.96 		1.59 		60
   @> Model/frame: 2
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	1041.57 		7.92 		504
   @> cavity 1: 	704.89 		5.96 		451
   @> cavity 2: 	409.74 		6.44 		222
   @> cavity 3: 	408.56 		8.82 		150
   @> cavity 4: 	282.52 		6.58 		161
   @> cavity 5: 	242.03 		6.5 		151
   @> cavity 6: 	135.95 		2.56 		110
   @> cavity 7: 	85.45 		2.02 		79
   @> cavity 8: 	74.31 		3.18 		74
   @> cavity 9: 	54.49 		1.64 		40
   @> Model/frame: 3
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	743.6 		11.04 		488
   @> cavity 1: 	645.63 		4.06 		422
   @> cavity 2: 	533.04 		7.76 		267
   @> cavity 3: 	498.7 		5.82 		313
   @> cavity 4: 	316.88 		7.25 		154
   @> cavity 5: 	186.54 		2.92 		120
   @> cavity 6: 	120.94 		5.19 		80
   @> cavity 7: 	113.18 		6.39 		61
   @> cavity 8: 	104.42 		3.87 		91
   @> cavity 9: 	67.45 		1.69 		78
   @> cavity 10: 	65.77 		2.11 		41
   @> cavity 11: 	53.98 		2.47 		37
   @> Model/frame: 4
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	1013.8 		7.65 		536
   @> cavity 1: 	538.61 		3.94 		445
   @> cavity 2: 	427.02 		5.2 		215
   @> cavity 3: 	316.25 		5.37 		154
   @> cavity 4: 	297.92 		5.47 		207
   @> cavity 5: 	207.6 		4.57 		164
   @> cavity 6: 	185.04 		10.36 		109
   @> cavity 7: 	121.18 		2.87 		103
   @> cavity 8: 	83.46 		1.78 		42
   @> cavity 9: 	74.88 		2.99 		61
   @> cavity 10: 	66.31 		2.78 		74
   @> cavity 11: 	53.17 		2.46 		38
   @> Model/frame: 5
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	1313.24 		7.53 		807
   @> cavity 1: 	420.76 		4.29 		302
   @> cavity 2: 	359.53 		6.31 		257
   @> cavity 3: 	285.6 		6.13 		197
   @> cavity 4: 	270.55 		6.46 		165
   @> cavity 5: 	105.65 		2.67 		112
   @> cavity 6: 	87.06 		2.05 		45
   @> cavity 7: 	86.81 		2.75 		35
   @> cavity 8: 	71.32 		3.64 		59
   @> Model/frame: 6
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	1894.14 		6.23 		1052
   @> cavity 1: 	473.35 		5.12 		268
   @> cavity 2: 	310.75 		5.07 		287
   @> cavity 3: 	165.67 		4.84 		99
   @> cavity 4: 	159.91 		8.57 		54
   @> cavity 5: 	135.33 		3.52 		78
   @> cavity 6: 	79.55 		4.07 		57
   @> cavity 7: 	64.64 		2.39 		49
   @> cavity 8: 	51.42 		2.75 		44
   @> Model/frame: 7
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	1032.11 		5.72 		593
   @> cavity 1: 	608.6 		7.25 		379
   @> cavity 2: 	423.08 		4.78 		214
   @> cavity 3: 	246.54 		7.99 		134
   @> cavity 4: 	239.05 		2.36 		141
   @> cavity 5: 	154.19 		3.43 		100
   @> cavity 6: 	150.54 		2.61 		99
   @> cavity 7: 	141.99 		4.32 		133
   @> cavity 8: 	115.31 		3.04 		119
   @> cavity 9: 	98.74 		2.56 		64
   @> cavity 10: 	73.49 		2.38 		58
   @> Model/frame: 8
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	813.97 		6.84 		650
   @> cavity 1: 	582.41 		6.27 		320
   @> cavity 2: 	440.57 		4.58 		216
   @> cavity 3: 	378.83 		6.9 		141
   @> cavity 4: 	326.96 		5.94 		208
   @> cavity 5: 	160.36 		5.11 		50
   @> cavity 6: 	108.76 		6.62 		58
   @> cavity 7: 	69.88 		2.47 		61
   @> cavity 8: 	67.22 		1.67 		47
   @> Model/frame: 9
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	882.04 		9.36 		364
   @> cavity 1: 	439.5 		3.4 		309
   @> cavity 2: 	421.07 		5.14 		266
   @> cavity 3: 	353.77 		7.89 		129
   @> cavity 4: 	336.98 		6.57 		215
   @> cavity 5: 	329.21 		5.73 		199
   @> cavity 6: 	253.94 		9.67 		192
   @> cavity 7: 	134.6 		3.55 		116
   @> cavity 8: 	104.39 		2.85 		37
   @> cavity 9: 	92.26 		2.23 		54
   @> cavity 10: 	67.0 		2.6 		61
   @> cavity 11: 	62.02 		3.27 		69
   @> cavity 12: 	59.87 		1.75 		47
   @> cavity 13: 	58.8 		2.42 		43
   @> cavity 14: 	58.26 		2.93 		46
   @> Model/frame: 10
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	1504.0 		7.49 		826
   @> cavity 1: 	795.94 		9.65 		441
   @> cavity 2: 	569.56 		4.85 		361
   @> cavity 3: 	225.63 		5.36 		149
   @> cavity 4: 	147.11 		2.94 		111
   @> cavity 5: 	69.97 		2.38 		70
   @> Model/frame: 11
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	1186.46 		8.9 		538
   @> cavity 1: 	621.2 		7.71 		364
   @> cavity 2: 	566.42 		7.08 		399
   @> cavity 3: 	281.48 		9.79 		205
   @> cavity 4: 	216.96 		2.89 		177
   @> cavity 5: 	212.12 		3.74 		132
   @> cavity 6: 	111.26 		2.34 		90
   @> cavity 7: 	89.38 		4.16 		53
   @> cavity 8: 	83.61 		2.76 		86
   @> cavity 9: 	62.0 		1.95 		49
   @> Model/frame: 12
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	1120.03 		10.15 		651
   @> cavity 1: 	882.51 		7.9 		516
   @> cavity 2: 	456.85 		5.46 		271
   @> cavity 3: 	215.66 		3.72 		122
   @> cavity 4: 	207.94 		10.73 		120
   @> cavity 5: 	172.32 		2.87 		133
   @> cavity 6: 	109.18 		2.52 		71
   @> cavity 7: 	91.33 		2.13 		50
   @> Model/frame: 13
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	924.55 		12.8 		533
   @> cavity 1: 	707.43 		4.08 		374
   @> cavity 2: 	442.47 		4.86 		261
   @> cavity 3: 	389.81 		5.27 		227
   @> cavity 4: 	203.99 		5.63 		147
   @> cavity 5: 	201.09 		4.35 		155
   @> cavity 6: 	131.22 		6.65 		79
   @> cavity 7: 	124.74 		2.54 		78
   @> cavity 8: 	56.42 		2.72 		45
   @> cavity 9: 	53.6 		3.18 		30
   @> Model/frame: 14
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	1510.06 		6.51 		720
   @> cavity 1: 	1375.79 		7.96 		682
   @> cavity 2: 	401.6 		5.37 		234
   @> cavity 3: 	182.55 		5.62 		66
   @> cavity 4: 	177.66 		2.57 		120
   @> cavity 5: 	87.75 		2.78 		67
   @> cavity 6: 	78.44 		9.66 		56
   @> Model/frame: 15
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	1306.11 		7.56 		810
   @> cavity 1: 	552.46 		8.72 		323
   @> cavity 2: 	526.45 		6.51 		396
   @> cavity 3: 	324.56 		5.79 		197
   @> cavity 4: 	277.43 		3.4 		148
   @> cavity 5: 	251.42 		4.28 		162
   @> cavity 6: 	176.09 		3.27 		135
   @> cavity 7: 	150.08 		3.24 		105
   @> cavity 8: 	76.55 		2.33 		47
   @> cavity 9: 	67.61 		2.09 		55
   @> Model/frame: 16
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	1523.36 		6.69 		874
   @> cavity 1: 	788.97 		5.88 		370
   @> cavity 2: 	310.64 		4.03 		183
   @> cavity 3: 	304.92 		5.95 		213
   @> cavity 4: 	136.76 		3.72 		123
   @> cavity 5: 	98.55 		2.16 		73
   @> cavity 6: 	93.83 		3.51 		104
   @> cavity 7: 	73.64 		1.64 		59
   @> cavity 8: 	72.21 		2.3 		54
   @> cavity 9: 	70.23 		2.41 		38
   @> Model/frame: 17
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	1136.94 		7.07 		599
   @> cavity 1: 	594.88 		6.76 		399
   @> cavity 2: 	333.92 		6.15 		240
   @> cavity 3: 	303.18 		6.35 		229
   @> cavity 4: 	276.13 		4.08 		174
   @> cavity 5: 	183.9 		5.73 		109
   @> cavity 6: 	169.33 		3.03 		127
   @> cavity 7: 	120.97 		3.08 		66
   @> cavity 8: 	112.92 		3.02 		90
   @> Model/frame: 18
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	1329.65 		7.74 		620
   @> cavity 1: 	942.27 		7.62 		585
   @> cavity 2: 	536.1 		6.4 		378
   @> cavity 3: 	400.23 		5.67 		248
   @> cavity 4: 	206.96 		12.19 		114
   @> cavity 5: 	153.19 		2.51 		107
   @> cavity 6: 	153.08 		2.18 		118
   @> cavity 7: 	104.02 		3.82 		43
   @> cavity 8: 	103.13 		3.3 		89
   @> cavity 9: 	85.65 		1.84 		68
   @> cavity 10: 	74.67 		2.51 		68
   @> Model/frame: 19
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	1462.62 		7.46 		873
   @> cavity 1: 	1243.19 		4.59 		614
   @> cavity 2: 	546.92 		6.01 		302
   @> cavity 3: 	273.83 		6.53 		105
   @> cavity 4: 	154.38 		5.37 		79
   @> cavity 5: 	151.99 		2.73 		135
   @> cavity 6: 	58.85 		3.6 		52
   @> cavity 7: 	57.84 		2.29 		47


.. ipython:: python
   :verbatim:

   residues = getSurfaceCavityResidueNamesMultipleFrames(atoms, cavities, 
			surface, residues_file_name=PDB_ID+'_resAA')

.. parsed-literal::

   @> Surface cavity residues were saved to: 2KID_resAA_model0_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: 2KID_resAA_model1_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: 2KID_resAA_model2_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: 2KID_resAA_model3_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: 2KID_resAA_model4_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: 2KID_resAA_model5_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: 2KID_resAA_model6_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: 2KID_resAA_model7_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: 2KID_resAA_model8_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: 2KID_resAA_model9_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: 2KID_resAA_model10_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: 2KID_resAA_model11_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: 2KID_resAA_model12_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: 2KID_resAA_model13_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: 2KID_resAA_model14_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: 2KID_resAA_model15_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: 2KID_resAA_model16_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: 2KID_resAA_model17_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: 2KID_resAA_model18_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: 2KID_resAA_model19_Residues_All_surface_cavities.txt


.. ipython:: python
   :verbatim:

   showSurfaceCavities(surface[0], model=model, show_surface=True)


.. figure:: images/cavitracer_figure23.jpg
   :scale: 50 %


Finally, the generated cavity PQR files are collected and passed to 
:func:`.calcSurfaceCavityOverlaps`. This step compares cavities detected
across different NMR models and produces an overlap representation, which
can be used to identify surface cavities that are consistently present
across the conformational ensemble.

.. ipython:: python
   :verbatim:

   import glob
   pqr_files_cavities = glob.glob(PDB_ID+"_CAV_?.pqr") + glob.glob(PDB_ID+"_CAV_??.pqr")
   calcSurfaceCavityOverlaps(pqr_files=pqr_files_cavities, 
	output_file_name=PDB_ID+'surface_cavity_overlap.pdb', max_proc=4)


.. parsed-literal::

   @> Number of PQR files: 20
   @> Resolution: 0.5
   @> max_proc: 4
   @> Calculating overlaps using 4 processes.
   @> 1942 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1979 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1988 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2148 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2147 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1751 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2034 atoms and 1 coordinate sets were parsed in 0.02s.
   @> 2152 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2066 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2012 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2033 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2207 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1958 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1934 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2091 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2093 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2438 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2378 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1929 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1945 atoms and 1 coordinate sets were parsed in 0.01s.
   @> Overlap written to: 2KIDsurface_cavity_overlap.pdb
   @> Number of occupied overlap voxels: 118035
   
   '2KIDsurface_cavity_overlap.pdb'


The final outcome can be displayed in VMD_. In the example shown below, 
the surface cavities present in at least 75% of the analyzed NMR models 
are displayed.

.. figure:: images/cavitracer_figure24.jpg
   :scale: 50 %


III. Detection of pores in multi-model PDBs
===============================================================================


In order to identify pores within protein structure, we will use a multi-model
PDB file with ~20 frames from the MD simulation. The structure belongs to
the vesicular monoamine transporter VMAT2, which contains 460 residues. The
full trajectory can be found in the WatFinder tutorial files.

First, we will upload the multi-model PDB ``case_study2_ev10_multi.pdb`` which
can be found in the tutorial files.

.. ipython:: python
   :verbatim:

   pdb_multi = parsePDB('case_study2_ev10_multi.pdb')

.. parsed-literal::

   @> 5986 atoms and 22 coordinate set(s) were parsed in 0.29s.

.. ipython:: python
   :verbatim:

   pdb_multi

.. parsed-literal::

   <AtomGroup: case_study2_ev10_multi (5986 atoms; active #0 of 22 coordsets)>

To identify pores, we should use :func:.`calcChannelsMultipleFrames`
function first and identify channels. Results will be saved with ``'ch_multi_'``
prefix, each channel in a separate file (``separate=True``). Additionally,
the requirement to obtain information about pores is
``return_details=True``. We will also apply ``max_proc=4`` to use four
processors for calculations.

.. ipython:: python
   :verbatim:

   channels, surface, details = calcChannelsMultipleFrames(pdb_multi,
				    inner_radius=0.85,
                                    output_path='ch_multi_', separate=True,
                                    return_details=True, max_proc=4)

.. parsed-literal::

   @> Frame/model: 0
   @> Frame/model: 2
   @> Frame/model: 4
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Frame/model: 6
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 A in 0.28s.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 A in 0.28s.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 A in 0.28s.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 A in 0.29s.
   @> Delaunay tessellation of 77434 points constructed in 3.97s.
   @> Delaunay tessellation of 77434 points constructed in 4.06s.
   @> Delaunay tessellation of 77434 points constructed in 4.08s.
   @> Delaunay tessellation of 77434 points constructed in 4.14s.
   @> Surface and inner simplices filtered in 3.77s.
   @> Surface and inner simplices filtered in 3.60s.
   @> Surface and inner simplices filtered in 3.72s.
   @> Surface and inner simplices filtered in 3.81s.
   @> 11 surface cavities detected and filtered in 0.60s.
   @> 9 surface cavities detected and filtered in 0.59s.
   @> 9 surface cavities detected and filtered in 0.65s.
   @> 8 surface cavities detected and filtered in 0.58s.
   @> Channel pathfinding (graph Dijkstra) over 8 cavities completed in 1.20s.
   @> Detected 32 channels.
   @> Saving multiple results to directory ..
   @> Channel calculation completed in 10.01s.
   @> Frame/model: 3
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 A in 0.27s.
   @> Channel pathfinding (graph Dijkstra) over 9 cavities completed in 1.66s.
   @> Detected 37 channels.
   @> Saving multiple results to directory ..
   @> Channel calculation completed in 10.43s.
   @> Frame/model: 5
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Channel pathfinding (graph Dijkstra) over 11 cavities completed in 1.97s.
   @> Detected 36 channels.
   @> Saving multiple results to directory ..
   @> Channel pathfinding (graph Dijkstra) over 9 cavities completed in 1.94s.
   @> Detected 42 channels.
   @> Saving multiple results to directory ..
   @> Channel calculation completed in 10.68s.
   @> Frame/model: 1
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 A in 0.27s.
   @> Channel calculation completed in 10.73s.
   @> Frame/model: 7
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 A in 0.30s.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 A in 0.30s.
   @> Delaunay tessellation of 77434 points constructed in 3.42s.
   @> Delaunay tessellation of 77434 points constructed in 3.67s.
   ..
   ..
   @> Frame/model: 19
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 A in 0.28s.
   @> Channel pathfinding (graph Dijkstra) over 10 cavities completed in 1.66s.
   @> Detected 32 channels.
   @> Saving multiple results to directory ..
   @> Channel calculation completed in 9.65s.
   @> Delaunay tessellation of 77434 points constructed in 3.26s.
   @> Delaunay tessellation of 77434 points constructed in 3.28s.
   @> Surface and inner simplices filtered in 3.06s.
   @> Surface and inner simplices filtered in 2.89s.
   @> 5 surface cavities detected and filtered in 0.56s.
   @> 7 surface cavities detected and filtered in 0.56s.
   @> Channel pathfinding (graph Dijkstra) over 5 cavities completed in 1.66s.
   @> Detected 31 channels.
   @> Saving multiple results to directory ..
   @> Channel calculation completed in 8.88s.
   @> Channel pathfinding (graph Dijkstra) over 7 cavities completed in 1.49s.
   @> Detected 31 channels.
   @> Saving multiple results to directory ..
   @> Channel calculation completed in 8.59s.

After channels, :func:.`calcPoresFromChannelsMultipleFrames` function can be
applied to reconstruct pores.

.. ipython:: python
   :verbatim:

   pores = calcPoresFromChannelsMultipleFrames(channels, details,
                                                min_end_to_end=45,
                                                output_path='pores_multi_',
                                                separate=True,
                                                min_bottleneck=0.65,
                                                max_proc=4)

.. parsed-literal::

   @> Frame/model: 0
   @> Frame/model: 2
   @> Frame/model: 4
   @> Frame/model: 3
   @> Frame/model: 5
   @> Frame/model: 6
   @> Frame/model: 1
   @> Frame/model: 8
   @> Frame/model: 9
   @> Frame/model: 7
   @> Frame/model: 10
   @> Frame/model: 12
   @> Frame/model: 11
   @> Frame/model: 14
   @> Frame/model: 15
   @> Frame/model: 16
   @> Frame/model: 18
   @> Frame/model: 17
   @> Frame/model: 13
   @> Frame/model: 20
   @> Frame/model: 19


.. ipython:: python
   :verbatim:

   pores

.. parsed-literal::

   [[<prody.proteins.channels.Channel at 0x77083df351b0>,
     <prody.proteins.channels.Channel at 0x77083df356c0>,
     <prody.proteins.channels.Channel at 0x77083df359f0>,
     <prody.proteins.channels.Channel at 0x77083df36050>,
     <prody.proteins.channels.Channel at 0x77083df36980>,
     <prody.proteins.channels.Channel at 0x77083df37280>,
     <prody.proteins.channels.Channel at 0x77083df37490>,
     <prody.proteins.channels.Channel at 0x77083df376d0>,
     <prody.proteins.channels.Channel at 0x77083df37760>,
     <prody.proteins.channels.Channel at 0x77083df34ac0>,
     <prody.proteins.channels.Channel at 0x77083df37a00>],
    [],
    [<prody.proteins.channels.Channel at 0x77089a9d5db0>,
     <prody.proteins.channels.Channel at 0x77089a9d4f10>,
     <prody.proteins.channels.Channel at 0x77089a9d6a70>,
     <prody.proteins.channels.Channel at 0x77089a9d7730>],
    [<prody.proteins.channels.Channel at 0x77089a9d77f0>,
     <prody.proteins.channels.Channel at 0x77089a9d5e10>,
     <prody.proteins.channels.Channel at 0x77089a9d72e0>,
     <prody.proteins.channels.Channel at 0x77083e0aa800>,
     <prody.proteins.channels.Channel at 0x77083e0aab90>,
     <prody.proteins.channels.Channel at 0x77083e0ab040>,
     <prody.proteins.channels.Channel at 0x77083e0ab100>,
     <prody.proteins.channels.Channel at 0x77083e0ab2b0>,
     <prody.proteins.channels.Channel at 0x77083e0ab490>,
     <prody.proteins.channels.Channel at 0x77083e0ab520>,
     <prody.proteins.channels.Channel at 0x77083e0a9540>,
     <prody.proteins.channels.Channel at 0x77083e0ab880>,
     <prody.proteins.channels.Channel at 0x77083e0abac0>,
     <prody.proteins.channels.Channel at 0x77083e0abbb0>,
     <prody.proteins.channels.Channel at 0x77083e0abc10>,
     <prody.proteins.channels.Channel at 0x77083e0abeb0>,
     <prody.proteins.channels.Channel at 0x77083e0a8490>,
     <prody.proteins.channels.Channel at 0x77083e0aa3b0>],
    [],
    [<prody.proteins.channels.Channel at 0x77083e857be0>,
     <prody.proteins.channels.Channel at 0x77083e857370>,
     <prody.proteins.channels.Channel at 0x77083e8314b0>,
     <prody.proteins.channels.Channel at 0x77083e833c40>,
     ..
     ..
    [<prody.proteins.channels.Channel at 0x77083e8319f0>,
     <prody.proteins.channels.Channel at 0x77083df35840>,
     <prody.proteins.channels.Channel at 0x77083df35cf0>,
     <prody.proteins.channels.Channel at 0x77083df36a40>,
     <prody.proteins.channels.Channel at 0x77083df37130>,
     <prody.proteins.channels.Channel at 0x77083df374f0>,
     <prody.proteins.channels.Channel at 0x77083df37ac0>,
     <prody.proteins.channels.Channel at 0x77083df349d0>,
     <prody.proteins.channels.Channel at 0x77083df37d00>,
     <prody.proteins.channels.Channel at 0x77083e857640>,
     <prody.proteins.channels.Channel at 0x77083e854ca0>,
     <prody.proteins.channels.Channel at 0x77083e8571c0>]]

To obtain information about residues that are forming pores, use
:func:.`getPoreResidueNamesMultipleFrames` function. When applying
``one_letter_aa=True`` all the residues will be written in a one-letter
code.

.. ipython:: python
   :verbatim:

   getPoreResidueNamesMultipleFrames(pdb_multi, pores, one_letter_aa=True)

.. parsed-literal::

   @> Model: 0
   @> 6651 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6671 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6641 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6656 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6661 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6736 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6736 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6751 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6756 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6756 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6771 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> Model: 1
   @> Model: 2
   @> 6506 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6641 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6651 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6656 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> Model: 3
   @> 6516 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6521 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6561 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6661 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6546 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6576 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6586 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6581 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6586 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6596 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6621 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6616 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6626 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6621 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6611 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6626 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6756 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6756 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> Model: 4
   @> Model: 5
   @> 6511 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6561 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6536 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6556 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6571 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6621 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6596 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6616 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6576 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6591 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6626 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6601 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6621 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6576 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6591 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6626 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6601 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6621 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6611 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6626 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6661 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6636 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6656 atoms and 1 coordinate set(s) were parsed in 0.06s.
   ..
   ..
   @> Model: 18
   @> 6531 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6531 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6671 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6531 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6476 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6476 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6686 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6471 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6471 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6681 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6486 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6486 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> Model: 19
   @> 6486 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6491 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6656 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6496 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> Model: 20
   @> 6546 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6566 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6551 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6571 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6551 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6571 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6576 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6596 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6621 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6641 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6621 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6641 atoms and 1 coordinate set(s) were parsed in 0.06s.

   [['pore0: L30, N34, L36, L37, V39, V40, I43, F135, K138, A139, Q142, F176, S179, S181, L185, R189, Q192, S196, S200, M204, D214, E215, R217, G218, M221, G222, L225, L228, K248, F252, G349, I350, H353, R357, D399, S400, S401, M403, P404, D426, F429, C430, Y433',
     'pore1: L30, N34, L36, L37, V39, V40, I43, F135, K138, A139, Q142, F176, S179, S181, L185, R189, Q192, S196, S200, M204, R217, G218, M221, G222, L225, L228, K248, F252, F348, G349, A352, H353, K354, M355, G356, R357, L359, C360, D399, S400, S401, M403, P404, D426, F429, C430, Y433',
     'pore2: L30, N34, L36, L37, V39, V40, I43, F135, K138, A139, Q142, F176, A177, F178, S179, L185, R189, Q192, S196, S200, M204, D214, E215, R217, G218, M221, G222, L225, L228, K248, F252, G349, I350, H353, R357, D399, S400, S401, M403, P404, D426, F429, C430, Y433',
     'pore3: L30, N34, L36, L37, V39, V40, I43, F135, K138, A139, Q142, F176, A177, F178, S179, L185, R189, Q192, S196, S200, M204, R217, G218, M221, G222, L225, L228, K248, F252, F348, G349, A352, H353, K354, G356, R357, L359, C360, L363, D399, S400, S401, M403, P404, D426, F429, C430, Y433',
     'pore4: L30, N34, L36, L37, V39, V40, I43, F135, K138, A139, Q142, F176, A177, F178, S179, L185, R189, Q192, S196, S200, M204, R217, G218, M221, G222, L225, L228, K248, F252, F348, G349, A352, H353, K354, M355, G356, R357, L359, C360, D399, S400, S401, M403, P404, D426, F429, C430, Y433',
     'pore5: L30, N34, L37, T38, V40, V41, P42, I44, P45, F135, K138, Q142, R189, S196, S200, M204, D214, E215, R217, G218, M221, G222, L225, L228, V232, P236, S240, Y243, E244, E312, I317, W318, K327, L330, G331, A333, F334, G349, I350, H353, R357, D399, S400, S401, M403, P404, D426, F429',
     'pore6: L30, N34, L37, T38, V40, V41, P42, I44, P45, F135, K138, Q142, R189, S196, S200, M204, D214, E215, R217, G218, M221, G222, L225, L228, V232, P236, S240, Y243, E312, I317, W318, M323, K327, L330, G331, A333, F334, G349, I350, H353, R357, D399, S400, S401, M403, P404, D426, F429',
     'pore7: L30, N34, L37, T38, V40, V41, P42, I43, I44, P45, S46, F135, K138, Q142, R189, S196, S200, M204, D214, E215, R217, G218, M221, G222, L225, L228, V232, P236, S240, Y243, E312, I317, W318, K327, L330, G331, A333, F334, G349, I350, H353, R357, D399, S400, S401, M403, P404, D426, F429',
     'pore8: L30, N34, L37, T38, V40, V41, P42, I44, P45, F135, K138, Q142, R189, S196, S200, M204, R217, G218, M221, G222, L225, L228, V232, P236, S240, Y243, E244, E312, I317, W318, K327, L330, G331, A333, F334, F348, G349, A352, H353, K354, M355, G356, R357, L359, C360, D399, S400, S401, M403, P404, D426, F429',
     'pore9: L30, N34, L37, T38, V40, V41, P42, I44, P45, F135, K138, Q142, R189, S196, S200, M204, R217, G218, M221, G222, L225, L228, V232, P236, S240, Y243, E312, I317, W318, M323, K327, L330, G331, A333, F334, F348, G349, A352, H353, K354, M355, G356, R357, L359, C360, D399, S400, S401, M403, P404, D426, F429',
     'pore10: L30, N34, L37, T38, V40, V41, P42, I43, I44, P45, S46, F135, K138, Q142, R189, S196, S200, M204, R217, G218, M221, G222, L225, L228, V232, P236, S240, Y243, E312, I317, W318, K327, L330, G331, A333, F334, F348, G349, A352, H353, K354, M355, G356, R357, L359, C360, D399, S400, S401, M403, P404, D426, F429'],
    [],
    ['pore0: I22, I25, V26, A29, L30, D33, N34, L37, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L225, L228, V232, L270, Q276, I308, E312, L315, I317, Q329, A333, F334, A337, S338, Y341, I381, I385, N388, M403, D426, F429',
     'pore1: L30, D33, N34, L36, L37, V40, I43, K138, Q142, F176, S179, L185, R189, Q192, S196, S200, M204, L225, L228, K248, F252, A297, S300, I301, Y341, R357, W358, A361, M397, V398, D399, S401, M402, M403, I405, Y408, D426, F429, C467, L470',
     'pore2: L30, D33, N34, L36, L37, V40, I43, K138, Q142, F176, S179, S180, S181, L185, R189, Q192, S196, S200, M204, L225, L228, K248, F252, A297, S300, I301, Y341, R357, W358, A361, M397, V398, D399, S401, M402, M403, I405, Y408, D426, F429, C467, L470',
     'pore3: L30, D33, N34, L36, L37, V40, I43, K138, Q142, F176, A177, S179, L185, R189, Q192, S196, S200, M204, L225, L228, K248, T249, F252, A297, S300, I301, Y341, R357, W358, A361, M397, V398, D399, S401, M402, M403, I405, Y408, D426, F429, C467, L470'],
    ['pore0: L30, L37, K138, A139, Q142, G150, T153, N154, Y158, S200, V201, A202, M204, G205, L225, L228, V232, G282, T283, P284, E312, L315, I317, Q329, A333, F334, P336, A337, S338, Y341, I381, I385, N388, M403, V417, G419, S420, Y422, A423, A425, D426, F429, C430, Y433',
     'pore1: L30, L37, K138, A139, Q142, G150, T153, N154, Y158, S200, V201, A202, M204, G205, L225, L228, V232, S279, Q280, K281, G282, E312, L315, I317, Q329, A333, F334, P336, A337, S338, Y341, I381, I385, N388, M403, V417, G419, S420, Y422, A423, A425, D426, F429, C430, Y433',
     'pore2: L30, L37, K138, A139, Q142, G150, T153, N154, R155, G157, Y158, S200, V201, A202, M204, G205, L225, L228, V232, E278, S279, Q280, K281, E312, L315, I317, Q329, A333, F334, P336, A337, S338, Y341, I381, I385, N388, M403, V417, Y418, G419, Y422, A423, A425, D426, F429, C430, Y433',
     'pore3: L30, L37, K138, A139, Q142, G150, T153, N154, Y158, S200, V201, A202, M204, G205, L225, L228, V232, G282, T283, T287, L288, K290, D291, I294, E312, L315, I317, Q329, A333, F334, P336, A337, S338, Y341, I381, I385, N388, M403, V410, H414, V417, G419, S420, V421, Y422, A423, A425, D426, F429, C430, Y433',
     'pore4: L30, L37, K138, A139, Q142, G150, T153, N154, R155, G157, Y158, S200, V201, A202, M204, G205, L225, L228, V232, E278, S279, Q280, K281, E312, L315, A333, F334, P336, A337, S338, Y341, I385, N388, F389, M403, V417, Y418, G419, Y422, A423, A425, D426, F429, C430, Y433',
     'pore5: L30, L37, K138, A139, Q142, G150, T153, N154, Y158, S200, V201, A202, M204, G205, L225, L228, V232, G282, T283, P284, L311, E312, A314, L315, F334, A337, S338, Y341, I375, P376, A378, K379, N380, G383, L384, P387, N388, M403, V417, G419, S420, Y422, A423, A425, D426, F429, C430, Y433',
     'pore6: L30, L37, K138, A139, Q142, G150, T153, N154, Y158, S200, V201, A202, M204, G205, S209, L225, L228, V232, S279, L311, E312, A314, L315, F334, A337, S338, Y341, I375, P376, A378, K379, N380, G383, L384, P387, N388, M403, S416, V417, Y418, G419, Y422, A423, A425, D426, F429, C430, Y433',
     'pore7: L30, L37, K138, A139, Q142, G150, T153, N154, Y158, S200, V201, A202, M204, G205, L225, L228, V232, S279, Q280, K281, G282, L311, E312, A314, L315, F334, A337, S338, Y341, I375, P376, A378, K379, N380, G383, L384, P387, N388, M403, V417, G419, S420, Y422, A423, A425, D426, F429, C430, Y433',
     'pore8: L30, L37, K138, A139, Q142, G150, T153, N154, G157, Y158, S200, V201, A202, M204, G205, S209, L225, L228, V232, E278, S279, L311, E312, A314, L315, F334, A337, S338, Y341, I375, P376, A378, K379, N380, G383, L384, P387, N388, M403, V417, Y418, G419, Y422, A423, A425, D426, F429, C430, Y433',
     'pore9: R17, L30, L37, K138, A139, Q142, G150, T153, N154, G157, Y158, S200, V201, A202, M204, G205, M206, L225, L228, V232, E278, S279, L311, E312, A314, L315, F334, A337, S338, Y341, I375, P376, A378, K379, N380, G383, L384, P387, N388, M403, V417, Y418, G419, Y422, A423, A425, D426, F429, C430, Y433',
     'pore10: L30, L37, K138, A139, Q142, G150, T153, N154, R155, G157, Y158, S200, V201, A202, M204, G205, L225, L228, V232, E278, S279, Q280, K281, L311, E312, A314, L315, F334, A337, S338, Y341, I375, P376, A378, K379, N380, G383, L384, P387, N388, M403, V417, Y418, G419, Y422, A423, A425, D426, F429, C430, Y433',
     'pore11: I25, V26, A29, L30, L37, K138, A139, Q142, P159, I160, I162, F163, S199, S200, G203, M204, L225, L228, V232, Q266, L270, Q271, P272, L311, E312, A314, L315, F334, A337, S338, Y341, I375, P376, A378, K379, N380, G383, L384, P387, N388, M403, D426, F429, C430, Y433',
     'pore12: I25, V26, A29, L30, L37, K138, A139, Q142, P159, I160, I162, F163, S199, S200, G203, M204, L225, L228, V232, Q266, L267, L270, Q271, P272, L311, E312, A314, L315, F334, A337, S338, Y341, I375, P376, A378, K379, N380, G383, L384, P387, N388, M403, D426, F429, C430, Y433',
     'pore13: I25, V26, A29, L30, L37, K138, A139, Q142, P159, I160, I162, F163, S199, S200, G203, M204, L225, L228, V232, Q266, L267, L270, Q271, P272, L311, E312, A314, L315, F334, A337, S338, Y341, I375, P376, A378, K379, N380, G383, L384, P387, N388, M403, D426, F429, C430, Y433',
     'pore14: I25, V26, A29, L30, L37, K138, A139, Q142, P159, I160, I162, F163, S199, S200, G203, M204, L225, L228, V232, Q266, L267, L270, P272, L311, E312, A314, L315, F334, A337, S338, Y341, I375, P376, A378, K379, N380, G383, L384, P387, N388, M403, D426, F429, C430, Y433',
     'pore15: A29, L30, L32, D33, N34, L36, K138, Q142, I162, G165, F166, M169, F170, T173, I174, A177, S196, S199, S200, M204, G218, N219, M221, G222, A224, L225, L228, F252, L253, L255, A256, V259, D262, Q266, T345, F348, G349, I350, H353, G356, R357, D399, S400, S401, M403, P404, D426, F429',
     'pore16: A29, L30, L32, D33, N34, L36, K138, Q142, G150, T153, N154, Y158, I162, G165, F166, M169, F170, T173, I174, A177, S196, S199, S200, V201, A202, M204, G205, L225, L228, F252, L253, L255, A256, V259, D262, Q266, G282, T283, T287, L288, K290, D291, I294, M403, V410, H414, V417, G419, S420, V421, Y422, A423, A425, D426, F429',
     'pore17: A29, L30, L32, D33, N34, L36, K138, Q142, G150, T153, N154, Y158, I162, G165, F166, M169, F170, T173, I174, A177, S196, S199, S200, V201, A202, M204, G205, L225, L228, F252, L253, L255, A256, V259, D262, Q266, G282, T283, T287, L288, K290, D291, I294, M403, V410, H414, V417, G419, S420, V421, Y422, A423, A425, D426, F429'],
    [],
    ['pore0: I22, I25, V26, A29, L30, D33, N34, L37, F135, K138, Q142, P159, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, Q271, P272, L311, E312, L315, P316, I317, M319, A333, L384, I385, N388, D426, F429, Y433',
     'pore1: I22, I25, V26, A29, L30, D33, N34, L37, F135, K138, Q142, P159, I160, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, P272, S273, R274, V275, L311, E312, L315, P316, I317, M319, A333, L384, I385, N388, D426, F429, Y433',
     'pore2: I22, I25, V26, A29, L30, D33, N34, L37, F135, K138, Q142, P159, I160, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, Q271, P272, R274, V275, Q276, L311, E312, L315, P316, I317, M319, A333, L384, I385, N388, D426, F429, Y433',
     'pore3: I22, I25, V26, A29, L30, D33, N34, L37, F135, K138, Q142, Y158, P159, I160, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, P272, R274, V275, Q276, L311, E312, L315, P316, I317, M319, A333, L384, I385, N388, D426, F429, Y433',
     'pore4: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, E120, L124, F135, K138, Q142, P159, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, Q271, P272, E312, P313, P316, I317, W318, F334, D426, F429, Y433',
     'pore5: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, E120, L124, F135, K138, Q142, P159, I160, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, P272, S273, R274, V275, E312, P313, P316, I317, W318, F334, D426, F429, Y433',
     'pore6: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, E120, L124, F135, K138, Q142, P159, I160, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, Q271, P272, R274, V275, Q276, E312, P313, P316, I317, W318, F334, D426, F429, Y433',
     'pore7: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, E120, L124, F135, K138, Q142, Y158, P159, I160, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, P272, R274, V275, Q276, E312, P313, P316, I317, W318, F334, D426, F429, Y433',
     'pore8: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, L124, F135, K138, Q142, P159, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, Q271, P272, E312, P313, A314, L315, P316, I317, W318, F334, D426, F429, Y433',
     'pore9: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, L124, F135, K138, Q142, P159, I160, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, P272, E312, P313, A314, L315, P316, I317, W318, F334, D426, F429, Y433',
     'pore10: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, L124, F135, K138, Q142, P159, I160, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, P272, S273, R274, V275, E312, P313, A314, L315, P316, I317, W318, F334, D426, F429, Y433',
     'pore11: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, L124, F135, K138, Q142, P159, I160, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, Q271, P272, R274, V275, Q276, E312, P313, A314, L315, P316, I317, W318, F334, D426, F429, Y433',
     'pore12: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, L124, F135, K138, Q142, Y158, P159, I160, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, P272, R274, V275, Q276, E312, P313, A314, L315, P316, I317, W318, F334, D426, F429, Y433',
     'pore13: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, L124, F135, K138, Q142, P159, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, Q271, P272, E312, P313, A314, L315, P316, I317, W318, F334, D426, F429, Y433',
     'pore14: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, L124, F135, K138, Q142, P159, I160, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, P272, E312, P313, A314, L315, P316, I317, W318, F334, D426, F429, Y433',
     'pore15: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, L124, F135, K138, Q142, P159, I160, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, P272, S273, R274, V275, E312, P313, A314, L315, P316, I317, W318, F334, D426, F429, Y433',
     'pore16: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, L124, F135, K138, Q142, P159, I160, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, Q271, P272, R274, V275, Q276, E312, P313, A314, L315, P316, I317, W318, F334, D426, F429, Y433',
     'pore17: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, L124, F135, K138, Q142, Y158, P159, I160, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, P272, R274, V275, Q276, E312, P313, A314, L315, P316, I317, W318, F334, D426, F429, Y433',
     'pore18: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, E120, D121, K122, L124, L125, E127, N128, V129, V131, G132, F135, K138, Q142, P159, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, Q271, P272, E312, P313, W318, F334, D426, F429, P437',
     'pore19: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, E120, D121, K122, L124, L125, E127, N128, V129, V131, G132, F135, K138, Q142, P159, I160, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, P272, E312, P313, W318, F334, D426, F429, Y433, P437',
     'pore20: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, E120, D121, K122, L124, L125, E127, N128, V129, V131, G132, F135, K138, Q142, P159, I160, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, P272, S273, R274, V275, E312, P313, W318, F334, D426, F429, Y433, P437',
     'pore21: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, E120, D121, K122, L124, L125, E127, N128, V129, V131, G132, F135, K138, Q142, P159, I160, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, Q271, P272, R274, V275, Q276, E312, P313, W318, F334, D426, F429, Y433, P437',
     'pore22: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, E120, D121, K122, L124, L125, E127, N128, V129, V131, G132, F135, K138, Q142, Y158, P159, I160, I162, F163, S196, S199, S200, G203, M206, L228, V232, Q266, L270, P272, R274, V275, Q276, E312, P313, W318, F334, D426, F429, Y433, P437'],
      ..
      ..
    ['pore0: L30, D33, N34, L37, T38, E120, N128, Q130, V131, K138, Q142, S196, S200, M204, D214, E215, G218, N219, M221, G222, L225, L228, V232, E312, I317, W318, Y341, T345, F348, G349, H353, R357, G396, D399, S400, M403, P404, D426, F429, Y433',
     'pore1: L30, D33, N34, L37, T38, E120, N128, Q130, V131, K138, Q142, S196, S200, M204, G218, N219, M221, G222, L225, L228, V232, E312, I317, W318, Y341, T345, F348, G349, H353, K354, M355, G356, R357, G396, D399, S400, M403, P404, D426, F429, Y433',
     'pore2: L30, D33, N34, L37, T38, E120, N128, Q130, V131, K138, Q142, S196, S200, M204, A208, Y211, T212, D214, E215, R217, G218, N219, M221, G222, L225, L228, V232, E312, I317, W318, Y341, T345, F348, G349, H353, R357, G396, D399, S400, M403, P404, G407, Y408, D411, Y418, Y422, D426, F429, Y433',
     'pore3: L30, D33, N34, L37, T38, V40, V41, I44, P45, V131, L134, F135, K138, Q142, S180, L185, R189, S196, S200, M204, G218, N219, M221, G222, L225, L228, V232, E312, Y341, T345, F348, G349, H353, K354, M355, G356, R357, G396, D399, S400, M403, P404, D426, F429, Y433',
     'pore4: I22, I25, V26, A29, L30, D33, N34, L37, T38, V41, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M206, L228, V232, L270, Q276, E312, L315, P316, I317, Q329, A333, F334, I381, L384, N388, D426, F429, Y433',
     'pore5: R17, I22, I25, V26, A29, L30, D33, N34, L37, T38, V41, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M206, L228, V232, L270, Q276, E312, L315, P316, I317, Q329, A333, F334, I381, L384, N388, D426, F429, Y433',
     'pore6: L30, D33, N34, L37, T38, V41, K138, Q142, S196, S200, M204, A208, Y211, T212, D214, E215, R217, G218, N219, M221, G222, L225, L228, V232, E312, L315, P316, I317, Q329, A333, F334, Y341, T345, F348, G349, H353, R357, I381, L384, N388, G396, D399, S400, M403, P404, G407, Y408, D411, Y418, Y422, D426, F429, Y433',
     'pore7: I22, I25, V26, A29, L30, D33, N34, L37, T38, K122, D123, L124, E127, N128, V131, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M206, L228, V232, L270, Q276, E312, P313, L315, P316, I317, W318, D426, F429, Y433',
     'pore8: R17, I22, I25, V26, A29, L30, D33, N34, L37, T38, K122, D123, L124, E127, N128, V131, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M206, L228, V232, L270, Q276, E312, P313, L315, P316, I317, W318, D426, F429, Y433',
     'pore9: L30, D33, N34, L37, T38, K122, D123, L124, E127, N128, V131, K138, Q142, S196, S200, M204, A208, Y211, T212, D214, E215, R217, G218, N219, M221, G222, L225, L228, V232, E312, P313, L315, P316, I317, W318, Y341, T345, F348, G349, H353, R357, G396, D399, S400, M403, P404, G407, Y408, D411, Y418, Y422, D426, F429, Y433',
     'pore10: I22, I25, V26, A29, L30, D33, N34, L37, T38, K122, D123, L124, E127, N128, V131, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M206, L228, V232, L270, Q276, E312, P313, A314, L315, P316, I317, W318, D426, F429, Y433',
     'pore11: R17, I22, I25, V26, A29, L30, D33, N34, L37, T38, K122, D123, L124, E127, N128, V131, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M206, L228, V232, L270, Q276, E312, P313, A314, L315, P316, I317, W318, D426, F429, Y433'],
    ['pore0: R17, I22, I25, V26, A29, L30, D33, N34, L37, K122, L124, N128, V131, K138, A139, Q142, Y158, P159, I162, S196, S199, S200, G203, M206, L228, L270, Q276, E312, P313, A314, L315, P316, I317, W318, F334, F429, C430, Y433',
     'pore1: R17, I22, I25, V26, A29, L30, D33, N34, L37, K122, N128, V131, K138, A139, Q142, Y158, P159, I162, S196, S199, S200, G203, M206, L228, L270, Q276, E312, P313, L315, P316, I317, W318, M319, F334, F429, C430, Y433',
     'pore2: L30, D33, N34, L37, V41, P42, I43, I44, P45, K138, A139, Q142, S196, S200, M204, D214, E215, R217, G218, M221, A224, L225, L228, Y243, K248, E312, I317, W318, M319, T322, M323, F334, R357, M403, P404, G407, Y422, D426, F429, C430, Y433',
     'pore3: R17, I22, I25, V26, A29, L30, D33, N34, L37, V41, K138, A139, Q142, Y158, P159, I162, S196, S199, S200, G203, M206, L228, L270, Q276, E312, L315, P316, I317, M319, Q329, L330, A333, F334, I381, F429, C430, Y433'],
    ['pore0: R17, S18, L21, I22, I25, V26, A29, L30, D33, N34, L37, T38, K122, L124, N128, V131, K138, Q142, I162, F166, M169, S196, S199, S200, M206, L228, V232, D262, Q266, V269, L270, E312, P313, A314, L315, P316, I317, W318, F334, D426, F429',
     'pore1: R17, S18, L21, I22, I25, V26, A29, L30, D33, N34, L37, T38, K122, L124, N128, V131, K138, Q142, I162, F166, M169, S196, S199, S200, M206, L228, V232, D262, Q266, V269, L270, E312, P313, A314, L315, P316, I317, W318, F334, D426, F429',
     'pore2: R17, S18, L21, I22, I25, V26, A29, L30, D33, N34, L37, T38, K122, D123, L124, E127, N128, V131, K138, Q142, I162, F166, M169, S196, S199, S200, M206, L228, V232, D262, Q266, V269, L270, E312, P313, I317, W318, F334, D426, F429',
     'pore3: R17, S18, L21, I22, I25, V26, A29, L30, D33, N34, L37, T38, K122, D123, L124, E127, N128, V131, K138, Q142, I162, F166, M169, S196, S199, S200, M206, L228, V232, D262, Q266, V269, L270, E312, P313, I317, W318, F334, D426, F429',
     'pore4: R17, S18, L21, I22, I25, V26, A29, L30, D33, N34, L37, T38, K122, N128, V131, K138, Q142, I162, F166, M169, S196, S199, S200, M206, L228, V232, D262, Q266, V269, L270, E312, P313, L315, P316, I317, W318, M319, F334, D426, F429',
     'pore5: R17, S18, L21, I22, I25, V26, A29, L30, D33, N34, L37, T38, K122, N128, V131, K138, Q142, I162, F166, M169, S196, S199, S200, M206, L228, V232, D262, Q266, V269, L270, E312, P313, L315, P316, I317, W318, M319, F334, D426, F429',
     'pore6: R17, S18, L21, I22, I25, V26, A29, L30, D33, N34, L37, T38, E120, D121, K122, E127, N128, V131, K138, Q142, I162, F166, M169, S196, S199, S200, M206, L228, V232, D262, Q266, V269, L270, E312, P313, I317, W318, F334, D426, F429',
     'pore7: R17, S18, L21, I22, I25, V26, A29, L30, D33, N34, L37, T38, E120, D121, K122, E127, N128, V131, K138, Q142, I162, F166, M169, S196, S199, S200, M206, L228, V232, D262, Q266, V269, L270, E312, P313, I317, W318, F334, D426, F429',
     'pore8: R17, S18, L21, I22, I25, V26, A29, L30, D33, N34, L37, T38, V41, P42, N128, V131, K138, Q142, I162, F166, M169, S196, S199, S200, M206, L228, V232, D262, Q266, V269, L270, E312, I317, W318, M319, E321, T322, R326, K327, W328, Q329, L330, F334, D426, F429',
     'pore9: R17, S18, L21, I22, I25, V26, A29, L30, D33, N34, L37, T38, V41, P42, N128, V131, K138, Q142, I162, F166, M169, S196, S199, S200, M206, L228, V232, D262, Q266, V269, L270, E312, I317, W318, M319, E321, T322, R326, K327, W328, Q329, L330, F334, D426, F429',
     'pore10: R17, S18, L21, I22, I25, V26, A29, L30, D33, N34, L37, K138, Q142, I162, F166, M169, S196, S199, S200, M206, L228, V232, D262, Q266, V269, L270, I308, E312, L315, P316, I317, W328, Q329, V332, A333, F334, L335, P336, A337, S338, I381, I385, N388, D426, F429, Y433',
     'pore11: R17, S18, L21, I22, I25, V26, A29, L30, D33, N34, L37, K138, Q142, I162, F166, M169, S196, S199, S200, M206, L228, V232, D262, Q266, V269, L270, I308, E312, L315, P316, I317, W328, Q329, V332, A333, F334, L335, P336, A337, S338, I381, I385, N388, D426, F429, Y433']]


To obtain information about pore's parameters, such as volume, length, or
bottlenck, use :func:.`getPoreParametersMultipleFrames` function.

.. ipython:: python
   :verbatim:

   getPoreParametersMultipleFrames(pores)


.. parsed-literal::

   @> Frame/model: 0
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	488.65 		69.46 		0.79
   @> pore 1: 	496.62 		71.97 		0.79
   @> pore 2: 	445.76 		69.13 		0.79
   @> pore 3: 	458.33 		71.27 		0.79
   @> pore 4: 	453.73 		71.64 		0.79
   @> pore 5: 	632.03 		82.55 		0.79
   @> pore 6: 	648.0 		83.52 		0.79
   @> pore 7: 	655.97 		84.3 		0.79
   @> pore 8: 	640.0 		85.07 		0.79
   @> pore 9: 	655.97 		86.03 		0.79
   @> pore 10: 	663.94 		86.81 		0.79
   @> Frame/model: 1
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 2
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	736.53 		76.26 		0.74
   @> pore 1: 	390.95 		70.74 		0.77
   @> pore 2: 	446.88 		74.25 		0.77
   @> pore 3: 	397.35 		73.78 		0.77
   @> Frame/model: 3
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	718.41 		66.0 		0.9
   @> pore 1: 	707.57 		67.72 		0.9
   @> pore 2: 	704.15 		70.05 		0.87
   @> pore 3: 	798.15 		82.9 		0.9
   @> pore 4: 	670.63 		68.01 		0.87
   @> pore 5: 	779.85 		74.07 		0.81
   @> pore 6: 	786.62 		75.23 		0.81
   @> pore 7: 	769.01 		75.8 		0.81
   @> pore 8: 	828.28 		77.69 		0.81
   @> pore 9: 	851.12 		79.58 		0.81
   @> pore 10: 	765.59 		78.12 		0.81
   @> pore 11: 	773.95 		78.67 		0.81
   @> pore 12: 	777.69 		79.7 		0.81
   @> pore 13: 	778.14 		79.9 		0.81
   @> pore 14: 	774.42 		79.68 		0.81
   @> pore 15: 	445.12 		74.17 		0.81
   @> pore 16: 	518.23 		90.65 		0.81
   @> pore 17: 	508.04 		90.06 		0.81
   @> Frame/model: 4
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 5
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	825.73 		65.49 		0.75
   @> pore 1: 	834.72 		69.75 		0.75
   @> pore 2: 	879.72 		71.52 		0.75
   @> pore 3: 	843.12 		72.77 		0.69
   @> pore 4: 	823.63 		70.8 		0.75
   @> pore 5: 	832.63 		75.06 		0.75
   @> pore 6: 	877.62 		76.83 		0.75
   @> pore 7: 	841.03 		78.09 		0.69
   @> pore 8: 	805.46 		72.68 		0.75
   @> pore 9: 	800.69 		73.07 		0.75
   @> pore 10: 	814.46 		76.94 		0.75
   @> pore 11: 	859.45 		78.71 		0.75
   @> pore 12: 	822.86 		79.97 		0.69
   @> pore 13: 	809.3 		73.65 		0.75
   @> pore 14: 	804.53 		74.04 		0.75
   @> pore 15: 	818.29 		77.91 		0.75
   @> pore 16: 	863.29 		79.68 		0.75
   @> pore 17: 	826.69 		80.93 		0.69
   @> pore 18: 	839.04 		79.23 		0.75
   @> pore 19: 	834.27 		79.62 		0.75
   @> pore 20: 	848.04 		83.48 		0.75
   @> pore 21: 	893.03 		85.25 		0.75
   @> pore 22: 	856.43 		86.51 		0.69
   @> Frame/model: 6
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	1081.0 		76.59 		0.83
   @> pore 1: 	1094.62 	79.18 		0.83
   @> pore 2: 	1071.63 	76.13 		0.83
   @> pore 3: 	1085.25 	78.72 		0.83
   ..
   ..
   @> Frame/model: 19
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	917.75 		66.68 		0.73
   @> pore 1: 	903.64 		68.75 		0.73
   @> pore 2: 	1005.38 	79.53 		0.75
   @> pore 3: 	837.51 		71.53 		0.73
   @> Frame/model: 20
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	827.19 		72.32 		0.74
   @> pore 1: 	825.69 		72.86 		0.74
   @> pore 2: 	837.73 		73.17 		0.74
   @> pore 3: 	836.23 		73.7 		0.74
   @> pore 4: 	825.38 		74.14 		0.74
   @> pore 5: 	823.88 		74.68 		0.74
   @> pore 6: 	841.65 		75.46 		0.74
   @> pore 7: 	840.15 		76.0 		0.74
   @> pore 8: 	850.76 		82.77 		0.74
   @> pore 9: 	849.26 		83.3 		0.74
   @> pore 10: 	661.25 		77.57 		0.7
   @> pore 11: 	659.75 		78.11 		0.7


   [([69.45502768105918,
      71.96649407548894,
      69.12750263241286,
      71.26708571765487,
      71.63961744090273,
      82.55451383619868,
      83.52396585265168,
      84.29783204005349,
      85.06541245274701,
      86.03486664862177,
      86.80869225174177],
     [0.787475287601636,
      0.787475287601636,
      0.787475287601636,
      0.787475287601636,
      0.787475287601636,
      0.787475287601636,
      0.787475287601636,
      0.787475287601636,
      0.787475287601636,
      0.787475287601636,
      0.787475287601636],
     [488.64964627041616,
      496.62316240682196,
      445.7552795979281,
      458.32782550479305,
      453.7287957343339,
      632.0290597965962,
      647.9970483717335,
      655.9696422330798,
      640.0025759330023,
      655.9705645081393,
      663.9431583694858]),
    ([], [], []),
    ([76.25554318296521, 70.7396012605854, 74.24741676738773, 73.7821638072704],
     [0.7377256742368303,
      0.7665948006419143,
      0.7665948006419143,
      0.7665948006419143],
     [736.5266497510764, 390.9548792283829, 446.8818284249258, 397.351310741719]),
    ([65.99612596892734,
      67.72240950857788,
      70.05324423093188,
      82.8991908134472,
      68.00508788536808,
      74.06987315935343,
      75.23354024582183,
      75.79567659473777,
      77.69478849148587,
      79.58425479138188,
      78.1216917701227,
      78.67470988459809,
      79.69602785270617,
      79.90483119958927,
      79.67745297454553,
      74.17406970340068,
      90.65437801038513,
      90.05595173996146],
     [0.9028329092896226,
      0.9028329092896226,
      0.8662049490899701,
      0.9028329092896226,
      0.8662049490899701,
      0.8061447412549654,
      0.8061447412549654,
      0.8061447412549654,
      0.8061447412549654,
      0.8061447412549654,
      0.8061447412549654,
      0.8061447412549654,
      0.8061447412549654,
      0.8061447412549654,
      0.8061447412549654,
      0.8077845691591083,
      0.8077845691591083,
      0.8077845691591083],
     [718.4089819932893,
      707.5721907118374,
      704.1519451225709,
      798.1523499571715,
      670.6317159810806,
      779.8499325383403,
      786.6155674519079,
      769.0131412568883,
      828.2774602538221,
      851.1161517501379,
      765.5928956676216,
      773.945018254014,
      777.6915703142158,
      778.1419792119323,
      774.4223480637129,
      445.12384283504224,
      518.2288512108375,
      508.03651031422163]),
    ([], [], []),
     ..
     ..
    ([66.68372279689746, 68.75221055480887, 79.53024687396352, 71.52985176307159],
     [0.734319287054954,
      0.734319287054954,
      0.7457030900913533,
      0.734319287054954],
     [917.7467921779041,
      903.6381808951435,
      1005.3774089996299,
      837.5051936376151]),
    ([72.32025313948844,
      72.85516519386799,
      73.1672543089378,
      73.70185140026537,
      74.14132642165086,
      74.67592272082813,
      75.46495842953746,
      75.99856845811442,
      82.76661156783351,
      83.30216472820338,
      77.57230691976758,
      78.10768359887855],
     [0.7364088729090561,
      0.7364088729090561,
      0.7364088729090561,
      0.7364088729090561,
      0.7364088729090561,
      0.7364088729090561,
      0.7364088729090561,
      0.7364088729090561,
      0.7364088729090561,
      0.7364088729090561,
      0.696807047060563,
      0.696807047060563],
     [827.1881069459962,
      825.6858201324296,
      837.7315469778725,
      836.229260164306,
      825.3801480802197,
      823.877861266653,
      841.6475524369033,
      840.1452656233367,
      850.7590533242167,
      849.2567665106501,
      661.2535770445566,
      659.75129023099])]
