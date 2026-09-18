.. _cavitracer_single:

I. Detection of channels in molecular dynamics (MD) trajectory
===============================================================================

Analysis of the trajectory will be performed on a short MD trajectory
containing over 200 frames of simulation performed for the vesicular monoamine 
transporter VMAT2. This protein contains 460 residues. During the analysis, 
only protein structure will be taken into consideration.

Before analyzing the trajectory, its need to be parsed (see more details
in `Trajectory tutorial`_).

.. ipython:: python
   :verbatim:

   PDBfile = 'caseStudy2.pdb'
   DCDfile = 'caseStudy2.dcd'
   atoms = parsePDB(PDBfile)
   dcd = Trajectory(DCDfile)
   dcd.link(atoms)
   dcd.setCoords(atoms)

.. parsed-literal::

   @> 56457 atoms and 1 coordinate set(s) were parsed in 0.69s.

To analyze only protein structure for the analysis, we will select protein
structure:

.. ipython:: python
   :verbatim:

   protein = atoms.select("protein")
   dcd.setAtoms(protein)

Next, to detect channels in MD trajectories, we need to use
:func:`.scalcChannelsMultipleFrames` function. To speed up the calculations,
``max_proc`` parameter can be set up to a higher number. In the example
below, we will use 4 processors to perform the calulations. The results will
saved with the prefix ``"chls_dcd"``, and all the channels will be storage
separately (``separate=True``). 

.. ipython:: python
   :verbatim:

   channels4, surfaces4=calcChannelsMultipleFrames(protein, dcd, 
		output_path='chls_dcd', separate=True, max_proc=4)

.. parsed-literal::

    @> Frame/model: 0
    @> Frame/model: 14
    @> Frame/model: 28
    @> Frame/model: 42
    @> The atoms supplied to calcChannels contain protein atoms only.
    @> The structure carries its hydrogens (105% of what a complete protein would hold), so inner_radius=1.20 is more conservative than it needs to be: the 1.2 Å floor exists only to keep a sub-water probe out of the space that missing hydrogens leave open, and here that space is filled. A smaller probe, down to about 0.9 Å, measures the narrow connections instead of reporting them closed.
    @> The atoms supplied to calcChannels contain protein atoms only.
    @> The structure carries its hydrogens (105% of what a complete protein would hold), so inner_radius=1.20 is more conservative than it needs to be: the 1.2 Å floor exists only to keep a sub-water probe out of the space that missing hydrogens leave open, and here that space is filled. A smaller probe, down to about 0.9 Å, measures the narrow connections instead of reporting them closed.
    @> The atoms supplied to calcChannels contain protein atoms only.
    @> The structure carries its hydrogens (105% of what a complete protein would hold), so inner_radius=1.20 is more conservative than it needs to be: the 1.2 Å floor exists only to keep a sub-water probe out of the space that missing hydrogens leave open, and here that space is filled. A smaller probe, down to about 0.9 Å, measures the narrow connections instead of reporting them closed.
    @> The atoms supplied to calcChannels contain protein atoms only.
    @> The structure carries its hydrogens (105% of what a complete protein would hold), so inner_radius=1.20 is more conservative than it needs to be: the 1.2 Å floor exists only to keep a sub-water probe out of the space that missing hydrogens leave open, and here that space is filled. A smaller probe, down to about 0.9 Å, measures the narrow connections instead of reporting them closed.
    @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 Å in 0.26s.
    @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 Å in 0.26s.
    @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 Å in 0.27s.
    @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 Å in 0.27s.
    @> Delaunay tessellation of 77434 points constructed in 3.84s.
    @> Delaunay tessellation of 77434 points constructed in 3.85s.
    @> Delaunay tessellation of 77434 points constructed in 3.94s.
    @> Delaunay tessellation of 77434 points constructed in 4.00s.
    @> Surface and inner simplices filtered in 4.13s.
    @> Surface and inner simplices filtered in 4.05s.
    @> Surface and inner simplices filtered in 4.21s.
    @> Surface and inner simplices filtered in 4.17s.
    @> Cavities: 178 found, 5 deeper than min_depth=5.0 Å and searched for channels, in 0.22s.
    @> Chambers (probe 1.40 Å): 3 of the 5 searched cavities have them; the other 3 are searched whole.
    @>     cavity 0: 4 chambers, 2 of them seeded.
    @>     cavity 1: 1 chamber, none of them deep and large enough to seed; searched whole.
    @>     cavity 3: 1 chamber, seeded.
    @> 6 search sites (sp) in 0.01s: one per seeded chamber, one per cavity searched whole.
    @> Cavities: 185 found, 11 deeper than min_depth=5.0 Å and searched for channels, in 0.34s.
    @> Cavities: 175 found, 7 deeper than min_depth=5.0 Å and searched for channels, in 0.27s.
    @> Channel search (Dijkstra) over 6 search sites in 5 cavities completed in 0.06s.
    @> Found 5 channels and 1 link (a link joins a deep chamber to a shallower one and never reaches the surface).
    @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
    @>     site  start_point [Å]             void                   volume [Å³]  depth [Å]  channels  links
    @>     sp0   [-3.683, -10.297, -19.324]  cavity 1, whole                548        5.3         1      -
    @>     sp1   [15.053, -1.441, -19.567]   cavity 2, whole                300        6.0         1      -
    @>     sp2   [-13.841, -3.510, 12.012]   cavity 4, whole                224        5.8         -      -  sealed
    @>     sp3   [5.789, 8.884, 11.904]      cavity 0, chamber 1/2          115        6.8         1      -
    @>     sp4   [-7.217, 7.866, 15.518]     cavity 3, chamber 1/1           78        5.7         2      -
    @>     sp5   [6.642, 3.825, 8.825]       cavity 0, chamber 2/2           59       13.2         -      1  -> sp3
    @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
    @> The 1 site marked sealed above report neither a channel nor a link: no route out of them survived - either narrower than bottleneck=1.20 Å, or dropped as a duplicate of a shallower site's, or the void is its own mouth and has nowhere to path to. Lower bottleneck to see how the first kind connect.
    @> Saving 5 channels and 1 links to directory ., one file per object named sp<site>_chl<n> and sp<site>_lnk<n>.
    @> Channel calculation completed in 8.63s.
    @> Frame/model: 1
    @> The atoms supplied to calcChannels contain protein atoms only.
    @> The structure carries its hydrogens (105% of what a complete protein would hold), so inner_radius=1.20 is more conservative than it needs to be: the 1.2 Å floor exists only to keep a sub-water probe out of the space that missing hydrogens leave open, and here that space is filled. A smaller probe, down to about 0.9 Å, measures the narrow connections instead of reporting them closed.
    @> Chambers (probe 1.40 Å): 6 of the 7 searched cavities have them; the other 4 are searched whole.
    @>     cavity 0: 2 chambers, all seeded.
    @>     cavity 1: 2 chambers, 1 of them seeded.
    @>     cavity 2: 3 chambers, none of them deep and large enough to seed; searched whole.
    @>     cavity 3: 1 chamber, none of them deep and large enough to seed; searched whole.
    @>     cavity 4: 1 chamber, none of them deep and large enough to seed; searched whole.
    @>     cavity 5: 1 chamber, seeded.
    @> 8 search sites (sp) in 0.05s: one per seeded chamber, one per cavity searched whole.
    @> Chambers (probe 1.40 Å): 7 of the 11 searched cavities have them; the other 8 are searched whole.
    @>     cavity 0: 3 chambers, 2 of them seeded.
    @>     cavity 1: 4 chambers, 1 of them seeded.
    @>     cavity 2: 2 chambers, none of them deep and large enough to seed; searched whole.
    @>     cavity 3: 3 chambers, none of them deep and large enough to seed; searched whole.
    @>     cavity 4: 2 chambers, none of them deep and large enough to seed; searched whole.
    @>     cavity 5: 3 chambers, 1 of them seeded.
    @>     cavity 6: 4 chambers, none of them deep and large enough to seed; searched whole.
    @> 12 search sites (sp) in 0.05s: one per seeded chamber, one per cavity searched whole.
    @> Channel search (Dijkstra) over 8 search sites in 7 cavities completed in 0.16s.
    @> Found 7 channels and 1 link (a link joins a deep chamber to a shallower one and never reaches the surface).
    @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
    @>     site  start_point [Å]             void                   volume [Å³]  depth [Å]  channels  links
    @>     sp0   [2.849, 3.805, 1.723]       cavity 0, chamber 1/2         1452       18.7         -      1  -> sp4
    @>     sp1   [-13.991, -3.043, 13.028]   cavity 2, whole                429        5.6         1      -
    @>     sp2   [-3.539, -10.671, -18.668]  cavity 3, whole                336        5.6         1      -
    @>     sp3   [14.449, -1.190, -18.655]   cavity 4, whole                251        5.2         1      -
    @>     sp4   [4.640, 9.767, 13.115]      cavity 0, chamber 2/2          251        5.8         1      -
    @>     sp5   [10.134, -1.611, 16.422]    cavity 6, whole                133        5.3         1      -
    @>     sp6   [13.306, 9.826, 14.612]     cavity 5, chamber 1/1           76       10.5         -      -  sealed
    @>     sp7   [-8.841, 7.348, 17.300]     cavity 1, chamber 1/1           58        6.0         2      -
    @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
    @> The 1 site marked sealed above report neither a channel nor a link: no route out of them survived - either narrower than bottleneck=1.20 Å, or dropped as a duplicate of a shallower site's, or the void is its own mouth and has nowhere to path to. Lower bottleneck to see how the first kind connect.
    @> Saving 7 channels and 1 links to directory ., one file per object named sp<site>_chl<n> and sp<site>_lnk<n>.
    @> Channel calculation completed in 8.77s.
    @> Cavities: 195 found, 5 deeper than min_depth=5.0 Å and searched for channels, in 0.36s.
    ..
    ..
   @> Frame/model: 209
    @> The atoms supplied to calcChannels contain protein atoms only.
    @> The structure carries its hydrogens (105% of what a complete protein would hold), so inner_radius=1.20 is more conservative than it needs to be: the 1.2 Å floor exists only to keep a sub-water probe out of the space that missing hydrogens leave open, and here that space is filled. A smaller probe, down to about 0.9 Å, measures the narrow connections instead of reporting them closed.
    @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 Å in 0.29s.
    @> Surface and inner simplices filtered in 3.06s.
    @> Cavities: 165 found, 6 deeper than min_depth=5.0 Å and searched for channels, in 0.35s.
    @> Delaunay tessellation of 77434 points constructed in 3.11s.
    @> Chambers (probe 1.40 Å): 6 of the 6 searched cavities have them; the other 4 are searched whole.
    @>     cavity 0: 10 chambers, 1 of them seeded.
    @>     cavity 1: 3 chambers, 1 of them seeded.
    @>     cavity 2: 2 chambers, none of them deep and large enough to seed; searched whole.
    @>     cavity 3: 2 chambers, none of them deep and large enough to seed; searched whole.
    @>     cavity 4: 2 chambers, none of them deep and large enough to seed; searched whole.
    @>     cavity 5: 1 chamber, none of them deep and large enough to seed; searched whole.
    @> 6 search sites (sp) in 0.08s: one per seeded chamber, one per cavity searched whole.
    @> Channel search (Dijkstra) over 6 search sites in 6 cavities completed in 0.19s.
    @> Found 8 channels.
    @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
    @>     site  start_point [Å]             void                   volume [Å³]  depth [Å]  channels  links
    @>     sp0   [0.516, -0.662, 11.707]     cavity 0, chamber 1/1         3079       11.1         3      -
    @>     sp1   [-3.599, -10.845, -19.331]  cavity 2, whole                488        5.4         2      -
    @>     sp2   [12.408, 6.748, 12.005]     cavity 3, whole                282        5.0         1      -
    @>     sp3   [16.923, 1.103, -8.312]     cavity 4, whole                145        9.4         -      -  sealed
    @>     sp4   [15.992, -2.265, 3.221]     cavity 5, whole                133        5.1         1      -
    @>     sp5   [10.351, 8.329, -13.347]    cavity 1, chamber 1/1          107        5.6         1      -
    @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
    @> The 1 site marked sealed above report neither a channel nor a link: no route out of them survived - either narrower than bottleneck=1.20 Å, or dropped as a duplicate of a shallower site's, or the void is its own mouth and has nowhere to path to. Lower bottleneck to see how the first kind connect.
    @> Saving 8 channels to directory ., one file per object named sp<site>_chl<n>.
    @> Channel calculation completed in 7.21s.
    @> Delaunay tessellation of 77434 points constructed in 3.11s.
    @> Surface and inner simplices filtered in 3.27s.
    @> Cavities: 143 found, 7 deeper than min_depth=5.0 Å and searched for channels, in 0.29s.
    @> Chambers (probe 1.40 Å): 6 of the 7 searched cavities have them; the other 3 are searched whole.
    @>     cavity 0: 9 chambers, 2 of them seeded.
    @>     cavity 1: 2 chambers, 1 of them seeded.
    @>     cavity 2: 1 chamber, seeded.
    @>     cavity 3: 2 chambers, 1 of them seeded.
    @>     cavity 4: 1 chamber, none of them deep and large enough to seed; searched whole.
    @>     cavity 5: 1 chamber, none of them deep and large enough to seed; searched whole.
    @> 8 search sites (sp) in 0.07s: one per seeded chamber, one per cavity searched whole.
    @> Channel search (Dijkstra) over 8 search sites in 7 cavities completed in 0.23s.
    @> Found 11 channels.
    @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
    @>     site  start_point [Å]             void                   volume [Å³]  depth [Å]  channels  links
    @>     sp0   [-0.638, 0.102, 11.113]     cavity 0, chamber 1/2         2384       11.8         3      -
    @>     sp1   [-4.464, 5.343, 18.136]     cavity 2, chamber 1/1          316       12.6         1      -
    @>     sp2   [-11.826, -3.055, 12.025]   cavity 0, chamber 2/2          245        5.0         1      -
    @>     sp3   [9.101, 1.844, -18.305]     cavity 4, whole                229        5.1         1      -
    @>     sp4   [7.010, 7.107, -13.851]     cavity 3, chamber 1/1          167        7.0         1      -
    @>     sp5   [-14.717, 3.685, -10.866]   cavity 5, whole                165        5.0         -      -  sealed
    @>     sp6   [-11.989, 10.577, 15.484]   cavity 6, whole                127        5.1         1      -
    @>     sp7   [-2.816, -11.821, -18.928]  cavity 1, chamber 1/1           68        7.1         3      -
    @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
    @> The 1 site marked sealed above report neither a channel nor a link: no route out of them survived - either narrower than bottleneck=1.20 Å, or dropped as a duplicate of a shallower site's, or the void is its own mouth and has nowhere to path to. Lower bottleneck to see how the first kind connect.
    @> Saving 11 channels to directory ., one file per object named sp<site>_chl<n>.
    @> Channel calculation completed in 7.20s.
    @> Surface and inner simplices filtered in 3.22s.
    @> Cavities: 183 found, 4 deeper than min_depth=5.0 Å and searched for channels, in 0.32s.
    @> Chambers (probe 1.40 Å): 3 of the 4 searched cavities have them; the other 2 are searched whole.
    @>     cavity 0: 8 chambers, 4 of them seeded.
    @>     cavity 1: 1 chamber, none of them deep and large enough to seed; searched whole.
    @>     cavity 2: 1 chamber, seeded.
    @> 7 search sites (sp) in 0.06s: one per seeded chamber, one per cavity searched whole.
    @> Channel search (Dijkstra) over 7 search sites in 4 cavities completed in 0.22s.
    @> Found 14 channels and 2 links (a link joins a deep chamber to a shallower one and never reaches the surface).
    @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
    @>     site  start_point [Å]             void                   volume [Å³]  depth [Å]  channels  links
    @>     sp0   [-0.147, -1.002, 11.590]    cavity 0, chamber 1/4         2548       10.2         6      1  -> sp2
    @>     sp1   [-3.179, -12.027, -18.873]  cavity 1, whole                898        5.5         1      -
    @>     sp2   [-9.973, -4.802, 11.669]    cavity 0, chamber 2/4          250        5.1         3      -
    @>     sp3   [-5.381, 5.295, 18.108]     cavity 0, chamber 3/4          203       11.0         -      1  -> sp0
    @>     sp4   [15.187, -2.137, -16.901]   cavity 3, whole                189        5.2         1      -
    @>     sp5   [7.630, 8.967, -16.557]     cavity 2, chamber 1/1          142        5.1         3      -
    @>     sp6   [-6.792, -3.004, -5.696]    cavity 0, chamber 4/4           60       33.4         -      -  sealed
    @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
    @> The 1 site marked sealed above report neither a channel nor a link: no route out of them survived - either narrower than bottleneck=1.20 Å, or dropped as a duplicate of a shallower site's, or the void is its own mouth and has nowhere to path to. Lower bottleneck to see how the first kind connect.
    @> Saving 14 channels and 2 links to directory ., one file per object named sp<site>_chl<n> and sp<site>_lnk<n>.
    @> Channel calculation completed in 7.21s.


All the details about the predicted channels can be displayed using
:func:`.getChannelParametersMultipleFrames`.

.. ipython:: python
   :verbatim:

   getChannelParametersMultipleFrames(channels4, param_file_name='DATA_chls_dcd')

.. parsed-literal::

    @> Frame/model: 0
    @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
    @> channel 0: 	133.94 		5.33 		1.57
    @> channel 1: 	122.67 		7.0 		1.94
    @> channel 2: 	109.21 		6.09 		1.23
    @> channel 3: 	60.68 		5.86 		1.2
    @> channel 4: 	72.68 		7.99 		1.25
    @> Frame/model: 1
    @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
    @> channel 0: 	145.47 		6.67 		2.09
    @> channel 1: 	127.23 		5.42 		1.44
    @> channel 2: 	131.92 		6.35 		1.5
    @> channel 3: 	61.92 		5.78 		1.32
    @> channel 4: 	42.95 		5.4 		1.21
    @> Frame/model: 2
    @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
    @> channel 0: 	118.38 		5.35 		1.37
    @> channel 1: 	107.01 		5.29 		1.36
    @> channel 2: 	129.91 		8.38 		1.77
    @> channel 3: 	73.14 		6.6 		1.32
    @> Frame/model: 3
    @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
    @> channel 0: 	137.98 		5.5 		1.51
    @> channel 1: 	108.45 		5.69 		1.91
    @> channel 2: 	121.44 		5.92 		1.43
    @> channel 3: 	66.54 		5.06 		1.5
    @> channel 4: 	59.88 		5.35 		1.27
    @> channel 5: 	65.03 		7.41 		1.34
    @> Frame/model: 4
    @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
    @> channel 0: 	102.95 		5.53 		1.78
    @> channel 1: 	129.82 		6.09 		1.42
    @> channel 2: 	114.55 		5.88 		1.3
    @> channel 3: 	54.64 		5.79 		1.21
    @> Frame/model: 5
    @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
    @> channel 0: 	126.05 		5.84 		1.46
    @> channel 1: 	109.5 		5.47 		1.32
    @> channel 2: 	90.93 		5.07 		1.3
    @> channel 3: 	135.77 		9.06 		1.8
    @> channel 4: 	63.41 		6.1 		1.28
    @> Frame/model: 6
    @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
    @> channel 0: 	105.43 		5.44 		1.5
    @> channel 1: 	142.46 		7.89 		1.91
    @> channel 2: 	100.91 		6.2 		1.39
    @> channel 3: 	86.71 		5.67 		1.21
    @> channel 4: 	74.72 		7.45 		1.44
    @> channel 5: 	81.93 		7.16 		1.27
    @> channel 6: 	93.11 		7.81 		1.25
    @> channel 7: 	101.84 		8.9 		1.35
    @> channel 8: 	56.39 		7.73 		1.23
    @> Frame/model: 7
    @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
    @> channel 0: 	111.33 		5.2 		1.52
    @> channel 1: 	110.31 		6.27 		1.74
    @> channel 2: 	121.22 		5.52 		1.37
    @> channel 3: 	102.25 		5.61 		1.49
    @> channel 4: 	69.79 		5.05 		1.21
    @> Frame/model: 8
    @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
    @> channel 0: 	103.97 		5.31 		1.55
    @> channel 1: 	122.22 		6.84 		1.91
    @> channel 2: 	82.63 		5.34 		1.38
    @> channel 3: 	55.86 		5.83 		1.27
    @> channel 4: 	53.21 		6.15 		1.3
    @> Frame/model: 9
    @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
    @> channel 0: 	82.2 		5.17 		1.46
    @> channel 1: 	123.87 		7.87 		1.63
    @> channel 2: 	51.53 		5.11 		1.24
    @> channel 3: 	170.7 		15.21 		1.28
    @> Frame/model: 10
    @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
    @> channel 0: 	94.01 		5.17 		1.67
    @> channel 1: 	91.23 		5.57 		1.45
    @> channel 2: 	58.61 		5.04 		1.31
    @> channel 3: 	73.58 		5.35 		1.49
    @> channel 4: 	48.58 		5.14 		1.29
    @> channel 5: 	55.45 		6.07 		1.22
    ..
    ..
    @> Frame/model: 210
    @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
    @> channel 0: 	167.15 		5.05 		2.01
    @> channel 1: 	83.6 		6.19 		1.33
    @> channel 2: 	67.07 		5.3 		1.28
    @> channel 3: 	289.49 		11.26 		1.66
    @> channel 4: 	308.26 		13.46 		1.68
    @> channel 5: 	252.38 		11.87 		1.47
    @> channel 6: 	314.15 		13.51 		1.3
    @> channel 7: 	343.86 		16.91 		1.34
    @> channel 8: 	324.33 		17.24 		1.37
    @> channel 9: 	363.67 		17.5 		1.31
    [([5.331767742552358,
       6.998387387741667,
       6.085527676462607,
       5.858854311444432,
       7.98983939686644],
      [1.572155685221137,
       1.9405313948234382,
       1.2296813140535539,
       1.20167477329606,
       1.2473063269041071],
      [133.9415516958124,
       122.67418882185356,
       109.20952421817532,
       60.68230076564663,
       72.6753897409704]),
     ([6.666208430331574,
       5.424479193692954,
       6.346111908810009,
       5.777844656631124,
       5.402867928697702],
      [2.090812708434658,
       1.4402877557423575,
       1.504851521494053,
       1.3243104022063263,
       1.2119138481263352],
      [145.47324696045195,
       127.22547601324104,
       131.92497925196335,
       61.91661146109752,
       42.95162370803827]),
     ..
     ..
     ([5.053279294414388,
       6.189079374115059,
       5.304283702411966,
       11.262593963001393,
       13.459036504111774,
       11.871591748508564,
       13.509516942638102,
       16.913899678851095,
       17.238638340624938,
       17.49961773787033],
      [2.010772399129248,
       1.3253521540721842,
       1.2797909473647724,
       1.6639913436549463,
       1.684645260508318,
       1.4660453547846706,
       1.2963158417703402,
       1.3370396034799514,
       1.3738494561152466,
       1.3114666587307637],
      [167.1493593307516,
       83.60335821671228,
       67.07295039682602,
       289.4915096364011,
       308.2623126013169,
       252.382705629095,
       314.1454320762822,
       343.8636015359858,
       324.33082262735957,
       363.66523396339335])]



Assigning the function's output to the ``results`` variable grants access 
to these parameters, as shown below.

.. ipython:: python
   :verbatim:

   results = getChannelParametersMultipleFrames(channels4)

.. parsed-literal::

    @> Frame/model: 0
    @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
    @> channel 0: 	133.94 		5.33 		1.57
    @> channel 1: 	122.67 		7.0 		1.94
    @> channel 2: 	109.21 		6.09 		1.23
    @> channel 3: 	60.68 		5.86 		1.2
    @> channel 4: 	72.68 		7.99 		1.25
    @> Frame/model: 1
    @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
    @> channel 0: 	145.47 		6.67 		2.09
    @> channel 1: 	127.23 		5.42 		1.44
    @> channel 2: 	131.92 		6.35 		1.5
    @> channel 3: 	61.92 		5.78 		1.32
    @> channel 4: 	42.95 		5.4 		1.21
    @> Frame/model: 2
    @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
    @> channel 0: 	118.38 		5.35 		1.37
    @> channel 1: 	107.01 		5.29 		1.36
    @> channel 2: 	129.91 		8.38 		1.77
    @> channel 3: 	73.14 		6.6 		1.32
    @> Frame/model: 3
    @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
    @> channel 0: 	137.98 		5.5 		1.51
    @> channel 1: 	108.45 		5.69 		1.91
    @> channel 2: 	121.44 		5.92 		1.43
    @> channel 3: 	66.54 		5.06 		1.5
    @> channel 4: 	59.88 		5.35 		1.27
    @> channel 5: 	65.03 		7.41 		1.34
    @> Frame/model: 4
    @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
    @> channel 0: 	102.95 		5.53 		1.78
    @> channel 1: 	129.82 		6.09 		1.42
    @> channel 2: 	114.55 		5.88 		1.3
    @> channel 3: 	54.64 		5.79 		1.21
    @> Frame/model: 5
    @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
    @> channel 0: 	126.05 		5.84 		1.46
    @> channel 1: 	109.5 		5.47 		1.32
    @> channel 2: 	90.93 		5.07 		1.3
    @> channel 3: 	135.77 		9.06 		1.8
    @> channel 4: 	63.41 		6.1 		1.28
    ..
    ..
    @> Frame/model: 210
    @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
    @> channel 0: 	167.15 		5.05 		2.01
    @> channel 1: 	83.6 		6.19 		1.33
    @> channel 2: 	67.07 		5.3 		1.28
    @> channel 3: 	289.49 		11.26 		1.66
    @> channel 4: 	308.26 		13.46 		1.68
    @> channel 5: 	252.38 		11.87 		1.47
    @> channel 6: 	314.15 		13.51 		1.3
    @> channel 7: 	343.86 		16.91 		1.34
    @> channel 8: 	324.33 		17.24 		1.37
    @> channel 9: 	363.67 		17.5 		1.31


Results for the first frame in the trajectory:

.. ipython:: python
   :verbatim:

   results[0]

.. parsed-literal::

    ([5.331767742552358,
      6.998387387741667,
      6.085527676462607,
      5.858854311444432,
      7.98983939686644],
     [1.572155685221137,
      1.9405313948234382,
      1.2296813140535539,
      1.20167477329606,
      1.2473063269041071],
     [133.9415516958124,
      122.67418882185356,
      109.20952421817532,
      60.68230076564663,
      72.6753897409704])

To obtain information about the lengths of the channels detected in the 
first frame in the trajectory (#0):

.. ipython:: python
   :verbatim:

   results[0][0]

.. parsed-literal::

    [5.331767742552358,
     6.998387387741667,
     6.085527676462607,
     5.858854311444432,
     7.98983939686644]

Bottlenecks of the channels in the first frame in the trajectory (#0):

.. ipython:: python
   :verbatim:

   results[0][1]

.. parsed-literal::

    [1.572155685221137,
     1.9405313948234382,
     1.2296813140535539,
     1.20167477329606,
     1.2473063269041071]

Volume of the channels detected in the first frame in the trajectory:

.. ipython:: python
   :verbatim:

   results[0][2]

.. parsed-literal::

    [133.9415516958124,
     122.67418882185356,
     109.20952421817532,
     60.68230076564663,
     72.6753897409704]

Once we have access to ``results``, we can display bottleneck data for all
the channels in the following way:

.. ipython:: python
   :verbatim:

   import matplotlib.pylab as plt
   all_Bottleneck = []
   for nr_i,i in enumerate(results):
      all_Bottleneck.extend(results[nr_i][1])

   plt.hist(all_Bottleneck)
   plt.show()


.. figure:: images/cavitracer_figure16.jpg
   :scale: 50 %

To obtain information about all residues that are participating in the
formation of channels, use :func:`.getChannelResidueNamesMultipleFrames`.
It may take some time because each channel is mapped to a particular
protein structure. Therefore, if only a paricualar region is interested for
us, we first might apply filtering of the channels which is described below.

.. ipython:: python
   :verbatim:

   residuesALL = getChannelResidueNamesMultipleFrames(protein, 
					channels4, dcd, 
					residues_file_name='DCD_res_ALL')

.. parsed-literal::

   @> Frame: 0
   @> Channel residues were saved to: DCD_res_ALL_frame0_Residues_All_channels.txt
   @> Frame: 1
   @> Channel residues were saved to: DCD_res_ALL_frame1_Residues_All_channels.txt
   @> Frame: 2
   @> Channel residues were saved to: DCD_res_ALL_frame2_Residues_All_channels.txt
   @> Frame: 3
   @> Channel residues were saved to: DCD_res_ALL_frame3_Residues_All_channels.txt
   @> Frame: 4
   @> Channel residues were saved to: DCD_res_ALL_frame4_Residues_All_channels.txt
   @> Frame: 5
   @> Channel residues were saved to: DCD_res_ALL_frame5_Residues_All_channels.txt
   @> Frame: 6
   @> Channel residues were saved to: DCD_res_ALL_frame6_Residues_All_channels.txt
   @> Frame: 7
   @> Channel residues were saved to: DCD_res_ALL_frame7_Residues_All_channels.txt
   @> Frame: 8
   @> Channel residues were saved to: DCD_res_ALL_frame8_Residues_All_channels.txt
   @> Frame: 9
   @> Channel residues were saved to: DCD_res_ALL_frame9_Residues_All_channels.txt
   @> Frame: 10
   @> Channel residues were saved to: DCD_res_ALL_frame10_Residues_All_channels.txt
   ..
   ..
   @> Frame: 208
   @> Channel residues were saved to: DCD_res_ALL_frame208_Residues_All_channels.txt
   @> Frame: 209
   @> Channel residues were saved to: DCD_res_ALL_frame209_Residues_All_channels.txt
   @> Frame: 210
   @> Channel residues were saved to: DCD_res_ALL_frame210_Residues_All_channels.txt


To have access to a particular frame, we need to use :meth:`.getFrame`.
Below, we will select third frame from the simulation (counting from 0):

.. ipython:: python
   :verbatim:

   frame3 = dcd.getFrame(2)
   frame3

.. parsed-literal::

   <Frame: 2 from caseStudy2 (selected 5986 of 56457 atoms)>


To obtain infromation about residues for a particular frame, the following
operations should be performed:

.. ipython:: python
   :verbatim:

   protein_frame3 = protein.copy()
   protein_frame3.setCoords(frame3.getCoords())

.. ipython:: python
   :verbatim:

   getChannelResidueNames(protein_frame3, channels4[2], 
				residues_file_name='DCD_fr3_res')

.. parsed-literal::

    @> Channel residues were saved to: DCD_fr3_res_Residues_All_channels.txt

    ['channel0: THR153:P, ASN154:P, GLY157:P, TYR158:P, SER209:P, GLN276:P, PRO277:P, GLU278:P, SER279:P, VAL417:P, TYR418:P, GLY419:P',
     'channel1: ASP214:P, HSP353:P, GLY356:P, ARG357:P, TRP358:P, ILE405:P, TYR408:P, LEU409:P, LEU470:P, ARG471:P, PRO473:P',
     'channel2: LEU315:P, ILE317:P, GLN329:P, ALA333:P, PRO336:P, ILE381:P, TYR382:P, LEU384:P, ILE385:P, ASN388:P',
     'channel3: PRO42:P, PRO45:P, SER240:P, TYR243:P, GLU244:P, LYS327:P, TRP328:P, LEU330:P, GLY331:P']


To display the results for a particular frame, we need to create a model
using :func:`.getVmdModel` function and provide a path to VMD_.

.. ipython:: python
   :verbatim:

   vmd_path = '/usr/local/bin/vmd'
   model3_traj = getVmdModel(vmd_path, frame3)

.. parsed-literal::

   @> Model created successfully.

Once the model is created, the channels can be displayed together with the
protein structure using :func:`.showChannels`.

.. ipython:: python
   :verbatim:

   showChannels(channels4[2], model=model3_traj)

.. figure:: images/cavitracer_figure15.jpg
   :scale: 50 %

To select channels with particular localization, we can further use
:func:`.selectChannelBySelection` function and specify which region will be
analyzed, for example, using residues, as shown below. In the example, we
are selecting channels that are generated from MD trajectory using
``"chls_dcd*_chl*.pqr"`` pattern. List with ``pqr`` files
called ``pqr_files_channels`` is further use by
:func:`.selectChannelBySelection`.

.. ipython:: python
   :verbatim:
   
   from pathlib import Path
   pqr_files_channels = [i.name for i in Path(".").glob("chls_dcd*_chl*.pqr") if i.is_file()]
   pqr_files_channels

.. parsed-literal::

    ['chls_dcd168_sp6_chl4.pqr',
     'chls_dcd0_sp1_chl2.pqr',
     'chls_dcd71_sp0_chl8.pqr',
     'chls_dcd132_sp0_chl6.pqr',
     'chls_dcd49_sp0_chl12.pqr',
     'chls_dcd44_sp4_chl11.pqr',
     'chls_dcd157_sp0_chl6.pqr',
     'chls_dcd107_sp0_chl10.pqr',
     'chls_dcd190_sp0_chl6.pqr',
     'chls_dcd154_sp1_chl0.pqr',
     'chls_dcd20_sp10_chl7.pqr',
     'chls_dcd127_sp1_chl9.pqr',
     'chls_dcd20_sp8_chl0.pqr',
     'chls_dcd34_sp1_chl4.pqr',
     'chls_dcd57_sp0_chl11.pqr',
     'chls_dcd80_sp0_chl7.pqr',
     'chls_dcd204_sp0_chl6.pqr',
     'chls_dcd139_sp0_chl4.pqr',
     'chls_dcd57_sp6_chl4.pqr',
     'chls_dcd162_sp4_chl4.pqr',
     'chls_dcd132_sp0_chl7.pqr',
     ..
     ..
     'chls_dcd174_sp0_chl9.pqr',
     'chls_dcd78_sp4_chl2.pqr',
     'chls_dcd185_sp2_chl0.pqr',
     'chls_dcd28_sp3_chl2.pqr',
     'chls_dcd133_sp0_chl13.pqr',
     'chls_dcd70_sp1_chl7.pqr',
     'chls_dcd102_sp1_chl4.pqr',
     'chls_dcd200_sp7_chl5.pqr',
     'chls_dcd164_sp2_chl5.pqr',
     'chls_dcd151_sp3_chl0.pqr',
     'chls_dcd96_sp0_chl6.pqr',
     'chls_dcd63_sp1_chl0.pqr',
     'chls_dcd69_sp3_chl3.pqr',
     'chls_dcd103_sp5_chl7.pqr',
     'chls_dcd12_sp5_chl0.pqr',
     'chls_dcd24_sp0_chl8.pqr',
     'chls_dcd200_sp4_chl4.pqr',
     'chls_dcd43_sp3_chl3.pqr',
     'chls_dcd90_sp4_chl2.pqr',
     ...]


.. ipython:: python
   :verbatim:

   selectChannelBySelection(atoms, pqr_files=pqr_files_channels, 
			residue_sele='resid 135 37 and backbone', 
                        folder_name="Selected_channel1", 
			distA=4.0)

.. parsed-literal::

    @> 95 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 65 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 170 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 95 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 200 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 110 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 125 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 160 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 150 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 50 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 65 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 95 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 55 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 55 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 235 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 80 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 170 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 80 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 60 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 65 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 100 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 40 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 210 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 40 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 40 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 50 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 85 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 40 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 125 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 115 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 195 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 70 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 45 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 80 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 140 atoms and 1 coordinate sets were parsed in 0.00s.
    ..
    ..
    @> 40 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 185 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 100 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 40 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 60 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 35 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 85 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 50 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 55 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 65 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 85 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 70 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 75 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 120 atoms and 1 coordinate sets were parsed in 0.00s.
    @> Filtered files are now in: Selected_channel1
    @> Selected files: 
    @> chls_dcd71_sp0_chl8.pqr chls_dcd132_sp0_chl7.pqr chls_dcd100_sp0_chl9.pqr chls_dcd92_sp0_chl7.pqr chls_dcd89_sp0_chl9.pqr chls_dcd141_sp0_chl4.pqr chls_dcd175_sp0_chl4.pqr chls_dcd93_sp0_chl7.pqr chls_dcd188_sp0_chl8.pqr chls_dcd145_sp0_chl2.pqr chls_dcd180_sp0_chl10.pqr chls_dcd165_sp0_chl10.pqr chls_dcd131_sp0_chl11.pqr chls_dcd187_sp0_chl6.pqr chls_dcd133_sp0_chl6.pqr chls_dcd75_sp0_chl9.pqr chls_dcd65_sp0_chl3.pqr chls_dcd74_sp0_chl5.pqr chls_dcd191_sp1_chl8.pqr chls_dcd68_sp0_chl5.pqr chls_dcd35_sp1_chl5.pqr chls_dcd149_sp0_chl2.pqr chls_dcd205_sp0_chl3.pqr chls_dcd86_sp0_chl7.pqr chls_dcd173_sp0_chl10.pqr chls_dcd120_sp0_chl6.pqr chls_dcd98_sp0_chl5.pqr chls_dcd177_sp0_chl4.pqr chls_dcd66_sp0_chl7.pqr chls_dcd191_sp1_chl11.pqr chls_dcd181_sp0_chl3.pqr chls_dcd114_sp0_chl9.pqr chls_dcd123_sp0_chl12.pqr chls_dcd67_sp0_chl8.pqr chls_dcd93_sp0_chl6.pqr chls_dcd59_sp0_chl5.pqr chls_dcd153_sp0_chl3.pqr chls_dcd210_sp0_chl5.pqr chls_dcd184_sp0_chl6.pqr chls_dcd107_sp0_chl9.pqr chls_dcd179_sp1_chl8.pqr chls_dcd111_sp0_chl3.pqr chls_dcd59_sp0_chl4.pqr chls_dcd206_sp0_chl4.pqr chls_dcd90_sp0_chl10.pqr chls_dcd102_sp0_chl11.pqr chls_dcd88_sp0_chl10.pqr chls_dcd113_sp0_chl6.pqr chls_dcd72_sp0_chl6.pqr chls_dcd91_sp0_chl9.pqr chls_dcd171_sp0_chl7.pqr chls_dcd66_sp0_chl8.pqr chls_dcd196_sp0_chl3.pqr chls_dcd143_sp2_chl5.pqr chls_dcd152_sp0_chl2.pqr chls_dcd114_sp0_chl8.pqr chls_dcd181_sp0_chl7.pqr chls_dcd147_sp0_chl7.pqr chls_dcd95_sp0_chl6.pqr chls_dcd101_sp0_chl5.pqr chls_dcd155_sp0_chl3.pqr chls_dcd165_sp0_chl9.pqr chls_dcd140_sp0_chl3.pqr chls_dcd135_sp0_chl2.pqr chls_dcd76_sp0_chl11.pqr chls_dcd60_sp0_chl5.pqr chls_dcd101_sp0_chl6.pqr chls_dcd209_sp0_chl4.pqr chls_dcd191_sp0_chl9.pqr chls_dcd124_sp11_chl10.pqr chls_dcd97_sp0_chl9.pqr chls_dcd73_sp0_chl10.pqr chls_dcd178_sp0_chl5.pqr chls_dcd116_sp0_chl6.pqr chls_dcd106_sp0_chl13.pqr chls_dcd88_sp0_chl11.pqr chls_dcd176_sp0_chl10.pqr chls_dcd149_sp0_chl7.pqr chls_dcd110_sp0_chl6.pqr chls_dcd105_sp0_chl7.pqr chls_dcd137_sp0_chl11.pqr chls_dcd122_sp0_chl7.pqr chls_dcd117_sp0_chl7.pqr chls_dcd96_sp0_chl7.pqr chls_dcd74_sp0_chl9.pqr chls_dcd112_sp0_chl6.pqr


Filtered pqr files will be stored in ``"Selected_channel1"``. Once the files
are filtered, the :func:`.calcChannelSurfaceOverlaps` function can be used to
display the overlapping surface that is shared by filtered channels, which
will be saved when using ``output_file_name`` option.

.. ipython:: python
   :verbatim:

   calcChannelSurfaceOverlaps(pqr_files="./Selected_channel1", 
			output_file_name='overlapping_surf_traj.pdb',
			max_proc=4)

.. parsed-literal::

    @> Number of PQR files: 86
    @> Resolution: 0.5
    @> max_proc: 4
    @> Calculating overlaps using 4 processes.
    @> 85 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 105 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 170 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 255 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 75 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 240 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 140 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 210 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 100 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 110 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 100 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 105 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 170 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 75 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 115 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 115 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 160 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 220 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 165 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 145 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 80 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 180 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 125 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 105 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 80 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 210 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 270 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 210 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 135 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 260 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 210 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 160 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 130 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 150 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 195 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 100 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 120 atoms and 1 coordinate sets were parsed in 0.00s.
    ..
    ..
    @> 260 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 170 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 275 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 110 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 125 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 110 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 90 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 170 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 220 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 175 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 180 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 225 atoms and 1 coordinate sets were parsed in 0.00s.
    @> Overlap written to: overlapping_surf_traj.pdb
    @> Number of occupied overlap voxels: 24776

    'overlapping_surf_traj.pdb'

.. figure:: images/cavitracer_figure17.jpg
   :scale: 50 %

.. figure:: images/cavitracer_figure18.jpg
   :scale: 50 %


II. Reconstruction of pores in molecular dynamics (MD) trajectory
===============================================================================

Now, we will use the same MD trajectory to show how to detect pores with
protein structure.

First, we need to upload the data in the same way as for the channels
calculations.

.. ipython:: python
   :verbatim:

   PDBfile = 'caseStudy2.pdb'
   DCDfile = 'caseStudy2.dcd'
   atoms = parsePDB(PDBfile)
   dcd = Trajectory(DCDfile)
   dcd.link(atoms)
   dcd.setCoords(atoms)
   protein = atoms.select("protein")
   dcd.setAtoms(protein)

.. parsed-literal::

   @> 56457 atoms and 1 coordinate set(s) were parsed in 0.77s.

Next, we use :func:`.calcChannelsMultipleFrames`, but this time with
``return_details=True``. Without it, it is not possible to obtain
information about pores. Additionally, this time, ``inner_radius=0.8`` 
to find narrower passages within protein structure.

.. ipython:: python
   :verbatim:

   channels, surface, details = calcChannelsMultipleFrames(protein, dcd,
				inner_radius=0.8, 
                                output_path='ch_dcd_', separate=True, 
				return_details=True, max_proc=4)

.. parsed-literal::

   @> Frame/model: 0
   @> Frame/model: 14
   @> Frame/model: 28
   @> Frame/model: 42
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 Å in 0.28s.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 Å in 0.29s.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 Å in 0.29s.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 Å in 0.29s.
   @> Delaunay tessellation of 77434 points constructed in 3.97s.
   @> Delaunay tessellation of 77434 points constructed in 3.99s.
   @> Delaunay tessellation of 77434 points constructed in 4.02s.
   @> Delaunay tessellation of 77434 points constructed in 4.11s.
   @> Surface and inner simplices filtered in 3.64s.
   @> Surface and inner simplices filtered in 3.77s.
   @> Surface and inner simplices filtered in 3.83s.
   @> Surface and inner simplices filtered in 3.86s.
   @> Surface cavities: 395 found, 10 deeper than min_depth=5.0 Å and searched for channels, in 0.66s.
   @> Surface cavities: 431 found, 7 deeper than min_depth=5.0 Å and searched for channels, in 0.73s.
   @> Surface cavities: 409 found, 6 deeper than min_depth=5.0 Å and searched for channels, in 0.69s.
   @> Surface cavities: 431 found, 14 deeper than min_depth=5.0 Å and searched for channels, in 0.71s.
   @> Chambers (probe 1.40 Å): 1 of the 10 searched cavities have them; the other 9 are searched whole.
   @>     cavity 0: 21 chambers, 8 of them seeded.
   @> 17 search sites (sp) in 0.05s: one per seeded chamber, one per cavity searched whole.
   @> Chambers (probe 1.40 Å): 3 of the 14 searched cavities have them; the other 12 are searched whole.
   @>     cavity 0: 9 chambers, 5 of them seeded.
   @>     cavity 1: 5 chambers, 3 of them seeded.
   @>     cavity 2: 1 chamber, none of them deep and large enough to seed; searched whole.
   @> 20 search sites (sp) in 0.06s: one per seeded chamber, one per cavity searched whole.
   @> Chambers (probe 1.40 Å): 1 of the 7 searched cavities have them; the other 6 are searched whole.
   @>     cavity 0: 26 chambers, 10 of them seeded.
   @> 16 search sites (sp) in 0.07s: one per seeded chamber, one per cavity searched whole.
   @> Chambers (probe 1.40 Å): 2 of the 6 searched cavities have them; the other 5 are searched whole.
   @>     cavity 0: 15 chambers, 8 of them seeded.
   @>     cavity 3: 1 chamber, none of them deep and large enough to seed; searched whole.
   @> 13 search sites (sp) in 0.09s: one per seeded chamber, one per cavity searched whole.
   @> Channel search (Dijkstra) over 20 search sites in 14 cavities completed in 2.09s.
   @> Found 35 channels and 3 links (a link joins a deep chamber to a shallower one and never reaches the surface).
   @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
   @>     site  void                   volume [Å³]  depth [Å]  channels  links
   @>     sp0   cavity 0, chamber 1/5         1452       17.0         2      1  -> sp9
   @>     sp1   cavity 2, whole                862        5.3         1      -
   @>     sp2   cavity 3, whole                286        9.4         2      -
   @>     sp3   cavity 4, whole                259        5.2         2      -
   @>     sp4   cavity 5, whole                254        6.6         2      -
   @>     sp5   cavity 6, whole                198        6.7         -      -  sealed
   @>     sp6   cavity 7, whole                158        6.2         1      -
   @>     sp7   cavity 8, whole                138        6.5         1      -
   @>     sp8   cavity 9, whole                126        6.0         -      -  sealed
   @>     sp9   cavity 0, chamber 2/5          108        5.3         2      -
   @>     sp10  cavity 10, whole               103        8.1         -      -  sealed
   @>     sp11  cavity 11, whole                93        8.2         -      -  sealed
   @>     sp12  cavity 12, whole                80        6.1         -      -  sealed
   @>     sp13  cavity 1, chamber 1/3           68        7.6         2      -
   @>     sp14  cavity 0, chamber 3/5           67        6.6         5      -
   @>     sp15  cavity 1, chamber 2/3           57       11.3         -      -  sealed
   @>     sp16  cavity 13, whole                56        6.0         1      -
   @>     sp17  cavity 0, chamber 4/5           34        5.8        10      -
   @>     sp18  cavity 0, chamber 5/5           31       11.9         2      1  -> sp0
   @>     sp19  cavity 1, chamber 3/3           31       12.9         2      1  -> sp13
   @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
   @> The 1 site marked sealed above report neither a channel nor a link: every route out of them is narrower than bottleneck=0.80 Å. Lower it to see how they connect.
   @> Saving 35 channels and 3 links to directory ., one file per object named sp<site>_chl<n> and sp<site>_lnk<n>.
   @> Channel calculation completed in 11.39s.
   @> Frame/model: 15
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Channel search (Dijkstra) over 13 search sites in 6 cavities completed in 2.19s.
   @> Found 40 channels and 5 links (a link joins a deep chamber to a shallower one and never reaches the surface).
   @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
   @>     site  void                   volume [Å³]  depth [Å]  channels  links
   @>     sp0   cavity 0, chamber 1/8         2748       13.7        14      2  -> sp8, sp6
   @>     sp1   cavity 1, whole                567       11.2         3      -
   @>     sp2   cavity 2, whole                427       24.2         -      -  sealed
   @>     sp3   cavity 3, whole                364        5.9         1      -
   @>     sp4   cavity 4, whole                362        5.4         1      -
   @>     sp5   cavity 0, chamber 2/8          256       12.1         2      -
   @>     sp6   cavity 0, chamber 3/8          154        5.3         3      -
   @>     sp7   cavity 5, whole                154        8.9         -      -  sealed
   @>     sp8   cavity 0, chamber 4/8          148       10.0         1      -
   @>     sp9   cavity 0, chamber 5/8          115        5.1         9      -
   @>     sp10  cavity 0, chamber 6/8           59        6.9         4      1  -> sp0
   @>     sp11  cavity 0, chamber 7/8           50       17.5         -      2  -> sp0, sp5
   @>     sp12  cavity 0, chamber 8/8           31        5.4         2      -
   @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
   @> Saving 40 channels and 5 links to directory ., one file per object named sp<site>_chl<n> and sp<site>_lnk<n>.
   @> Channel calculation completed in 11.55s.
   ..
   ..
   @> Frame/model: 209
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 Å in 0.27s.
   @> Surface and inner simplices filtered in 2.72s.
   @> Surface cavities: 429 found, 10 deeper than min_depth=5.0 Å and searched for channels, in 0.72s.
   @> Channel search (Dijkstra) over 17 search sites in 10 cavities completed in 2.07s.
   @> Found 32 channels and 6 links (a link joins a deep chamber to a shallower one and never reaches the surface).
   @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
   @>     site  void                   volume [Å³]  depth [Å]  channels  links
   @>     sp0   cavity 0, chamber 1/8         3008       11.1        13      -
   @>     sp1   cavity 2, whole                374        5.1         1      -
   @>     sp2   cavity 3, whole                316        7.0         -      -  sealed
   @>     sp3   cavity 4, whole                306        7.2         -      -  sealed
   @>     sp4   cavity 0, chamber 2/8          303       12.7         1      2  -> sp0, sp13
   @>     sp5   cavity 5, whole                300        8.6         2      -
   @>     sp6   cavity 7, whole                139        5.2         1      -
   @>     sp7   cavity 8, whole                120        6.3         1      -
   @>     sp8   cavity 0, chamber 3/8           87        5.0         3      -
   @>     sp9   cavity 9, whole                 64        6.2         -      -  sealed
   @>     sp10  cavity 0, chamber 4/8           50       25.5         -      1  -> sp0
   @>     sp11  cavity 0, chamber 5/8           48       15.1         3      1  -> sp0
   @>     sp12  cavity 1, chamber 1/1           42        5.4         6      -
   @>     sp13  cavity 0, chamber 6/8           36        7.7         1      -
   @>     sp14  cavity 6, chamber 1/1           35        6.9         -      -  sealed
   @>     sp15  cavity 0, chamber 7/8           31       26.7         -      1  -> sp0
   @>     sp16  cavity 0, chamber 8/8           30       20.7         -      1  -> sp0
   @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
   @> The 1 site marked sealed above report neither a channel nor a link: every route out of them is narrower than bottleneck=0.80 Å. Lower it to see how they connect.
   @> Saving 32 channels and 6 links to directory ., one file per object named sp<site>_chl<n> and sp<site>_lnk<n>.
   @> Channel calculation completed in 9.70s.
   @> Chambers (probe 1.40 Å): 3 of the 10 searched cavities have them; the other 8 are searched whole.
   @>     cavity 0: 17 chambers, 6 of them seeded.
   @>     cavity 1: 3 chambers, 1 of them seeded.
   @>     cavity 2: 1 chamber, none of them deep and large enough to seed; searched whole.
   @> 15 search sites (sp) in 0.09s: one per seeded chamber, one per cavity searched whole.
   @> Delaunay tessellation of 77434 points constructed in 2.95s.
   @> Channel search (Dijkstra) over 15 search sites in 10 cavities completed in 1.91s.
   @> Found 42 channels and 6 links (a link joins a deep chamber to a shallower one and never reaches the surface).
   @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
   @>     site  void                   volume [Å³]  depth [Å]  channels  links
   @>     sp0   cavity 0, chamber 1/6         2362       11.8        16      -
   @>     sp1   cavity 0, chamber 2/6          312        7.2         3      1  -> sp0
   @>     sp2   cavity 2, whole                277        7.3         -      -  sealed
   @>     sp3   cavity 0, chamber 3/6          252       12.2         4      3  -> sp0, sp8, sp6
   @>     sp4   cavity 3, whole                251       14.2         -      -  sealed
   @>     sp5   cavity 4, whole                216        5.0         2      -
   @>     sp6   cavity 0, chamber 4/6          198        8.3         1      1  -> sp0
   @>     sp7   cavity 1, chamber 1/1          153        7.0         6      -
   @>     sp8   cavity 0, chamber 5/6          104        7.9         3      1  -> sp0
   @>     sp9   cavity 5, whole                103        5.6         -      -  sealed
   @>     sp10  cavity 6, whole                 77        5.3         -      -  sealed
   @>     sp11  cavity 7, whole                 75        5.1         -      -  sealed
   @>     sp12  cavity 0, chamber 6/6           68        7.1         7      -
   @>     sp13  cavity 8, whole                 68        5.8         -      -  sealed
   @>     sp14  cavity 9, whole                 65        5.8         -      -  sealed
   @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
   @> Saving 42 channels and 6 links to directory ., one file per object named sp<site>_chl<n> and sp<site>_lnk<n>.
   @> Channel calculation completed in 9.42s.
   @> Surface and inner simplices filtered in 2.73s.
   @> Surface cavities: 467 found, 10 deeper than min_depth=5.0 Å and searched for channels, in 0.66s.
   @> Chambers (probe 1.40 Å): 2 of the 10 searched cavities have them; the other 8 are searched whole.
   @>     cavity 0: 14 chambers, 8 of them seeded.
   @>     cavity 1: 1 chamber, seeded.
   @> 17 search sites (sp) in 0.07s: one per seeded chamber, one per cavity searched whole.
   @> Channel search (Dijkstra) over 17 search sites in 10 cavities completed in 1.58s.
   @> Found 46 channels and 9 links (a link joins a deep chamber to a shallower one and never reaches the surface).
   @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
   @>     site  void                   volume [Å³]  depth [Å]  channels  links
   @>     sp0   cavity 0, chamber 1/8         2514       10.2        15      2  -> sp3, sp7
   @>     sp1   cavity 2, whole                353        6.1         1      -
   @>     sp2   cavity 3, whole                322        5.3         1      -
   @>     sp3   cavity 0, chamber 2/8          250        5.1         5      -
   @>     sp4   cavity 4, whole                207        5.6         1      -
   @>     sp5   cavity 0, chamber 3/8          202       26.7         -      1  -> sp13
   @>     sp6   cavity 0, chamber 4/8          195        7.6         6      2  -> sp0, sp0
   @>     sp7   cavity 0, chamber 5/8          141        6.5         7      -
   @>     sp8   cavity 5, whole                116        6.3         1      -
   @>     sp9   cavity 6, whole                 69        5.2         1      -
   @>     sp10  cavity 7, whole                 64        5.1         -      -  sealed
   @>     sp11  cavity 8, whole                 62        5.0         -      -  sealed
   @>     sp12  cavity 9, whole                 62        5.1         1      -
   @>     sp13  cavity 0, chamber 6/8           60       20.4         1      2  -> sp0, sp0
   @>     sp14  cavity 0, chamber 7/8           48       15.1         2      1  -> sp0
   @>     sp15  cavity 1, chamber 1/1           45        5.5         4      -
   @>     sp16  cavity 0, chamber 8/8           33       31.3         -      1  -> sp5
   @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
   @> Saving 46 channels and 9 links to directory ., one file per object named sp<site>_chl<n> and sp<site>_lnk<n>.
   @> Channel calculation completed in 8.70s.


Once the channels are identified, a function called
:func:`.calcPoresFromChannelsMultipleFrames` can be applied to reconstrct
pores. To eliminate certain pores an additional filters will be applied,
such as ``min_end_to_end`` and ``min_bottleneck``.

.. ipython:: python
   :verbatim:

   pores = calcPoresFromChannelsMultipleFrames(channels, details, 
					min_end_to_end=45, 
                                	output_path='pores_dcd_', 
					separate=True, 
					min_bottleneck=0.6, 
					max_proc=4)

.. parsed-literal::
   
   @> Frame/model: 0
   @> Frame/model: 1
   ..
   ..
   @> Frame/model: 204
   @> Frame/model: 205
   @> Frame/model: 206
   @> Frame/model: 207
   @> Frame/model: 208
   @> Frame/model: 209


.. ipython:: python
   :verbatim:

   pores

.. parsed-literal::

   [[<prody.proteins.channels.Channel at 0x723a25804250>,
     <prody.proteins.channels.Channel at 0x723a258043d0>,
     <prody.proteins.channels.Channel at 0x723a258044c0>,
     <prody.proteins.channels.Channel at 0x723a25804610>,
     <prody.proteins.channels.Channel at 0x723a25804700>,
     <prody.proteins.channels.Channel at 0x723a25804820>,
     <prody.proteins.channels.Channel at 0x723a258048e0>,
     <prody.proteins.channels.Channel at 0x723a25804a00>,
     <prody.proteins.channels.Channel at 0x723a25804af0>,
     <prody.proteins.channels.Channel at 0x723a25804c10>,
     <prody.proteins.channels.Channel at 0x723a25804d00>,
     <prody.proteins.channels.Channel at 0x723a25804e50>,
     <prody.proteins.channels.Channel at 0x723a25804f10>,
     <prody.proteins.channels.Channel at 0x723a25805060>,
     <prody.proteins.channels.Channel at 0x723a25805120>,
     <prody.proteins.channels.Channel at 0x723a25805240>,
     <prody.proteins.channels.Channel at 0x723a25805390>],
    [<prody.proteins.channels.Channel at 0x723a25805480>,
     <prody.proteins.channels.Channel at 0x723a258055a0>,
     <prody.proteins.channels.Channel at 0x723a25805690>,
     <prody.proteins.channels.Channel at 0x723a25805780>,
     <prody.proteins.channels.Channel at 0x723a258058a0>,
     <prody.proteins.channels.Channel at 0x723a25805990>,
     <prody.proteins.channels.Channel at 0x723a25805a80>,
     ..
     <prody.proteins.channels.Channel at 0x723a76160c40>,
     <prody.proteins.channels.Channel at 0x723a76161ed0>,
     <prody.proteins.channels.Channel at 0x723a76161bd0>,
     <prody.proteins.channels.Channel at 0x723a76162ce0>,
     <prody.proteins.channels.Channel at 0x723a76160130>,
     <prody.proteins.channels.Channel at 0x723a76160430>]]


.. ipython:: python
   :verbatim:

   getPoreParametersMultipleFrames(pores, param_file_name='pores_DATA')

.. parsed-literal::

   @> Frame/model: 0
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 1
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	707.21 		83.72 		0.82
   @> pore 1: 	724.46 		85.46 		0.82
   @> Frame/model: 2
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 3
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 4
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 5
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 6
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 7
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 8
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 9
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 10
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 11
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 12
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 13
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 14
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 15
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 16
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 17
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 18
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 19
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 20
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 21
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 22
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 23
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 24
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 25
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 26
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 27
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 28
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 29
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 30
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 31
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 32
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	701.25 		69.59 		0.9
   @> pore 1: 	709.49 		70.43 		0.9
   @> pore 2: 	709.35 		72.61 		0.9
   @> pore 3: 	724.86 		72.81 		0.9
   @> pore 4: 	656.23 		68.74 		0.9
   @> pore 5: 	664.47 		69.57 		0.9
   @> pore 6: 	664.33 		71.75 		0.9
   @> pore 7: 	667.89 		70.18 		0.9
   @> pore 8: 	676.13 		71.02 		0.9
   @> pore 9: 	675.99 		73.2 		0.9
   @> pore 10: 	667.69 		71.78 		0.83
   @> pore 11: 	675.93 		72.61 		0.83
   @> pore 12: 	675.79 		74.79 		0.83
   @> pore 13: 	678.16 		74.01 		0.9
   @> pore 14: 	686.26 		77.03 		0.9
   @> Frame/model: 33
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 34
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	689.31 		66.51 		0.84
   @> pore 1: 	855.15 		76.99 		0.84
   @> pore 2: 	877.82 		78.06 		0.84
   @> Frame/model: 35
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 36
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 37
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	861.31 		72.16 		0.84
   @> pore 1: 	858.31 		73.62 		0.84
   @> pore 2: 	847.59 		73.02 		0.84
   @> pore 3: 	865.16 		75.06 		0.84
   @> Frame/model: 38
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 39
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 40
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 41
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 42
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 43
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 44
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 45
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 46
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 47
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	820.01 		65.24 		0.93
   @> Frame/model: 48
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 49
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 50
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 51
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 52
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 53
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 54
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 55
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 56
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 57
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	893.65 		60.69 		0.86
   @> Frame/model: 58
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 59
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 60
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 61
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	930.42 		66.35 		0.83
   @> Frame/model: 62
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 63
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 64
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 65
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 66
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 67
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 68
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 69
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	1176.59 	70.54 		0.91
   @> pore 1: 	1164.71 	70.27 		0.82
   @> pore 2: 	1108.46 	75.52 		0.81
   @> pore 3: 	1186.28 	72.98 		0.91
   @> pore 4: 	1174.4 		72.7 		0.82
   @> pore 5: 	1057.83 	77.63 		0.81
   @> pore 6: 	1167.48 	72.58 		0.89
   @> pore 7: 	1155.6 		72.3 		0.82
   @> pore 8: 	1099.35 	77.56 		0.81
   @> pore 9: 	1076.07 	79.94 		0.81
   @> Frame/model: 70
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 71
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 72
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 73
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 74
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 75
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 76
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 77
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 78
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 79
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 80
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 81
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	795.21 		67.59 		0.81
   @> Frame/model: 82
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 83
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	1108.72 	65.41 		0.88
   @> Frame/model: 84
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 85
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 86
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 87
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 88
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 89
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 90
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 91
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 92
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 93
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	1144.18 	68.3 		0.86
   @> pore 1: 	1060.45 	68.18 		0.86
   @> pore 2: 	1085.82 	69.92 		0.86
   @> pore 3: 	942.95 		68.95 		0.86
   @> Frame/model: 94
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 95
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 96
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 97
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	830.66 		60.01 		0.86
   @> pore 1: 	833.23 		60.23 		0.83
   @> pore 2: 	858.98 		66.3 		0.85
   @> pore 3: 	805.77 		60.0 		0.87
   @> pore 4: 	808.33 		60.23 		0.83
   @> pore 5: 	834.08 		66.3 		0.85
   @> pore 6: 	782.83 		63.5 		0.86
   @> pore 7: 	785.39 		63.73 		0.83
   @> pore 8: 	801.0 		67.76 		0.84
   @> pore 9: 	811.14 		69.8 		0.85
   @> pore 10: 	798.51 		66.21 		0.87
   @> pore 11: 	801.08 		66.44 		0.83
   @> pore 12: 	816.69 		70.47 		0.84
   @> pore 13: 	826.83 		72.5 		0.85
   @> Frame/model: 98
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 99
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	1238.29 	69.38 		0.99
   @> pore 1: 	1213.77 	69.32 		1.05
   @> pore 2: 	1154.92 	72.43 		0.9
   @> Frame/model: 100
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	989.59 		63.76 		0.91
   @> pore 1: 	860.23 		64.99 		0.89
   @> pore 2: 	862.52 		66.36 		0.8
   @> Frame/model: 101
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	968.66 		70.49 		0.9
   @> pore 1: 	934.0 		71.32 		0.84
   @> Frame/model: 102
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 103
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 104
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 105
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 106
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 107
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 108
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 109
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 110
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	966.59 		65.25 		1.02
   @> pore 1: 	962.3 		66.72 		1.04
   @> pore 2: 	967.57 		68.76 		1.05
   @> pore 3: 	943.48 		69.54 		0.91
   @> pore 4: 	1016.2 		68.52 		1.04
   @> pore 5: 	1021.47 	70.56 		1.05
   @> pore 6: 	980.77 		68.13 		1.04
   @> pore 7: 	986.04 		70.16 		1.05
   @> pore 8: 	960.28 		67.83 		1.04
   @> pore 9: 	965.55 		69.86 		1.05
   @> Frame/model: 111
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	1034.87 	62.91 		1.06
   @> pore 1: 	1065.49 	63.72 		1.06
   @> pore 2: 	1053.66 	63.97 		1.06
   @> pore 3: 	1041.08 	68.13 		0.81
   @> Frame/model: 112
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	948.16 		64.22 		0.81
   @> pore 1: 	864.89 		64.85 		0.81
   @> Frame/model: 113
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 114
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 115
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 116
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	876.28 		69.2 		0.82
   @> pore 1: 	889.76 		74.78 		0.82
   @> pore 2: 	898.28 		75.65 		0.82
   @> pore 3: 	911.33 		76.67 		0.82
   @> Frame/model: 117
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	998.54 		64.7 		0.88
   @> pore 1: 	983.81 		65.91 		0.88
   @> pore 2: 	973.74 		72.61 		0.88
   @> Frame/model: 118
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	1046.85 	68.35 		0.92
   @> pore 1: 	1040.18 	68.24 		0.81
   @> pore 2: 	980.96 		71.31 		0.92
   @> pore 3: 	974.29 		71.19 		0.81
   @> pore 4: 	1054.73 	68.35 		0.99
   @> pore 5: 	1004.6 		68.69 		0.99
   @> pore 6: 	1032.15 	74.8 		0.92
   @> pore 7: 	1025.48 	74.69 		0.81
   @> pore 8: 	1040.33 	75.56 		0.9
   @> pore 9: 	1062.8 		69.78 		0.99
   @> pore 10: 	1012.67 	70.13 		0.99
   @> pore 11: 	1040.22 	76.24 		0.92
   @> pore 12: 	1033.55 	76.12 		0.81
   @> pore 13: 	1048.4 		77.0 		0.9
   @> pore 14: 	742.19 		68.97 		0.81
   @> pore 15: 	884.19 		67.34 		0.85
   @> pore 16: 	834.06 		67.69 		0.85
   @> pore 17: 	861.61 		73.78 		0.85
   @> pore 18: 	854.94 		73.67 		0.81
   @> pore 19: 	869.79 		74.55 		0.85
   @> pore 20: 	791.32 		66.88 		0.85
   @> pore 21: 	741.2 		67.22 		0.85
   @> pore 22: 	768.74 		73.33 		0.85
   @> pore 23: 	762.08 		73.21 		0.81
   @> pore 24: 	776.93 		74.09 		0.85
   @> Frame/model: 119
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 120
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 121
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 122
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 123
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 124
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 125
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 126
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	943.17 		68.97 		0.82
   @> pore 1: 	870.09 		71.6 		0.82
   @> pore 2: 	875.44 		76.61 		0.82
   @> pore 3: 	877.86 		78.0 		0.82
   @> pore 4: 	880.76 		78.43 		0.82
   @> Frame/model: 127
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 128
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 129
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 130
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 131
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 132
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	712.66 		76.07 		0.81
   @> Frame/model: 133
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 134
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 135
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 136
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 137
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 138
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 139
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 140
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 141
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 142
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 143
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 144
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 145
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 146
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 147
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 148
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 149
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	850.56 		66.88 		0.88
   @> pore 1: 	846.51 		65.5 		0.89
   @> pore 2: 	843.03 		68.95 		0.88
   @> pore 3: 	823.02 		69.18 		0.88
   @> pore 4: 	805.46 		69.65 		0.88
   @> Frame/model: 150
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	826.42 		65.89 		0.9
   @> pore 1: 	808.95 		65.19 		0.9
   @> pore 2: 	856.57 		67.05 		0.9
   @> pore 3: 	815.74 		66.1 		0.9
   @> pore 4: 	776.43 		65.71 		0.9
   @> pore 5: 	758.96 		65.01 		0.9
   @> pore 6: 	806.58 		66.87 		0.9
   @> pore 7: 	765.75 		65.92 		0.9
   @> pore 8: 	819.02 		68.64 		0.9
   @> pore 9: 	801.55 		67.94 		0.9
   @> pore 10: 	849.17 		69.8 		0.9
   @> pore 11: 	808.35 		68.85 		0.9
   @> Frame/model: 151
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 152
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 153
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 154
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 155
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	886.15 		67.7 		0.88
   @> pore 1: 	852.29 		62.77 		0.99
   @> pore 2: 	887.65 		63.75 		0.99
   @> pore 3: 	870.92 		65.88 		0.99
   @> pore 4: 	861.77 		69.53 		0.99
   @> pore 5: 	860.44 		69.12 		0.91
   @> pore 6: 	851.76 		69.59 		0.88
   @> pore 7: 	854.81 		70.82 		0.88
   @> pore 8: 	887.56 		66.04 		0.99
   @> pore 9: 	922.92 		67.02 		0.99
   @> pore 10: 	906.19 		69.15 		0.99
   @> pore 11: 	897.04 		72.8 		0.99
   @> pore 12: 	895.7 		72.4 		0.91
   @> pore 13: 	887.02 		72.86 		0.88
   @> pore 14: 	985.69 		73.86 		0.99
   @> pore 15: 	1021.04 	74.84 		0.99
   @> pore 16: 	1004.31 	76.97 		0.99
   @> pore 17: 	995.17 		80.62 		0.99
   @> pore 18: 	993.83 		80.21 		0.91
   @> pore 19: 	985.15 		80.68 		0.88
   @> pore 20: 	988.2 		81.91 		0.88
   @> pore 21: 	947.58 		73.15 		0.99
   @> pore 22: 	982.93 		74.13 		0.99
   @> pore 23: 	966.2 		76.26 		0.99
   @> pore 24: 	957.06 		79.91 		0.99
   @> pore 25: 	955.72 		79.5 		0.91
   @> pore 26: 	947.04 		79.97 		0.88
   @> pore 27: 	950.1 		81.2 		0.88
   @> pore 28: 	906.84 		71.65 		0.99
   @> pore 29: 	942.2 		72.63 		0.99
   @> pore 30: 	925.47 		74.76 		0.99
   @> pore 31: 	916.32 		78.41 		0.99
   @> pore 32: 	914.98 		78.0 		0.91
   @> pore 33: 	906.31 		78.47 		0.88
   @> pore 34: 	909.36 		79.7 		0.88
   @> pore 35: 	856.95 		76.3 		0.88
   @> pore 36: 	856.12 		75.86 		0.88
   @> pore 37: 	680.2 		72.13 		0.81
   @> pore 38: 	678.86 		71.72 		0.81
   @> pore 39: 	670.19 		72.18 		0.81
   @> pore 40: 	915.94 		80.73 		0.96
   @> pore 41: 	906.79 		84.38 		0.96
   @> pore 42: 	905.46 		83.98 		0.91
   @> pore 43: 	896.78 		84.44 		0.88
   @> pore 44: 	688.33 		73.15 		0.81
   @> pore 45: 	686.99 		72.74 		0.81
   @> pore 46: 	678.31 		73.2 		0.81
   @> pore 47: 	885.53 		83.95 		0.86
   ..
   ..


   [([], [], []),
    ([83.72181023903833, 85.46162472147293],
     [0.8227300033137241, 0.8227300033137241],
     [707.2123815585367, 724.4646694011919]),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([69.5949898851407,
      70.42852451372993,
      72.60822067706627,
      72.81313607212425,
      68.73784984146138,
      69.57114702715292,
      71.75016583222246,
      70.18230030971134,
      71.01620131367125,
      73.19713589126198,
      71.78021275467378,
      72.61364881055758,
      74.7936466868367,
      74.00589255027344,
      77.026290787233],
     [0.9034884486107735,
      0.9034884486107735,
      0.9034884486107735,
      0.9034884486107735,
      0.9034884486107735,
      0.9034884486107735,
      0.9034884486107735,
      0.9034884486107735,
      0.9034884486107735,
      0.9034884486107735,
      0.825794417769907,
      0.825794417769907,
      0.825794417769907,
      0.9034884486107735,
      0.9034884486107735],
     [701.2501227800298,
      709.4917858287764,
      709.3492774145426,
      724.8622564896269,
      656.2258665272245,
      664.4675295759711,
      664.3250211617374,
      667.8930749725629,
      676.1347380213094,
      675.9922296070758,
      667.6876949427416,
      675.9293579914882,
      675.7868495772544,
      678.1626699366484,
      686.2618245711614]),
    ([], [], []),
    ([66.50699517343948, 76.99357439825909, 78.05759150154641],
     [0.8411005602059334, 0.8411005602059334, 0.8411005602059334],
     [689.3078797713896, 855.1473499011635, 877.821764326892]),
    ([], [], []),
    ([], [], []),
    ([72.1598121093621, 73.6212787291321, 73.01609384001452, 75.06090159160428],
     [0.8383209859042193,
      0.8383209859042193,
      0.8383209859042193,
      0.8383209859042193],
     [861.3065768018067,
      858.3125545315595,
      847.5870240474648,
      865.1578695365567]),
     ..
     ..


To visualize the results directly in ProDy create the model and use
:func:`.showPores` function.

.. ipython:: python
   :verbatim:

   vmd_path = '/usr/local/bin/vmd'
   model = getVmdModel(vmd_path, protein)

.. parsed-literal::

   @> Model created successfully.

	
To display all the pores from frame #1:

.. ipython:: python
   :verbatim:

   showPores(pores[1], model=model)


.. figure:: images/cavitracer_figure29.jpg
   :scale: 50 %

To display first pore in frame #1:

.. ipython:: python
   :verbatim:

   showPores(pores[1][0], model=model)

.. figure:: images/cavitracer_figure30.jpg
   :scale: 50 %


To display second pore in frame #1:

.. ipython:: python
   :verbatim:

   showPores(pores[1][1], model=model)

.. figure:: images/cavitracer_figure31.jpg
   :scale: 50 %

Below are also the results for frames #101 and #201 (counting from 0). 
As we can see, pores are changing along the MD trajectory.

.. ipython:: python
   :verbatim:

   showPores(pores[100], model=model)

.. figure:: images/cavitracer_figure33.jpg
   :scale: 50 %


.. ipython:: python
   :verbatim:

   showPores(pores[150], model=model)

.. figure:: images/cavitracer_figure32.jpg
   :scale: 50 %

Next, the residues that are forming the pores can be identified using
:func:`.getPoreResidueNamesMultipleFrames` function.

.. ipython:: python
   :verbatim:

   getPoreResidueNamesMultipleFrames(protein, pores, dcd,
		residues_file_name='pores_Residues')

.. parsed-literal::

   @> Frame: 0
   @> Pore residues were saved to: pores_Residues_frame0_Residues_All_pores.txt
   @> Frame: 1
   @> Pore residues were saved to: pores_Residues_frame1_Residues_All_pores.txt
   @> Frame: 2
   @> Pore residues were saved to: pores_Residues_frame2_Residues_All_pores.txt
   @> Frame: 3
   @> Pore residues were saved to: pores_Residues_frame3_Residues_All_pores.txt
   @> Frame: 4
   @> Pore residues were saved to: pores_Residues_frame4_Residues_All_pores.txt
   @> Frame: 5
   @> Pore residues were saved to: pores_Residues_frame5_Residues_All_pores.txt
   @> Frame: 6
   @> Pore residues were saved to: pores_Residues_frame6_Residues_All_pores.txt
   @> Frame: 7
   @> Pore residues were saved to: pores_Residues_frame7_Residues_All_pores.txt
   @> Frame: 8
   @> Pore residues were saved to: pores_Residues_frame8_Residues_All_pores.txt
   @> Frame: 9
   @> Pore residues were saved to: pores_Residues_frame9_Residues_All_pores.txt
   @> Frame: 10
   @> Pore residues were saved to: pores_Residues_frame10_Residues_All_pores.txt
   ..
   ..


    ['pore0: ASP33:P, ASN34:P, LEU37:P, THR38:P, VAL41:P, ILE44:P, GLU120:P, ASP121:P, LYS122:P, ASP123:P, GLU127:P, ASN128:P, VAL131:P, GLY132:P, PHE135:P, LYS138:P, GLN142:P, SER196:P, SER200:P, MET204:P, LEU207:P, ALA208:P, TYR211:P, THR212:P, ARG217:P, VAL220:P, MET221:P, LEU225:P, LEU228:P, VAL232:P, GLU312:P, PRO313:P, LEU315:P, PRO316:P, ILE317:P, TRP318:P, PHE334:P, TYR341:P, MET403:P, GLY407:P, ASP411:P, TYR418:P, TYR422:P, ALA425:P, ASP426:P, PHE429:P, TYR433:P',
     'pore1: ASP33:P, ASN34:P, LEU37:P, THR38:P, VAL41:P, ILE44:P, GLU120:P, ASP121:P, LYS122:P, ASP123:P, GLU127:P, ASN128:P, VAL131:P, GLY132:P, PHE135:P, LYS138:P, GLN142:P, SER196:P, SER200:P, MET204:P, LEU207:P, ALA208:P, TYR211:P, ASP214:P, GLU215:P, ARG217:P, GLY218:P, VAL220:P, MET221:P, LEU225:P, LEU228:P, VAL232:P, GLU312:P, PRO313:P, LEU315:P, PRO316:P, ILE317:P, TRP318:P, PHE334:P, TYR341:P, ARG357:P, MET403:P, PRO404:P, GLY407:P, TYR422:P, ALA425:P, ASP426:P, PHE429:P, TYR433:P',
     'pore2: ASP33:P, ASN34:P, LEU37:P, THR38:P, VAL41:P, ILE44:P, GLU120:P, ASP121:P, LYS122:P, ASP123:P, GLU127:P, ASN128:P, VAL131:P, GLY132:P, PHE135:P, LYS138:P, GLN142:P, SER196:P, SER200:P, MET204:P, LEU207:P, ALA208:P, TYR211:P, ASP214:P, GLU215:P, ARG217:P, GLY218:P, VAL220:P, MET221:P, LEU225:P, LEU228:P, VAL232:P, GLU312:P, PRO313:P, LEU315:P, PRO316:P, ILE317:P, TRP318:P, PHE334:P, TYR341:P, ARG357:P, MET403:P, PRO404:P, GLY407:P, TYR422:P, ALA425:P, ASP426:P, PHE429:P, TYR433:P',
     'pore3: ASP33:P, ASN34:P, LEU37:P, THR38:P, VAL41:P, PRO42:P, ILE44:P, SER46:P, ASN128:P, VAL131:P, GLY132:P, PHE135:P, LYS138:P, GLN142:P, SER196:P, SER200:P, MET204:P, LEU207:P, ALA208:P, TYR211:P, THR212:P, ARG217:P, VAL220:P, MET221:P, LEU225:P, LEU228:P, VAL232:P, GLU312:P, PRO313:P, ILE317:P, TRP318:P, MET319:P, GLU321:P, THR322:P, MET323:P, ARG326:P, LYS327:P, TRP328:P, GLN329:P, LEU330:P, PHE334:P, TYR341:P, MET403:P, GLY407:P, ASP411:P, TYR418:P, TYR422:P, ALA425:P, ASP426:P, PHE429:P, TYR433:P'],
    [],
    [],
    [],
    ['pore0: LEU30:P, ASN34:P, LEU37:P, THR38:P, VAL41:P, ILE44:P, GLU120:P, ASN128:P, GLN130:P, VAL131:P, GLY132:P, PHE135:P, LYS138:P, GLN142:P, SER200:P, MET204:P, ALA208:P, TYR211:P, THR212:P, ARG217:P, MET221:P, ALA224:P, LEU225:P, LEU228:P, VAL232:P, ILE308:P, GLU312:P, PRO313:P, ILE317:P, TRP318:P, PHE334:P, TYR341:P, MET403:P, PRO404:P, GLY407:P, ASP411:P, TYR418:P, TYR422:P, ASP426:P, PHE429:P, TYR433:P',
     'pore1: LEU30:P, ASN34:P, LEU37:P, THR38:P, VAL41:P, PRO42:P, ILE43:P, ILE44:P, SER46:P, VAL131:P, GLY132:P, PHE135:P, LYS138:P, GLN142:P, SER200:P, MET204:P, ALA208:P, TYR211:P, THR212:P, ARG217:P, MET221:P, ALA224:P, LEU225:P, LEU228:P, VAL232:P, TYR243:P, GLU244:P, ILE308:P, GLU312:P, PRO313:P, ILE317:P, TRP318:P, MET319:P, MET320:P, THR322:P, MET323:P, SER325:P, LEU330:P, PHE334:P, TYR341:P, MET403:P, PRO404:P, GLY407:P, ASP411:P, TYR418:P, TYR422:P, ASP426:P, PHE429:P, TYR433:P',
     'pore2: LEU30:P, ASN34:P, LEU37:P, THR38:P, VAL41:P, PRO42:P, ILE44:P, VAL131:P, GLY132:P, PHE135:P, LYS138:P, GLN142:P, SER200:P, MET204:P, ALA208:P, TYR211:P, THR212:P, ARG217:P, MET221:P, ALA224:P, LEU225:P, LEU228:P, VAL232:P, ILE308:P, GLU312:P, PRO313:P, ILE317:P, TRP318:P, MET319:P, GLU321:P, THR322:P, ARG326:P, LYS327:P, TRP328:P, GLN329:P, LEU330:P, PHE334:P, TYR341:P, MET403:P, PRO404:P, GLY407:P, ASP411:P, TYR418:P, TYR422:P, ASP426:P, PHE429:P, TYR433:P',
     'pore3: LEU30:P, ASN34:P, LEU37:P, THR38:P, VAL41:P, PRO42:P, ILE44:P, VAL131:P, GLY132:P, PHE135:P, LYS138:P, GLN142:P, SER200:P, MET204:P, ALA208:P, TYR211:P, THR212:P, ARG217:P, MET221:P, ALA224:P, LEU225:P, LEU228:P, VAL232:P, ILE308:P, GLU312:P, PRO313:P, ILE317:P, TRP318:P, MET319:P, GLU321:P, THR322:P, ARG326:P, LYS327:P, TRP328:P, GLN329:P, LEU330:P, PHE334:P, TYR341:P, MET403:P, PRO404:P, GLY407:P, ASP411:P, TYR418:P, TYR422:P, ASP426:P, PHE429:P, TYR433:P'],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    []]


III. Identification of surface cavities in molecular dynamics (MD) trajectory
===============================================================================


In this example, surface cavities are analyzed for a trajectory prepared from
the structure with PDB ID 5KQM. This structure corresponds to human low molecular
weight phosphotyrosine protein phosphatase (LMW-PTP, ACP1) in complex with MES.

The trajectory used here contains both the protein without the ligand (MES). 
The protein selection is applied to the trajectory using :meth:`Trajectory.setAtoms`, 
so that only protein coordinates are passed to the surface-cavity calculation.

.. ipython:: python
   :verbatim:
   
   PDBfile = '5kqm_all_sci_wWat.pdb'
   DCDfile = '5kqm_all_sci_wWat.dcd'
   atoms = parsePDB(PDBfile)
   dcd = Trajectory(DCDfile)
   dcd.link(atoms)
   dcd.setCoords(atoms)

.. ipython:: python
   :verbatim:

   protein = atoms.select("protein")
   dcd.setAtoms(protein)


.. parsed-literal::

   @> 2425 atoms and 1 coordinate set(s) were parsed in 0.02s.


The selected protein and the linked trajectory are then used as input for
:func:`.calcSurfaceCavitiesMultipleFrames`. The function calculates surface
cavities independently for each trajectory frame. With ``separate=True``, each
detected cavity is saved as an individual PQR file for each frame, using the
prefix provided by ``output_path``.

.. ipython:: python
   :verbatim:

   cavities, surfaces=calcSurfaceCavitiesMultipleFrames(protein, dcd, 
		output_path='cav_dcd', separate=True)
		

.. parsed-literal::

   @> Frame/model: 0
   @> Frame/model: 7
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 2425 atoms with 31833 homogeneous balls of radius 1.20 Å in 0.17s.
   @> Substituted 2425 atoms with 31833 homogeneous balls of radius 1.20 Å in 0.17s.
   @> Delaunay tessellation of 31833 points constructed in 1.40s.
   @> Delaunay tessellation of 31833 points constructed in 1.40s.
   @> Surface and inner simplices filtered in 0.41s.
   @> Surface and inner simplices filtered in 0.44s.
   @> Surface cavities: 228 found, 25 deeper than min_depth=1.5 Å and kept, in 0.23s.
   @> Surface cavities: 222 found, 22 deeper than min_depth=1.5 Å and kept, in 0.21s.
   @> Returning surface cavities
   @> Saving multiple surface cavities to directory ..
   @> Surface cavity calculation completed in 2.29s.
   @> Frame/model: 8
   @> Returning surface cavities
   @> Saving multiple surface cavities to directory ..
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Surface cavity calculation completed in 2.31s.
   @> Frame/model: 1
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 2425 atoms with 31833 homogeneous balls of radius 1.20 Å in 0.11s.
   @> Substituted 2425 atoms with 31833 homogeneous balls of radius 1.20 Å in 0.12s.
   @> Delaunay tessellation of 31833 points constructed in 1.35s.
   @> Delaunay tessellation of 31833 points constructed in 1.37s.
   @> Surface and inner simplices filtered in 0.39s.
   @> Surface and inner simplices filtered in 0.49s.
   @> Surface cavities: 209 found, 22 deeper than min_depth=1.5 Å and kept, in 0.23s.
   @> Returning surface cavities
   @> Saving multiple surface cavities to directory ..
   @> Surface cavity calculation completed in 2.16s.
   ..
   ..
   @> Frame/model: 47
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 2425 atoms with 31833 homogeneous balls of radius 1.20 Å in 0.10s.
   @> Delaunay tessellation of 31833 points constructed in 1.15s.
   @> Surface and inner simplices filtered in 0.40s.
   @> Surface cavities: 225 found, 25 deeper than min_depth=1.5 Å and kept, in 0.21s.
   @> Returning surface cavities
   @> Saving multiple surface cavities to directory ..
   @> Surface cavity calculation completed in 1.94s.
   @> Frame/model: 48
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 2425 atoms with 31833 homogeneous balls of radius 1.20 Å in 0.11s.
   @> Delaunay tessellation of 31833 points constructed in 1.15s.
   @> Surface and inner simplices filtered in 0.43s.
   @> Surface cavities: 251 found, 26 deeper than min_depth=1.5 Å and kept, in 0.24s.
   @> Returning surface cavities
   @> Saving multiple surface cavities to directory ..
   @> Surface cavity calculation completed in 1.99s.


The parameters of the detected surface cavities can be extracted with
:func:`.getSurfaceCavityParametersMultipleFrames`. This function analyzes the
surface cavities returned for all trajectory frames and reports, for each frame,
the volume, depth, and number of tetrahedra assigned to each detected cavity.
If ``param_file_name`` is provided, the results are also saved to text files
using this name as a prefix.

.. ipython:: python
   :verbatim:

   parameters = getSurfaceCavityParametersMultipleFrames(cavities,
   		                         param_file_name='cavi_param')

.. parsed-literal::

   @> Model/frame: 0
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	402.75 		4.87 		132
   @> cavity 1: 	392.11 		4.12 		161
   @> cavity 2: 	318.46 		6.12 		136
   @> cavity 3: 	214.1 		2.99 		110
   @> cavity 4: 	169.45 		3.64 		48
   @> cavity 5: 	116.88 		5.59 		27
   @> cavity 6: 	114.71 		1.97 		61
   @> cavity 7: 	103.98 		2.34 		47
   @> cavity 8: 	92.73 		1.84 		29
   @> cavity 9: 	89.05 		1.62 		48
   @> cavity 10: 	83.4 		1.93 		27
   @> cavity 11: 	77.01 		3.19 		60
   @> cavity 12: 	74.52 		1.99 		34
   @> cavity 13: 	67.8 		1.96 		38
   @> cavity 14: 	53.49 		2.6 		22
   @> Model/frame: 1
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	355.13 		4.06 		111
   @> cavity 1: 	273.6 		5.54 		137
   @> cavity 2: 	247.31 		2.26 		105
   @> cavity 3: 	188.24 		2.05 		78
   @> cavity 4: 	176.18 		7.41 		107
   @> cavity 5: 	162.63 		2.34 		56
   @> cavity 6: 	149.92 		3.13 		54
   @> cavity 7: 	142.73 		1.79 		61
   @> cavity 8: 	133.09 		5.38 		42
   @> cavity 9: 	94.12 		2.03 		54
   @> cavity 10: 	84.77 		1.75 		48
   @> cavity 11: 	84.44 		2.76 		26
   @> cavity 12: 	73.41 		2.48 		27
   @> cavity 13: 	53.78 		2.42 		40
   ..
   ..
   @> Model/frame: 50
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	503.78 		4.84 		145
   @> cavity 1: 	214.13 		3.85 		108
   @> cavity 2: 	188.51 		2.75 		95
   @> cavity 3: 	173.52 		2.4 		62
   @> cavity 4: 	158.52 		3.03 		71
   @> cavity 5: 	158.33 		1.82 		64
   @> cavity 6: 	130.09 		2.08 		38
   @> cavity 7: 	126.88 		2.82 		58
   @> cavity 8: 	94.47 		4.01 		51
   @> cavity 9: 	89.66 		1.78 		42
   @> cavity 10: 	86.84 		1.89 		28
   @> cavity 11: 	80.71 		2.55 		54
   @> cavity 12: 	68.29 		2.36 		22
   @> cavity 13: 	54.54 		2.17 		40
   @> cavity 14: 	54.13 		1.71 		21
   	

Residues lining the detected surface cavities can be identified with
:func:`.getSurfaceCavityResidueNamesMultipleFrames`. The function analyzes each
trajectory frame using the corresponding set of surface cavities and surface
vertices returned by :func:`.calcSurfaceCavitiesMultipleFrames`.

For each cavity, residues are selected based on their distance from the cavity
surface. The calculation uses the protein coordinates from the matching
trajectory frame, so the supplied trajectory should be the same trajectory that
was used for cavity detection. If ``residues_file_name`` is provided, the
identified residues are also saved to text files using this name as a prefix.

.. ipython:: python
   :verbatim:

   residues = getSurfaceCavityResidueNamesMultipleFrames(protein, cavities,
                             surfaces, dcd, residues_file_name='cavi_resAA')


.. parsed-literal::

   @> Surface cavity residues were saved to: cavi_resAA_frame0_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: cavi_resAA_frame1_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: cavi_resAA_frame2_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: cavi_resAA_frame3_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: cavi_resAA_frame4_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: cavi_resAA_frame5_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: cavi_resAA_frame6_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: cavi_resAA_frame7_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: cavi_resAA_frame8_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: cavi_resAA_frame9_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: cavi_resAA_frame10_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: cavi_resAA_frame11_Residues_All_surface_cavities.txt
   ..
   ..
   @> Surface cavity residues were saved to: cavi_resAA_frame43_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: cavi_resAA_frame44_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: cavi_resAA_frame45_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: cavi_resAA_frame46_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: cavi_resAA_frame47_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: cavi_resAA_frame48_Residues_All_surface_cavities.txt
   @> Surface cavity residues were saved to: cavi_resAA_frame49_Residues_All_surface_cavities.txt


The returned ``residues`` object contains residue lists for all analyzed frames.
By default, chain identifiers are included in the residue labels, which makes it
possible to distinguish residues from different protein chains.

The residue lists obtained for all trajectory frames can be further analyzed with
:func:`.calcFrequentObjectResidues`. This function counts how often individual
residues occur among the residues lining the detected objects. The counting can
be performed separately for each chain, which is useful for multichain systems
or oligomeric proteins.

In this example, the function is applied to residues lining surface cavities.
The ``object_type`` argument is used only for labeling the output and can be set
to ``'channel'``, ``'pore'``, or ``'surface_cavity'`` depending on the analyzed
object.

.. ipython:: python
   :verbatim:

   frequent_residues = calcFrequentObjectResidues(residues, 
					output_file_name='cavi_freq_res')

.. parsed-literal::

   @> Residue counts by chain were saved to: cavi_freq_res_ResCounts.txt


.. ipython:: python
   :verbatim:

   frequent_residues

.. parsed-literal::

   'P': Counter({'GLU128': 51,
             'HSE72': 51,
             'ARG75': 51,
             'LYS6': 50,
             'LEU13': 50,
             'HSE157': 50,
             'GLN76': 50,
             'THR78': 50,
             'LYS79': 50,
             'ARG40': 50,
             'PRO130': 50,
             'ASP129': 50,
             'LYS155': 49,
             'THR140': 49,
             'SER71': 49,
             'ASP137': 49,
             'TYR131': 49,
             'ARG27': 49,
             'GLU80': 48,
             'ASP42': 48,
             'TYR119': 48,
             'THR84': 48,
             'LYS102': 48,
             'ILE16': 48,
             'ILE51': 48,
             'THR31': 48,
             'TRP39': 48,
             'ALA156': 48,
             'ILE126': 47,
             'SER94': 47,
             'ILE77': 47,
             'TYR49': 47,
             'GLU154': 47,
             'LYS28': 46,
             'THR5': 46,
             'PHE85': 46,
             'LEU153': 46,
             'LYS110': 46,
             'ASP86': 46,
             'GLN124': 45,
             'VAL73': 45,
             'THR46': 45,
             'VAL106': 45,
             'GLU23': 44,
             'SER118': 44,
             'LYS123': 44,
             'ILE68': 44,
             'PRO54': 44,
             'ARG101': 44,
             'ARG150': 44,
             'VAL41': 44,
             'ALA83': 44,
             'LYS112': 44,
             'GLN60': 43,
             'GLU50': 43,
             'PRO69': 43,
             'ILE35': 43,
             'MET70': 43,
             'ASP92': 43,
             'ARG18': 43,
             'TYR87': 43,
             'GLY48': 42,
             'SER36': 42,
             'GLN33': 42,
             'ASP32': 42,
             'SER47': 41,
             'TYR57': 41,
             'SER136': 41,
             'PRO55': 41,
             'GLU37': 41,
             'LYS107': 41,
             'GLU93': 40,
             'ALA151': 40,
             'SER7': 40,
             'ASP56': 40,
             'ASP98': 40,
             'GLN143': 40,
             'THR108': 40,
             'ASN53': 39,
             'GLN105': 39,
             'ARG58': 39,
             'ARG97': 39,
             'TYR132': 39,
             'GLY14': 38,
             'VAL146': 38,
             'ASN38': 38,
             'ARG147': 37,
             'SER43': 37,
             'ASP120': 37,
             'GLY52': 36,
             'LEU29': 35,
             'GLN144': 35,
             'VAL141': 35,
             'ALA74': 34,
             'SER61': 33,
             'TYR142': 33,
             'ALA111': 33,
             'GLU114': 33,
             'ASN15': 32,
             'LEU125': 32,
             'LYS64': 32,
             'PHE82': 32,
             'GLN122': 31,
             'ASP81': 30,
             'PRO121': 30,
             ..
             ..
             'PHE26': 1})}


.. ipython:: python
   :verbatim:

   import matplotlib.pyplot as plt
   showFrequentObjectResidues(frequent_residues)
   plt.show()


.. figure:: images/cavitracer_figure37.png
   :scale: 50 %

If ``count_residue_names=True``, residue types are counted instead of individual
residue positions. For example, residues such as ``ARG75:P`` and ``ARG150:P``
are counted together as ``ARG`` for chain ``R``. This option is useful when
we want to identify which amino-acid types most frequently line the detected
surface cavities.

.. ipython:: python
   :verbatim:

   frequent_residue_names = calcFrequentObjectResidues(residues, 
						count_residue_names=True)
   frequent_residue_names

.. parsed-literal::

   {'P': Counter({'SER': 51,
             'TYR': 51,
             'ASP': 51,
             'ALA': 51,
             'LYS': 51,
             'VAL': 51,
             'GLN': 51,
             'GLU': 51,
             'GLY': 51,
             'PRO': 51,
             'ILE': 51,
             'LEU': 51,
             'THR': 51,
             'ARG': 51,
             'HSE': 51,
             'PHE': 50,
             'ASN': 50,
             'TRP': 48,
             'MET': 43,
             'CYS': 35})}

.. ipython:: python
   :verbatim:

   import matplotlib.pyplot as plt
   showFrequentObjectResidues(frequent_residue_names)
   plt.show()


.. figure:: images/cavitracer_figure38.png
   :scale: 50 %


If ``count_once_per_frame=False``, every occurrence of a residue is counted. In
this mode, if the same residue appears in more than one detected surface cavity
within the same frame, it contributes more than once to the final count. This
option emphasizes repeated object-level participation rather than simple
frame-level occurrence.

.. ipython:: python
   :verbatim:

   frequent_residue_occurrences = calcFrequentObjectResidues(
       residues,
       count_once_per_frame=False,
       output_file_name='cavi_frequent_residue_occurrences')

   showFrequentObjectResidues(frequent_residue_occurrences, top=20)
   plt.show()


.. parsed-literal::
   
   @> Residue counts by chain were saved to: cavi_frequent_residue_occurrences_ResCounts.txt


.. figure:: images/cavitracer_figure39.png
   :scale: 50 %


In general, ``count_residue_names=True`` changes what is counted, from residue
positions to residue types, whereas ``count_once_per_frame=False`` changes how
often repeated appearances are counted.


Finally, surface cavities detected across trajectory frames can be combined into
a spatial overlap map using :func:`.calcSurfaceCavityOverlaps`. In this example,
the PQR files generated by :func:`.calcSurfaceCavitiesMultipleFrames` are first
collected with :mod:`glob`. The patterns ``"cav_dcd?.pqr"`` and
``"cav_dcd??.pqr"`` select files with one- and two-digit frame numbers,
respectively.

The overlap calculation voxelizes the FIL pseudoatoms from all selected PQR
files and writes a PDB file in which the occupancy column stores the normalized
frequency of each voxel. Regions with higher occupancy correspond to surface
cavity regions that are detected in more frames.


.. ipython:: python
   :verbatim:

   import glob
   pqr_files_cavities = glob.glob("cav_dcd?.pqr") + glob.glob("cav_dcd??.pqr")
   calcSurfaceCavityOverlaps(pqr_files=pqr_files_cavities, 
			     output_file_name='surface_cavity_overlap.pdb', 
			     max_proc=4)

.. parsed-literal::

   @> Number of PQR files: 51
   @> Resolution: 0.5
   @> max_proc: 4
   @> Calculating overlaps using 4 processes.
   @> 686 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 803 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 946 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1072 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 818 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 966 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1140 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 772 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 858 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 1043 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 723 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 980 atoms and 1 coordinate sets were parsed in 0.01s.
   ..
   @> 1171 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 749 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 1020 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 962 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 691 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 936 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 899 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 869 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 885 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Overlap written to: surface_cavity_overlap.pdb
   @> Number of occupied overlap voxels: 138222

   'surface_cavity_overlap.pdb'


The results can be visualized using VMD_ as shown below.

The resulting overlap map can be visualized at different occupancy thresholds.
Below, we show examples for occupancy values greater than 0.2 and 0.5.

At the lower occupancy threshold, the surface-cavity region overlaps with the
ligand-binding site, indicating that this region is detected in part of the
trajectory. At the higher threshold, only a few recurrent cavity regions remain,
and these regions do not coincide with the ligand-binding site.

This suggests that the ligand-binding cavity may be sensitive to the selected
surface-cavity parameters in this system. However, it may also reflect the
dynamic character of this region. Since the ligand is not present in the
trajectory (and only displayed based on initial PDB structure) 
The binding-site cavity may not remain open throughout the
simulation and may instead alternate between more open and more closed
conformations.

Therefore, the default settings should be treated as a starting point rather
than a universal choice. Nevertheless, this can be further examined using
:func:`.scanSurfaceCavityParameters`, which samples different combinations of
surface-cavity parameters and helps determine whether the observed cavity is
robustly detected under alternative settings.


.. figure:: images/cavitracer_figure36.jpg
   :scale: 50 %


.. figure:: images/cavitracer_figure36b.jpg
   :scale: 50 %


.. _Trajectory tutorial: http://www.bahargroup.org/prody/tutorials/trajectory_analysis/

