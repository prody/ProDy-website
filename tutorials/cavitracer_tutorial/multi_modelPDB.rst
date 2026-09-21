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
    @> WARNING The atoms supplied to calcChannels() contain components other than protein: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein, provide an appropriate selection, for example atoms.select('protein').
    @> The structure carries its hydrogens (99% of what a complete protein would hold), so inner_radius=1.20 is more conservative than it needs to be: the 1.2 Å floor exists only to keep a sub-water probe out of the space that missing hydrogens leave open, and here that space is filled. A smaller probe, down to about 0.9 Å, measures the narrow connections instead of reporting them closed.
    @> Frame/model: 6
    @> WARNING The atoms supplied to calcChannels() contain components other than protein: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein, provide an appropriate selection, for example atoms.select('protein').
    @> The structure carries its hydrogens (99% of what a complete protein would hold), so inner_radius=1.20 is more conservative than it needs to be: the 1.2 Å floor exists only to keep a sub-water probe out of the space that missing hydrogens leave open, and here that space is filled. A smaller probe, down to about 0.9 Å, measures the narrow connections instead of reporting them closed.
    @> WARNING The atoms supplied to calcChannels() contain components other than protein: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein, provide an appropriate selection, for example atoms.select('protein').
    @> The structure carries its hydrogens (99% of what a complete protein would hold), so inner_radius=1.20 is more conservative than it needs to be: the 1.2 Å floor exists only to keep a sub-water probe out of the space that missing hydrogens leave open, and here that space is filled. A smaller probe, down to about 0.9 Å, measures the narrow connections instead of reporting them closed.
    @> WARNING The atoms supplied to calcChannels() contain components other than protein: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein, provide an appropriate selection, for example atoms.select('protein').
    @> The structure carries its hydrogens (99% of what a complete protein would hold), so inner_radius=1.20 is more conservative than it needs to be: the 1.2 Å floor exists only to keep a sub-water probe out of the space that missing hydrogens leave open, and here that space is filled. A smaller probe, down to about 0.9 Å, measures the narrow connections instead of reporting them closed.
    @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 Å in 0.17s.
    @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 Å in 0.16s.
    @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 Å in 0.16s.
    @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 Å in 0.17s.
    @> Delaunay tessellation of 48645 points constructed in 2.23s.
    @> Delaunay tessellation of 48645 points constructed in 2.25s.
    @> Delaunay tessellation of 48645 points constructed in 2.27s.
    @> Delaunay tessellation of 48645 points constructed in 2.31s.
    @> Surface and inner simplices filtered in 3.55s.
    @> Surface and inner simplices filtered in 3.65s.
    @> Surface and inner simplices filtered in 3.74s.
    @> Cavities: 124 found, 6 deeper than min_depth=5.0 Å and searched for channels, in 0.30s.
    @> Chambers (probe 1.40 Å): 6 of the 6 searched cavities have them; the other 2 are searched whole.
    @>     cavity 0: 9 chambers, 4 of them seeded.
    @>     cavity 1: 1 chamber, seeded.
    @>     cavity 2: 2 chambers, all seeded.
    @>     cavity 3: 1 chamber, none of them deep and large enough to seed; searched whole.
    @>     cavity 4: 1 chamber, seeded.
    @>     cavity 5: 1 chamber, none of them deep and large enough to seed; searched whole.
    @> 10 search sites (sp) in 0.04s: one per seeded chamber, one per cavity searched whole.
    @> Surface and inner simplices filtered in 3.89s.
    @> Channel search (Dijkstra) over 10 search sites in 6 cavities completed in 0.22s.
    @> Found 7 channels and 1 link (a link joins a deep chamber to a shallower one and never reaches the surface).
    @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
    @>     site  start_point [Å]            void                   volume [Å³]  depth [Å]  channels  links
    @>     sp0   [21.926, -10.124, 5.079]   cavity 0, chamber 1/4          522       10.9         1      1  -> sp4
    @>     sp1   [20.081, -17.989, 6.672]   cavity 0, chamber 2/4          492       14.0         -      -  sealed
    @>     sp2   [6.773, -11.822, 9.205]    cavity 3, whole                386        5.3         1      -
    @>     sp3   [7.699, -22.580, 5.815]    cavity 1, chamber 1/1          222        5.3         2      -
    @>     sp4   [21.805, -11.913, 12.999]  cavity 0, chamber 3/4          182        7.6         1      -
    @>     sp5   [24.355, -4.057, 10.820]   cavity 5, whole                165        6.2         1      -
    @>     sp6   [23.684, -13.703, -1.318]  cavity 0, chamber 4/4          135       12.1         -      -  sealed
    @>     sp7   [16.259, -26.465, -0.311]  cavity 2, chamber 1/2          102       15.2         -      -  sealed
    @>     sp8   [17.979, -31.068, -6.015]  cavity 2, chamber 2/2           86        5.3         1      -
    @>     sp9   [23.826, -27.926, -9.146]  cavity 4, chamber 1/1           77        6.5         -      -  sealed
    @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
    @> The 4 sites marked sealed above report neither a channel nor a link: no route out of them survived - either narrower than bottleneck=1.20 Å, or dropped as a duplicate of a shallower site's, or the void is its own mouth and has nowhere to path to. Lower bottleneck to see how the first kind connect.
    @> No output path given.
    @> Channel calculation completed in 6.47s.
    @> Frame/model: 5
    @> Cavities: 142 found, 8 deeper than min_depth=5.0 Å and searched for channels, in 0.31s.
    @> Cavities: 170 found, 12 deeper than min_depth=5.0 Å and searched for channels, in 0.39s.
    @> WARNING The atoms supplied to calcChannels() contain components other than protein: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein, provide an appropriate selection, for example atoms.select('protein').
    @> The structure carries its hydrogens (99% of what a complete protein would hold), so inner_radius=1.20 is more conservative than it needs to be: the 1.2 Å floor exists only to keep a sub-water probe out of the space that missing hydrogens leave open, and here that space is filled. A smaller probe, down to about 0.9 Å, measures the narrow connections instead of reporting them closed.
    @> Chambers (probe 1.40 Å): 4 of the 8 searched cavities have them; the other 5 are searched whole.
    @>     cavity 0: 11 chambers, 4 of them seeded.
    @>     cavity 1: 6 chambers, 1 of them seeded.
    @>     cavity 2: 9 chambers, 2 of them seeded.
    @>     cavity 3: 4 chambers, none of them deep and large enough to seed; searched whole.
    @> 12 search sites (sp) in 0.03s: one per seeded chamber, one per cavity searched whole.
    @> Chambers (probe 1.40 Å): 8 of the 12 searched cavities have them; the other 8 are searched whole.
    @>     cavity 0: 6 chambers, 2 of them seeded.
    @>     cavity 1: 3 chambers, 1 of them seeded.
    @>     cavity 2: 8 chambers, 2 of them seeded.
    @>     cavity 3: 6 chambers, 2 of them seeded.
    @>     cavity 4: 3 chambers, none of them deep and large enough to seed; searched whole.
    @>     cavity 6: 1 chamber, none of them deep and large enough to seed; searched whole.
    @>     cavity 7: 1 chamber, none of them deep and large enough to seed; searched whole.
    @>     cavity 10: 1 chamber, none of them deep and large enough to seed; searched whole.
    @> 15 search sites (sp) in 0.04s: one per seeded chamber, one per cavity searched whole.
    @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 Å in 0.16s.
    @> Cavities: 148 found, 7 deeper than min_depth=5.0 Å and searched for channels, in 0.30s.
    @> Channel search (Dijkstra) over 12 search sites in 8 cavities completed in 0.27s.
    @> Found 17 channels and 3 links (a link joins a deep chamber to a shallower one and never reaches the surface).
    @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
    @>     site  start_point [Å]            void                   volume [Å³]  depth [Å]  channels  links
    @>     sp0   [25.010, -17.684, 4.313]   cavity 0, chamber 1/4          601       15.2         -      1  -> sp5
    @>     sp1   [6.843, -7.297, 9.076]     cavity 3, whole                401        5.1         2      -
    @>     sp2   [14.981, -6.159, 5.257]    cavity 4, whole                230        5.2         1      -
    @>     sp3   [15.196, 2.827, 12.407]    cavity 5, whole                221        5.0         1      -
    @>     sp4   [11.799, -25.863, 1.271]   cavity 2, chamber 1/2          168       10.4         1      1  -> sp9
    @>     sp5   [24.246, -10.204, 4.146]   cavity 0, chamber 2/4          126        6.6         2      -
    @>     sp6   [12.953, -4.132, 16.569]   cavity 6, whole                122        5.1         1      -
    @>     sp7   [2.420, -11.935, 14.852]   cavity 7, whole                122        5.0         -      -  sealed
    @>     sp8   [23.812, -27.541, -6.464]  cavity 1, chamber 1/1           95        5.1         4      -
    @>     sp9   [8.380, -22.057, 5.421]    cavity 2, chamber 2/2           62        5.5         2      -
    @>     sp10  [27.160, -4.114, 11.059]   cavity 0, chamber 3/4           61        9.6         2      1  -> sp5
    @>     sp11  [15.807, -9.286, 18.768]   cavity 0, chamber 4/4           60        5.4         1      -
    @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
    @> The 1 site marked sealed above report neither a channel nor a link: no route out of them survived - either narrower than bottleneck=1.20 Å, or dropped as a duplicate of a shallower site's, or the void is its own mouth and has nowhere to path to. Lower bottleneck to see how the first kind connect.
    @> No output path given.
    @> Channel calculation completed in 6.75s.
    ..
    ..
    @> Frame/model: 19
    @> WARNING The atoms supplied to calcChannels() contain components other than protein: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein, provide an appropriate selection, for example atoms.select('protein').
    @> The structure carries its hydrogens (99% of what a complete protein would hold), so inner_radius=1.20 is more conservative than it needs to be: the 1.2 Å floor exists only to keep a sub-water probe out of the space that missing hydrogens leave open, and here that space is filled. A smaller probe, down to about 0.9 Å, measures the narrow connections instead of reporting them closed.
    @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 Å in 0.15s.
    @> Delaunay tessellation of 48645 points constructed in 1.88s.
    @> Delaunay tessellation of 48645 points constructed in 1.98s.
    @> Surface and inner simplices filtered in 2.82s.
    @> Cavities: 121 found, 5 deeper than min_depth=5.0 Å and searched for channels, in 0.31s.
    @> Surface and inner simplices filtered in 2.65s.
    @> Chambers (probe 1.40 Å): 4 of the 5 searched cavities have them; the other 2 are searched whole.
    @>     cavity 0: 11 chambers, 2 of them seeded.
    @>     cavity 1: 8 chambers, 2 of them seeded.
    @>     cavity 2: 6 chambers, 2 of them seeded.
    @>     cavity 3: 2 chambers, none of them deep and large enough to seed; searched whole.
    @> 8 search sites (sp) in 0.04s: one per seeded chamber, one per cavity searched whole.
    @> Channel search (Dijkstra) over 8 search sites in 5 cavities completed in 0.20s.
    @> Found 11 channels and 2 links (a link joins a deep chamber to a shallower one and never reaches the surface).
    @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
    @>     site  start_point [Å]            void                   volume [Å³]  depth [Å]  channels  links
    @>     sp0   [23.144, -18.236, -1.439]  cavity 0, chamber 1/2         1186       17.4         2      1  -> sp3
    @>     sp1   [23.412, -3.331, 12.911]   cavity 1, chamber 1/2          234        6.3         3      -
    @>     sp2   [22.092, -3.963, 6.241]    cavity 3, whole                230        8.3         -      -  sealed
    @>     sp3   [23.202, -26.571, -8.196]  cavity 0, chamber 2/2          185        5.1         2      -
    @>     sp4   [14.319, -7.314, 15.321]   cavity 1, chamber 2/2          124        6.8         1      -
    @>     sp5   [12.818, -26.860, -9.587]  cavity 4, whole                 98        5.4         1      -
    @>     sp6   [11.214, -25.136, 2.572]   cavity 2, chamber 1/2           90        9.2         1      -
    @>     sp7   [6.788, -21.069, 4.057]    cavity 2, chamber 2/2           69        5.2         1      1  -> sp6
    @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
    @> The 1 site marked sealed above report neither a channel nor a link: no route out of them survived - either narrower than bottleneck=1.20 Å, or dropped as a duplicate of a shallower site's, or the void is its own mouth and has nowhere to path to. Lower bottleneck to see how the first kind connect.
    @> No output path given.
    @> Channel calculation completed in 5.35s.
    @> Cavities: 158 found, 8 deeper than min_depth=5.0 Å and searched for channels, in 0.25s.
    @> Chambers (probe 1.40 Å): 5 of the 8 searched cavities have them; the other 4 are searched whole.
    @>     cavity 0: 1 chamber, seeded.
    @>     cavity 1: 2 chambers, 1 of them seeded.
    @>     cavity 2: 3 chambers, none of them deep and large enough to seed; searched whole.
    @>     cavity 3: 2 chambers, 1 of them seeded.
    @>     cavity 4: 4 chambers, 1 of them seeded.
    @> 8 search sites (sp) in 0.02s: one per seeded chamber, one per cavity searched whole.
    @> Channel search (Dijkstra) over 8 search sites in 8 cavities completed in 0.11s.
    @> Found 9 channels.
    @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
    @>     site  start_point [Å]            void                   volume [Å³]  depth [Å]  channels  links
    @>     sp0   [21.895, -2.059, 11.024]   cavity 2, whole                828        5.2         2      -
    @>     sp1   [25.921, -10.322, 5.407]   cavity 0, chamber 1/1          453        6.6         3      -
    @>     sp2   [10.964, -28.343, -1.779]  cavity 1, chamber 1/1          348       14.6         1      -
    @>     sp3   [23.261, -21.361, -7.474]  cavity 5, whole                209        5.2         -      -  sealed
    @>     sp4   [14.525, -28.962, -8.787]  cavity 6, whole                142        5.0         1      -
    @>     sp5   [13.634, -16.251, -0.681]  cavity 7, whole                114        6.5         -      -  sealed
    @>     sp6   [25.428, -27.690, -7.615]  cavity 4, chamber 1/1           72        5.0         1      -
    @>     sp7   [8.461, -13.246, 10.599]   cavity 3, chamber 1/1           55        6.6         1      -
    @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
    @> The 2 sites marked sealed above report neither a channel nor a link: no route out of them survived - either narrower than bottleneck=1.20 Å, or dropped as a duplicate of a shallower site's, or the void is its own mouth and has nowhere to path to. Lower bottleneck to see how the first kind connect.
    @> No output path given.
    @> Channel calculation completed in 5.16s.


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
    @> channel 0: 	143.02 		5.12 		2.16
    @> channel 1: 	156.56 		5.63 		2.13
    @> channel 2: 	82.11 		5.2 		1.22
    @> channel 3: 	59.58 		5.45 		1.31
    @> channel 4: 	57.34 		5.43 		1.42
    @> channel 5: 	138.17 		8.91 		1.83
    @> channel 6: 	116.31 		7.89 		1.55
    @> channel 7: 	57.24 		5.86 		1.39
    @> channel 8: 	45.48 		5.08 		1.3
    @> channel 9: 	103.51 		7.71 		1.47
    @> channel 10: 	118.17 		7.87 		1.37
    @> channel 11: 	47.46 		5.17 		1.22
    @> channel 12: 	98.62 		7.1 		1.26
    @> channel 13: 	114.23 		8.01 		1.32
    @> channel 14: 	98.95 		8.16 		1.38
    @> channel 15: 	96.32 		8.24 		1.27
    @> channel 16: 	208.28 		13.18 		1.39
    @> channel 17: 	107.31 		9.96 		1.22
    @> channel 18: 	91.08 		9.67 		1.25
    @> channel 19: 	178.54 		17.39 		1.34
    @> Frame 1
    @> channel 0: 	154.79 		7.64 		1.56
    @> channel 1: 	106.57 		5.2 		1.38
    @> channel 2: 	94.52 		6.12 		1.4
    @> channel 3: 	89.31 		7.1 		1.43
    @> channel 4: 	139.76 		8.05 		1.57
    @> channel 5: 	208.2 		9.99 		1.6
    @> channel 6: 	97.08 		9.32 		1.35
    @> channel 7: 	138.63 		14.03 		1.26
    @> channel 8: 	183.16 		18.57 		1.32
    @> Frame 2
    @> channel 0: 	172.94 		5.56 		2.04
    @> channel 1: 	100.49 		5.13 		1.55
    @> channel 2: 	81.77 		5.13 		1.49
    @> channel 3: 	140.25 		6.26 		1.65
    @> channel 4: 	196.48 		8.15 		1.82
    @> channel 5: 	149.87 		8.01 		1.85
    @> channel 6: 	56.32 		5.09 		1.31
    @> channel 7: 	57.37 		5.27 		1.38
    @> channel 8: 	178.63 		9.62 		1.87
    @> channel 9: 	71.43 		5.91 		1.24
    @> channel 10: 	83.39 		5.89 		1.22
    @> channel 11: 	88.78 		6.8 		1.42
    @> channel 12: 	188.33 		11.16 		1.49
    @> channel 13: 	124.67 		10.48 		1.52
    @> channel 14: 	159.12 		10.79 		1.28
    @> channel 15: 	192.52 		13.06 		1.28
    @> channel 16: 	88.43 		9.6 		1.25
    @> Frame 3
    @> channel 0: 	131.72 		5.05 		2.07
    @> channel 1: 	156.53 		5.91 		2.0
    @> channel 2: 	100.83 		5.11 		1.67
    @> channel 3: 	81.27 		5.26 		1.59
    @> channel 4: 	91.33 		5.72 		1.64
    @> channel 5: 	75.27 		5.82 		1.52
    @> channel 6: 	143.37 		7.55 		1.78
    @> channel 7: 	57.28 		5.61 		1.33
    @> channel 8: 	59.67 		5.15 		1.23
    @> channel 9: 	175.3 		10.95 		1.52
    @> channel 10: 	126.15 		8.6 		1.23
    @> channel 11: 	110.28 		8.7 		1.37
    @> channel 12: 	122.91 		9.33 		1.39
    @> channel 13: 	127.08 		9.82 		1.35
    @> channel 14: 	101.8 		11.01 		1.2
    @> channel 15: 	128.14 		12.99 		1.25
    @> Frame 4
    @> channel 0: 	123.07 		5.29 		1.85
    @> channel 1: 	146.08 		6.12 		1.81
    @> channel 2: 	83.73 		5.33 		1.41
    @> channel 3: 	59.39 		6.23 		1.35
    @> channel 4: 	112.83 		7.64 		1.27
    @> channel 5: 	132.86 		10.63 		1.49
    @> channel 6: 	157.45 		11.0 		1.54
    @> Frame 5
    @> channel 0: 	139.11 		6.15 		1.95
    @> channel 1: 	74.29 		5.63 		1.64
    @> channel 2: 	56.11 		5.13 		1.37
    @> channel 3: 	77.77 		6.33 		1.54
    @> channel 4: 	127.36 		7.32 		1.43
    @> channel 5: 	63.46 		5.73 		1.38
    @> channel 6: 	223.51 		12.23 		1.43
    @> channel 7: 	191.22 		13.15 		1.43
    @> channel 8: 	105.26 		10.82 		1.25
    @> Frame 6
    @> channel 0: 	153.07 		5.58 		2.05
    @> channel 1: 	120.38 		6.37 		1.8
    @> channel 2: 	48.96 		5.19 		1.24
    @> channel 3: 	57.51 		5.78 		1.35
    @> channel 4: 	62.44 		5.96 		1.23
    @> channel 5: 	53.95 		6.25 		1.21
    @> channel 6: 	185.81 		13.92 		1.54
    @> channel 7: 	139.53 		12.52 		1.49
    @> Frame 7
    @> channel 0: 	135.7 		6.14 		2.0
    @> channel 1: 	128.58 		6.44 		1.9
    @> channel 2: 	110.73 		6.82 		1.64
    @> channel 3: 	149.23 		7.41 		1.47
    @> channel 4: 	97.58 		6.56 		1.59
    @> channel 5: 	93.89 		6.93 		1.51
    @> channel 6: 	62.01 		5.05 		1.32
    @> channel 7: 	158.78 		8.78 		1.41
    @> channel 8: 	40.53 		5.18 		1.22
    @> channel 9: 	175.45 		11.14 		1.39
    @> channel 10: 	108.39 		10.77 		1.27
    @> channel 11: 	200.1 		17.9 		1.38
    @> channel 12: 	212.12 		16.48 		1.28
    @> channel 13: 	283.11 		21.41 		1.22
    @> channel 14: 	308.93 		23.28 		1.22
    ..
    ..
    @> Frame 16
    @> channel 0: 	160.79 		5.51 		2.19
    @> channel 1: 	151.01 		5.73 		2.26
    @> channel 2: 	99.34 		5.79 		1.6
    @> channel 3: 	58.35 		5.43 		1.23
    @> channel 4: 	63.46 		6.03 		1.41
    @> channel 5: 	72.01 		6.5 		1.47
    @> channel 6: 	50.92 		5.57 		1.21
    @> channel 7: 	144.15 		9.3 		1.26
    @> channel 8: 	107.2 		9.63 		1.51
    @> channel 9: 	95.96 		9.86 		1.2
    @> channel 10: 	137.89 		13.12 		1.37
    @> Frame 17
    @> channel 0: 	177.18 		5.26 		2.33
    @> channel 1: 	282.9 		7.61 		2.57
    @> channel 2: 	96.33 		5.98 		1.68
    @> channel 3: 	176.65 		6.82 		1.45
    @> channel 4: 	128.82 		7.33 		1.38
    @> channel 5: 	214.99 		9.7 		1.48
    @> channel 6: 	47.79 		5.51 		1.28
    @> channel 7: 	90.09 		9.01 		1.29
    @> channel 8: 	124.24 		8.99 		1.28
    @> channel 9: 	203.14 		19.07 		1.32
    @> channel 10: 	220.9 		20.28 		1.34
    @> Frame 18
    @> channel 0: 	142.26 		5.23 		2.06
    @> channel 1: 	177.98 		5.7 		2.27
    @> channel 2: 	92.57 		5.39 		1.37
    @> channel 3: 	179.77 		8.33 		1.48
    @> channel 4: 	124.85 		6.15 		1.52
    @> channel 5: 	111.96 		5.62 		1.44
    @> channel 6: 	89.74 		6.21 		1.59
    @> channel 7: 	60.02 		5.35 		1.31
    @> channel 8: 	53.7 		5.16 		1.23
    @> channel 9: 	207.54 		11.13 		1.8
    @> channel 10: 	58.82 		6.37 		1.36
    @> channel 11: 	115.62 		8.38 		1.35
    @> channel 12: 	206.72 		11.47 		1.59
    @> channel 13: 	110.18 		10.64 		1.26
    @> channel 14: 	218.2 		17.33 		1.35
    @> channel 15: 	235.63 		19.34 		1.38
    @> Frame 19
    @> channel 0: 	107.73 		5.05 		1.81
    @> channel 1: 	136.68 		6.41 		1.62
    @> channel 2: 	54.64 		5.22 		1.31
    @> channel 3: 	60.05 		5.38 		1.27
    @> channel 4: 	110.02 		6.73 		1.31
    @> channel 5: 	146.61 		9.13 		1.55
    @> channel 6: 	69.13 		6.7 		1.32
    @> channel 7: 	181.61 		12.24 		1.57
    @> channel 8: 	207.35 		15.95 		1.23
    [([5.1159961833333965,
       5.627233267369147,
       5.195727657739742,
       5.448758815170745,
       5.4342705574128685,
       8.914863774483388,
       7.893319750827118,
       5.858472451134283,
       5.075313425942609,
       7.712208761869687,
       7.867250085408082,
       5.173823112220275,
       7.10093008540346,
       8.008011053455906,
       8.156038475279328,
       8.235104378712549,
       13.177331465945684,
       9.96446376274854,
       9.666753290902463,
       17.38806847944559],
      [2.1567777650960682,
       2.1263460795719897,
       1.2225037003157293,
       1.3146182253324865,
       1.4208793751816609,
       1.8338274580135379,
       1.5454877356054963,
       1.3932112907545389,
       1.300583861012459,
       1.4713319711686845,
       1.3690313673904264,
       1.2168517054007506,
       1.264148408689324,
       1.3203717383382563,
       1.378929287436679,
       1.2717473585146666,
       1.3905361568756553,
       1.2187824258754498,
       1.248574734574169,
       1.342019728498342],
      [143.01938089982764,
       156.55949552547258,
       82.10540362090218,
       59.58486200649433,
       57.33920986457909,
       138.16748738383325,
       116.30741093339563,
       57.24053264558316,
       45.47728598631714,
       103.5060252778599,
       118.16564640566253,
       47.45974379153111,
       98.62223929044059,
       114.22905984975547,
       98.95348718793016,
       96.3235430857895,
       208.27967858539594,
       107.30974859958471,
       91.07624460304994,
       178.54472441659416]),
     ([7.637175533932657,
       5.198226273241115,
       6.1202562529259374,
       7.09785151917118,
       8.045414340792902,
       9.988857932075868,
       9.3248161421848,
       14.033065550178495,
       18.567802015599554],
      [1.5601494686064619,
       1.3767172331905473,
       1.396267082044916,
       1.4284237637831512,
       1.5658677317956504,
       1.6029776227073833,
       1.350213699846994,
       1.2608484709557646,
       1.3209182777512203],
      [154.78550417343658,
       106.57367513483952,
       94.52136455364138,
       89.31387928861822,
       139.75898360524664,
       208.19582161514597,
       97.07859977306268,
       138.62786733595408,
       183.16315631016775]),
     ([5.564768586770382,
       5.125592837642575,
       5.132347147735286,
       6.257068432764722,
       8.149521134951737,
       8.007123358414487,
       5.086103516156363,
       5.274594489689208,
       9.623714975261752,
       5.9139518312774975,
       5.889984568056052,
       6.801042709097185,
       11.15775350525378,
       10.482390989525939,
       10.792289800886987,
       13.059968229382683,
       9.604060489646617],
      [2.0436517377300296,
       1.5547935046061492,
       1.4891256704445446,
       1.6535939598156253,
       1.8151141150300918,
       1.8511633456153114,
       1.309977556279841,
       1.3800434532746404,
       1.870197280118645,
       1.2432104412467455,
       1.2226089559395221,
       1.4214477666594607,
       1.493602183399627,
       1.516331471359109,
       1.2815724347258317,
       1.2815724347258317,
       1.2507665880725696],
      [172.93678656814484,
       100.48917270529004,
       81.77402977045877,
       140.2521178144088,
       196.4770027439564,
       149.8675814687977,
       56.316067985357,
       57.365740153520655,
       178.6271622482581,
       71.43413129722485,
       83.38671728390946,
       88.78179468487173,
       188.33130134861608,
       124.6717346137086,
       159.1176396149233,
       192.51539363268603,
       88.43270566437143]),
       ..
       ..
     ([5.258814934303968,
       7.605733776094616,
       5.975818122746541,
       6.821570913596652,
       7.330579344866887,
       9.698950245933322,
       5.512250959836164,
       9.012343912361064,
       8.987751594149437,
       19.07062744142514,
       20.283693624695125],
      [2.3300009054509756,
       2.5715037905131197,
       1.684856898326969,
       1.4519718419488827,
       1.3792787483511126,
       1.4820406846925873,
       1.2751981481515553,
       1.2905790752200803,
       1.279487485982471,
       1.317007322252801,
       1.343393256412998],
      [177.1843597018072,
       282.8962476967863,
       96.32607701073246,
       176.65024211765223,
       128.81603141393205,
       214.98867763107393,
       47.79156372743498,
       90.08685897110804,
       124.2389696171072,
       203.1362470767652,
       220.89605848796035]),
     ([5.23091456960509,
       5.701560589910243,
       5.390283488981547,
       8.332532345767264,
       6.145972348252318,
       5.623731361463549,
       6.210711121115695,
       5.34630528664903,
       5.161388981141756,
       11.125681179624504,
       6.367069016340153,
       8.37610256207036,
       11.470239296338878,
       10.635889912461767,
       17.33289175108473,
       19.343934495380523],
      [2.0648574448259875,
       2.273908917574122,
       1.3656257482244853,
       1.480415060441138,
       1.5220880820135527,
       1.4394601747258473,
       1.588968041403174,
       1.3057836380450885,
       1.2265880296810254,
       1.7987375406554766,
       1.3587098133472473,
       1.3496821937207135,
       1.5939716505163937,
       1.2565407530514239,
       1.3476175878099073,
       1.3797681454934472],
      [142.25902644518777,
       177.98229183317824,
       92.57336596469592,
       179.76562081830508,
       124.84832334058711,
       111.96219455353156,
       89.74393175483978,
       60.019829973979256,
       53.698656532431556,
       207.54448810603398,
       58.82081828828874,
       115.62351397605242,
       206.7223729193255,
       110.18399869120566,
       218.20249936375973,
       235.62984831530122]),
     ([5.045736488701314,
       6.409882870352796,
       5.215578946984696,
       5.37842948280604,
       6.726471648000702,
       9.13100622731233,
       6.6982197682861795,
       12.243694703104588,
       15.952173832088667],
      [1.8070795784810165,
       1.616931848718034,
       1.311088541070305,
       1.2660582399424871,
       1.3116047511790423,
       1.5545485726276655,
       1.3245688120948775,
       1.5680374935807395,
       1.2315669424777755],
      [107.73182750342251,
       136.67619529707667,
       54.63779846892617,
       60.048378838556516,
       110.01672587876756,
       146.61100323224144,
       69.13194151204993,
       181.61442134605466,
       207.34898023929765])]


Access to the residues that are forming the channels is provided by
:func:`.getChannelResidueNames` function. Below, the example on how to
obtain information for frame #0.

.. ipython:: python
   :verbatim:

   p3.setACSIndex(0)
   getChannelResidueNames(p3, channels3[0])

.. parsed-literal::

    ['channel0: PHE137:A, MET140:A, GLY141:A, GLY144:A, ALA147:A, ALA148:A, ALA151:A, TYR200:A, TYR204:A, LEU219:A, TYR223:A, RET301:A',
     'channel1: SER43:A, THR44:A, PHE47:A, GLY66:A, THR69:A, GLY70:A, PHE73:A, LYS231:A, ILE232:A, GLY235:A, LEU236:A',
     'channel2: SER32:A, LEU35:A, VAL36:A, ALA39:A, LEU221:A, ASN224:A, LEU225:A, PHE228:A',
     'channel3: GLY196:A, TRP197:A, ILE199:A, TYR200:A, LEU225:A, ALA226:A, VAL229:A, ASN230:A',
     'channel4: ALA148:A, TRP149:A, PHE152:A, PRO201:A, TYR204:A, PHE205:A, TYR208:A',
     'channel5: TYR110:A, LEU113:A, ALA114:A, ALA116:A, THR117:A, ALA120:A, GLY171:A, ALA174:A, ALA178:A, SER179:A, VAL182:A, TYR186:A',
     'channel6: TYR110:A, LEU113:A, ALA114:A, ALA120:A, LEU166:A, TRP167:A, ALA168:A, GLY169:A, GLU170:A, GLY171:A, LYS172:A, VAL182:A, TYR186:A',
     'channel7: GLY21:A, GLY207:A, MET210:A, GLY211:A, ASP212:A, GLY214:A, SER215:A, ASN218:A, LEU219:A',
     'channel8: GLY207:A, TYR208:A, LEU209:A, MET210:A, GLY211:A, ASP212:A, GLY213:A, GLY214:A, SER215:A, ASN218:A, LEU219:A',
     'channel9: TYR110:A, LEU113:A, ALA114:A, THR117:A, VAL119:A, ALA120:A, SER122:A, GLU170:A, GLY171:A, ALA174:A, VAL182:A, TYR186:A',
     'channel10: VAL133:A, PHE152:A, GLY155:A, CYS156:A, ALA158:A, TRP159:A, TRP197:A, ALA198:A, PRO201:A, RET301:A',
     'channel11: TYR28:A, THR29:A, SER32:A, ILE84:A, LEU217:A, ASN220:A, LEU221:A',
     'channel12: LYS57:A, TRP58:A, SER61:A, PHE109:A, ILE112:A, LEU113:A, ALA116:A, VAL182:A',
     'channel13: THR44:A, PHE47:A, LEU62:A, SER65:A, GLY66:A, THR69:A, GLY70:A, PHE73:A, LYS231:A, PHE234:A, GLY235:A, ILE238:A',
     'channel14: TYR110:A, LEU113:A, ALA114:A, ALA120:A, GLY121:A, SER122:A, LYS125:A, GLY169:A, GLU170:A, GLY171:A, VAL182:A, TYR186:A',
     'channel15: LEU132:A, VAL133:A, VAL136:A, ILE154:A, GLY155:A, ALA158:A, TRP159:A, TRP197:A, RET301:A',
     'channel16: TYR76:A, ARG80:A, TRP83:A, THR91:A, ARG94:A, TRP98:A, PHE137:A, GLY141:A, GLY144:A, TYR200:A, TYR223:A, ALA226:A, ASP227:A, ASN230:A, RET301:A',
     'channel17: TYR110:A, LEU113:A, ALA114:A, ALA120:A, GLY171:A, LYS172:A, ALA178:A, VAL182:A, GLN183:A, TYR186:A',
     'channel18: LEU105:A, ILE106:A, MET189:A, ILE192:A, ILE193:A, PHE195:A, GLY196:A, TRP197:A, ILE199:A, VAL229:A, ASN230:A, LEU233:A, PHE234:A, RET301:A',
     'channel19: TRP98:A, LEU99:A, VAL102:A, PRO103:A, LEU105:A, ILE106:A, LYS126:A, LEU127:A, GLY130:A, SER131:A, VAL133:A, MET134:A, MET189:A, ILE193:A, TRP197:A, ASN230:A, PHE234:A, RET301:A']


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

    @> Channel residues were saved to: results_frame1_Residues_All_channels.txt
    ['channel0: L40:A, S43:A, T44:A, F47:A, G66:A, T69:A, G70:A, F73:A, K231:A, I232:A, G235:A, L236:A',
     'channel1: Y110:A, L113:A, A114:A, A120:A, G121:A, S122:A, G169:A, E170:A, G171:A, A174:A, V182:A, Y186:A',
     'channel2: A151:A, F152:A, G155:A, C156:A, W159:A, W197:A, A198:A, P201:A, V202:A, RET301:A',
     'channel3: S43:A, F46:A, F47:A, G66:A, T69:A, K231:A, I232:A, G235:A, L236:A',
     'channel4: Y110:A, L113:A, A114:A, T117:A, A120:A, E170:A, G171:A, S173:A, A174:A, T177:A, A178:A, V182:A, Y186:A',
     'channel5: W83:A, V136:A, F137:A, M140:A, G141:A, G144:A, M146:A, A147:A, P150:A, I154:A, Y204:A, Y208:A, L219:A, Y223:A, RET301:A',
     'channel6: W98:A, L99:A, V102:A, P103:A, L127:A, G130:A, S131:A, V133:A, M134:A, W197:A, RET301:A',
     'channel7: F47:A, K57:A, W58:A, S61:A, L62:A, S65:A, G66:A, F109:A, I112:A, L113:A, A115:A, A116:A, V182:A, F234:A, G235:A, I238:A, W239:A',
     'channel8: Y76:A, T91:A, V92:A, R94:A, Y95:A, W98:A, G130:A, V133:A, M134:A, F137:A, G138:A, E142:A, W197:A, RET301:A']


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

To identify pores, we should use :func:`.calcChannelsMultipleFrames`
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
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Frame/model: 4
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Frame/model: 6
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 Å in 0.29s.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 Å in 0.30s.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 Å in 0.29s.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 Å in 0.30s.
   @> Delaunay tessellation of 77434 points constructed in 3.91s.
   @> Delaunay tessellation of 77434 points constructed in 3.94s.
   @> Delaunay tessellation of 77434 points constructed in 3.97s.
   @> Delaunay tessellation of 77434 points constructed in 4.04s.
   @> Surface and inner simplices filtered in 3.82s.
   @> Surface and inner simplices filtered in 3.74s.
   @> Surface and inner simplices filtered in 3.74s.
   @> Surface and inner simplices filtered in 3.95s.
   @> Surface cavities: 380 found, 11 deeper than min_depth=5.0 Å and searched for channels, in 0.61s.
   @> Surface cavities: 396 found, 9 deeper than min_depth=5.0 Å and searched for channels, in 0.66s.
   @> Surface cavities: 384 found, 9 deeper than min_depth=5.0 Å and searched for channels, in 0.63s.
   @> Surface cavities: 386 found, 8 deeper than min_depth=5.0 Å and searched for channels, in 0.60s.
   @> Chambers (probe 1.40 Å): 2 of the 11 searched cavities have them; the other 10 are searched whole.
   @>     cavity 0: 20 chambers, 8 of them seeded.
   @>     cavity 5: 1 chamber, none of them deep and large enough to seed; searched whole.
   @> 18 search sites (sp) in 0.05s: one per seeded chamber, one per cavity searched whole.
   @> Chambers (probe 1.40 Å): 3 of the 9 searched cavities have them; the other 7 are searched whole.
   @>     cavity 0: 18 chambers, 9 of them seeded.
   @>     cavity 1: 3 chambers, 1 of them seeded.
   @>     cavity 3: 1 chamber, none of them deep and large enough to seed; searched whole.
   @> 17 search sites (sp) in 0.09s: one per seeded chamber, one per cavity searched whole.
   @> Chambers (probe 1.40 Å): 4 of the 9 searched cavities have them; the other 6 are searched whole.
   @>     cavity 0: 12 chambers, 5 of them seeded.
   @>     cavity 1: 2 chambers, all seeded.
   @>     cavity 2: 3 chambers, none of them deep and large enough to seed; searched whole.
   @>     cavity 4: 1 chamber, seeded.
   @> 14 search sites (sp) in 0.09s: one per seeded chamber, one per cavity searched whole.
   @> Chambers (probe 1.40 Å): 4 of the 8 searched cavities have them; the other 5 are searched whole.
   @>     cavity 0: 23 chambers, 11 of them seeded.
   @>     cavity 1: 2 chambers, 1 of them seeded.
   @>     cavity 2: 1 chamber, none of them deep and large enough to seed; searched whole.
   @>     cavity 3: 1 chamber, seeded.
   @> 18 search sites (sp) in 0.06s: one per seeded chamber, one per cavity searched whole.
   @> Channel search (Dijkstra) over 14 search sites in 9 cavities completed in 1.40s.
   @> Found 35 channels and 6 links (a link joins a deep chamber to a shallower one and never reaches the surface).
   @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
   @>     site  void                   volume [Å³]  depth [Å]  channels  links
   @>     sp0   cavity 0, chamber 1/5         2747       16.5        14      -
   @>     sp1   cavity 2, whole                680        5.3         6      -
   @>     sp2   cavity 3, whole                645        5.1         1      -
   @>     sp3   cavity 5, whole                236       10.4         1      -
   @>     sp4   cavity 0, chamber 2/5          220        5.3         3      2  -> sp0, sp0
   @>     sp5   cavity 0, chamber 3/5          151       34.6         -      1  -> sp6
   @>     sp6   cavity 0, chamber 4/5          111       26.0         -      1  -> sp0
   @>     sp7   cavity 6, whole                 92        5.1         1      -
   @>     sp8   cavity 7, whole                 73        5.2         -      -  sealed
   @>     sp9   cavity 8, whole                 57        5.4         1      -
   @>     sp10  cavity 1, chamber 1/2           42        6.8         3      1  -> sp11
   @>     sp11  cavity 1, chamber 2/2           41        5.1         3      -
   @>     sp12  cavity 4, chamber 1/1           39        6.5         1      -
   @>     sp13  cavity 0, chamber 5/5           34        7.4         1      1  -> sp0
   @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
   @> Saving 35 channels and 6 links to directory ., one file per object named sp<site>_chl<n> and sp<site>_lnk<n>.
   @> Channel calculation completed in 10.64s.
   ..
   @> Found 25 channels and 6 links (a link joins a deep chamber to a shallower one and never reaches the surface).
   @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
   @>     site  void                   volume [Å³]  depth [Å]  channels  links
   @>     sp0   cavity 0, chamber 1/7         2608       12.4        10      -
   @>     sp1   cavity 0, chamber 2/7          231       13.7         2      1  -> sp0
   @>     sp2   cavity 1, whole                182        6.3         1      -
   @>     sp3   cavity 0, chamber 3/7          138        5.5         6      1  -> sp0
   @>     sp4   cavity 2, whole                 83        5.1         -      -  sealed
   @>     sp5   cavity 3, whole                 78        5.1         1      -
   @>     sp6   cavity 4, whole                 74        5.1         -      -  sealed
   @>     sp7   cavity 0, chamber 4/7           54       10.8         -      2  -> sp1, sp0
   @>     sp8   cavity 0, chamber 5/7           48        5.8         2      1  -> sp0
   @>     sp9   cavity 0, chamber 6/7           44       11.4         -      1  -> sp8
   @>     sp10  cavity 0, chamber 7/7           31        5.1         3      -
   @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
   @> Saving 25 channels and 6 links to directory ., one file per object named sp<site>_chl<n> and sp<site>_lnk<n>.
   @> Channel calculation completed in 8.90s.
   @> Channel search (Dijkstra) over 18 search sites in 7 cavities completed in 1.75s.
   @> Found 33 channels and 10 links (a link joins a deep chamber to a shallower one and never reaches the surface).
   @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
   @>     site  void                     volume [Å³]  depth [Å]  channels  links
   @>     sp0   cavity 0, chamber 1/12          2871       15.8        15      1  -> sp3
   @>     sp1   cavity 1, whole                  350        5.2         1      -
   @>     sp2   cavity 2, whole                  334       10.5         -      -  sealed
   @>     sp3   cavity 0, chamber 2/12           315       10.8         5      -
   @>     sp4   cavity 0, chamber 3/12           293       12.3         1      2  -> sp3, sp0
   @>     sp5   cavity 3, whole                  269        9.4         -      -  sealed
   @>     sp6   cavity 4, whole                  267        5.4         2      -
   @>     sp7   cavity 5, whole                  148        9.0         -      -  sealed
   @>     sp8   cavity 6, whole                   84        5.3         -      -  sealed
   @>     sp9   cavity 0, chamber 4/12            57        9.4         -      -  sealed
   @>     sp10  cavity 0, chamber 5/12            56       12.6         -      3  -> sp9, sp4, sp0
   @>     sp11  cavity 0, chamber 6/12            53       24.2         -      -  sealed
   @>     sp12  cavity 0, chamber 7/12            51        5.1         5      -
   @>     sp13  cavity 0, chamber 8/12            45       20.3         -      1  -> sp0
   @>     sp14  cavity 0, chamber 9/12            42        5.5         1      1  -> sp0
   @>     sp15  cavity 0, chamber 10/12           40        9.8         1      -
   @>     sp16  cavity 0, chamber 11/12           31        5.0         2      -
   @>     sp17  cavity 0, chamber 12/12           31       18.3         -      2  -> sp0, sp14
   @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
   @> The 2 sites marked sealed above report neither a channel nor a link: every route out of them is narrower than bottleneck=0.85 Å. Lower it to see how they connect.
   @> Saving 33 channels and 10 links to directory ., one file per object named sp<site>_chl<n> and sp<site>_lnk<n>.
   @> Channel calculation completed in 9.47s.


After channels, :func:`.calcPoresFromChannelsMultipleFrames` function can be
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
:func:`.getPoreResidueNamesMultipleFrames` function. When applying
``one_letter_aa=True`` all the residues will be written in a one-letter
code.

.. ipython:: python
   :verbatim:

   getPoreResidueNamesMultipleFrames(pdb_multi, pores, one_letter_aa=True)

.. parsed-literal::

   @> Model: 0
   @> Model: 1
   @> Model: 2
   @> Model: 3
   @> Model: 4
   @> Model: 5
   @> Model: 6
   @> Model: 7
   @> Model: 8
   @> Model: 9
   @> Model: 10
   @> Model: 11
   @> Model: 12
   @> Model: 13
   @> Model: 14
   @> Model: 15
   @> Model: 16
   @> Model: 17
   @> Model: 18
   @> Model: 19
   @> Model: 20

   [[],
    [],
    [],
    [],
    [],
    [],
    [],
    ['pore0: I22:P, I25:P, V26:P, A29:P, L30:P, D33:P, N34:P, L37:P, T38:P, V41:P, I44:P, P45:P, S46:P, S119:P, E120:P, L124:P, N128:P, V131:P, F135:P, K138:P, Q142:P, Y158:P, P159:P, I162:P, S196:P, S199:P, S200:P, G203:P, M206:P, L228:P, V232:P, L270:P, I308:P, L311:P, E312:P, I317:P, W318:P, M319:P, F334:P, A337:P, S338:P, Y341:P, N388:P, G392:P, I395:P, D426:P, F429:P, Y433:P',
     'pore1: I22:P, I25:P, V26:P, A29:P, L30:P, D33:P, N34:P, L37:P, T38:P, V41:P, I44:P, P45:P, S46:P, E120:P, L124:P, N128:P, V131:P, F135:P, K138:P, Q142:P, Y158:P, P159:P, I162:P, S196:P, S199:P, S200:P, G203:P, M206:P, L228:P, V232:P, L270:P, I308:P, L311:P, E312:P, I317:P, W318:P, M319:P, F334:P, A337:P, S338:P, Y341:P, N388:P, G392:P, I395:P, D426:P, F429:P, Y433:P',
     'pore2: I22:P, I25:P, V26:P, A29:P, L30:P, D33:P, N34:P, L37:P, T38:P, V41:P, I44:P, P45:P, S46:P, S119:P, E120:P, L124:P, N128:P, V131:P, F135:P, K138:P, Q142:P, Y158:P, P159:P, I162:P, S196:P, S199:P, S200:P, G203:P, M206:P, L228:P, V232:P, L270:P, I308:P, L311:P, E312:P, I317:P, W318:P, F334:P, A337:P, S338:P, Y341:P, N388:P, G392:P, I395:P, D426:P, F429:P, Y433:P'],
    [],
    [],
    ['pore0: I22:P, I25:P, V26:P, A29:P, L30:P, D33:P, N34:P, L37:P, T38:P, V40:P, V41:P, I44:P, K122:P, D123:P, L124:P, N128:P, V131:P, G132:P, F135:P, K138:P, Q142:P, I162:P, F166:P, M169:P, S196:P, S199:P, S200:P, G203:P, M206:P, L228:P, V232:P, D262:P, Q266:P, V269:P, L270:P, Q276:P, I308:P, E312:P, P313:P, A314:P, L315:P, P316:P, I317:P, W318:P, F334:P, Y341:P, K379:P, D426:P, F429:P, Y433:P, P437:P',
     'pore1: I22:P, I25:P, V26:P, A29:P, L30:P, D33:P, N34:P, L37:P, T38:P, V40:P, V41:P, I44:P, L124:P, N128:P, V131:P, G132:P, F135:P, K138:P, Q142:P, I162:P, F166:P, M169:P, S196:P, S199:P, S200:P, G203:P, M206:P, L228:P, V232:P, D262:P, Q266:P, V269:P, L270:P, Q276:P, I308:P, E312:P, P313:P, L315:P, P316:P, I317:P, W318:P, F334:P, Y341:P, K379:P, D426:P, F429:P, Y433:P, P437:P',
     'pore2: I22:P, I25:P, V26:P, A29:P, L30:P, D33:P, N34:P, L37:P, T38:P, V40:P, V41:P, I44:P, S46:P, V131:P, F135:P, K138:P, Q142:P, I162:P, F166:P, M169:P, S196:P, S199:P, S200:P, G203:P, M206:P, L228:P, V232:P, D262:P, Q266:P, V269:P, L270:P, Q276:P, I308:P, E312:P, I317:P, W318:P, M319:P, E321:P, T322:P, R326:P, K327:P, W328:P, Q329:P, L330:P, F334:P, Y341:P, D426:P, F429:P, Y433:P'],
    [],
    [],
    [],
    [],
    ['pore0: L30:P, N34:P, L37:P, T38:P, V41:P, D121:P, K122:P, D123:P, L124:P, E127:P, N128:P, V131:P, G132:P, F135:P, K138:P, Q142:P, G150:P, T153:P, N154:P, S200:P, V201:P, A202:P, M204:P, G205:P, L225:P, L228:P, V232:P, S279:P, Q280:P, K281:P, G282:P, E312:P, P313:P, A314:P, L315:P, P316:P, I317:P, W318:P, F334:P, M403:P, V417:P, G419:P, S420:P, Y422:P, A423:P, A425:P, D426:P, F429:P, C430:P, Y433:P',
     'pore1: L30:P, N34:P, L37:P, T38:P, V41:P, S119:P, E120:P, D121:P, K122:P, D123:P, L124:P, E127:P, N128:P, V131:P, G132:P, F135:P, K138:P, Q142:P, G150:P, T153:P, N154:P, S200:P, V201:P, A202:P, M204:P, G205:P, L225:P, L228:P, V232:P, Q280:P, G282:P, T283:P, E312:P, P313:P, A314:P, L315:P, P316:P, I317:P, W318:P, M319:P, F334:P, M403:P, V417:P, G419:P, S420:P, Y422:P, A423:P, A425:P, D426:P, F429:P, C430:P, Y433:P',
     'pore2: L30:P, N34:P, L37:P, T38:P, V41:P, S119:P, E120:P, D121:P, K122:P, D123:P, L124:P, E127:P, N128:P, V131:P, G132:P, F135:P, K138:P, Q142:P, G150:P, T153:P, N154:P, S200:P, V201:P, A202:P, M204:P, G205:P, L225:P, L228:P, V232:P, S279:P, Q280:P, K281:P, G282:P, E312:P, P313:P, A314:P, L315:P, P316:P, I317:P, W318:P, M319:P, F334:P, M403:P, V417:P, G419:P, S420:P, Y422:P, A423:P, A425:P, D426:P, F429:P, C430:P, Y433:P',
     'pore3: L30:P, N34:P, L37:P, T38:P, V41:P, D121:P, K122:P, D123:P, L124:P, E127:P, N128:P, V131:P, G132:P, F135:P, K138:P, Q142:P, G150:P, T153:P, N154:P, S200:P, V201:P, A202:P, M204:P, G205:P, L225:P, L228:P, V232:P, S279:P, Q280:P, K281:P, G282:P, E312:P, P313:P, A314:P, L315:P, P316:P, I317:P, W318:P, M319:P, F334:P, M403:P, V417:P, G419:P, S420:P, Y422:P, A423:P, A425:P, D426:P, F429:P, C430:P, Y433:P',
     'pore4: L30:P, N34:P, L37:P, T38:P, V41:P, E120:P, D121:P, K122:P, D123:P, L124:P, E127:P, N128:P, V131:P, G132:P, F135:P, K138:P, Q142:P, G150:P, T153:P, N154:P, S200:P, V201:P, A202:P, M204:P, G205:P, L225:P, L228:P, V232:P, S279:P, Q280:P, K281:P, G282:P, E312:P, P313:P, I317:P, W318:P, F334:P, M403:P, V417:P, G419:P, S420:P, Y422:P, A423:P, A425:P, D426:P, F429:P, C430:P, Y433:P'],
    [],
    [],
    [],
    [],
    []]


To obtain information about pore's parameters, such as volume, length, or
bottlenck, use :func:`.getPoreParametersMultipleFrames` function.

.. ipython:: python
   :verbatim:

   getPoreParametersMultipleFrames(pores)


.. parsed-literal::

   @> Frame/model: 0
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 1
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
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
   @> pore 0: 	1176.59 		70.54 		0.91
   @> pore 1: 	1186.28 		72.98 		0.91
   @> pore 2: 	1167.48 		72.58 		0.89
   @> Frame/model: 8
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 9
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 10
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	1238.29 		69.38 		0.99
   @> pore 1: 	1213.77 		69.32 		1.05
   @> pore 2: 	1154.92 		72.43 		0.9
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
   @> pore 0: 	850.56 		66.88 		0.88
   @> pore 1: 	846.51 		65.5 		0.89
   @> pore 2: 	843.03 		68.95 		0.88
   @> pore 3: 	823.02 		69.18 		0.88
   @> pore 4: 	805.46 		69.65 		0.88
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

   [([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([70.54471000644224, 72.97668159155492, 72.57642405027067],
     [0.9118922696158223, 0.9118922696158223, 0.8935285215508568],
     [1176.5915219465767, 1186.2812882068201, 1167.4824441528995]),
    ([], [], []),
    ([], [], []),
    ([69.38102226316408, 69.32201029179326, 72.42687193686625],
     [0.9935284974640439, 1.052007326607179, 0.8999179901412739],
     [1238.2922478727692, 1213.7679963582323, 1154.9207891309438]),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([66.88389586204593,
      65.50350047091752,
      68.94966156122996,
      69.17585872830173,
      69.65214108106966],
     [0.8785444091924826,
      0.8880679059215846,
      0.8785444091924826,
      0.8785444091924826,
      0.8785444091924826],
     [850.5576852359668,
      846.5127103297621,
      843.0339610299748,
      823.0212056231321,
      805.4572654140516]),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], []),
    ([], [], [])]


