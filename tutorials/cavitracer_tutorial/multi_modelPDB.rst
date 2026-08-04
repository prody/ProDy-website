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
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 A in 0.19s.
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 A in 0.18s.
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 A in 0.18s.
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 A in 0.18s.
   @> Delaunay tessellation of 48645 points constructed in 2.36s.
   @> Delaunay tessellation of 48645 points constructed in 2.35s.
   @> Delaunay tessellation of 48645 points constructed in 2.39s.
   @> Delaunay tessellation of 48645 points constructed in 2.39s.
   @> Surface and inner simplices filtered in 3.00s.
   @> Surface and inner simplices filtered in 3.00s.
   @> Surface and inner simplices filtered in 3.02s.
   @> Surface and inner simplices filtered in 3.08s.
   @> 2 surface cavities detected and filtered in 0.38s.
   @> 3 surface cavities detected and filtered in 0.38s.
   @> 4 surface cavities detected and filtered in 0.37s.
   @> 2 surface cavities detected and filtered in 0.42s.
   @> Channel pathfinding (graph Dijkstra) over 4 cavities completed in 1.65s.
   @> Detected 21 channels.
   @> No output path given.
   @> Channel calculation completed in 7.65s.
   @> Frame/model: 7
   @> WARNING The atoms supplied to calcChannels() contain non-protein components: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein structure, provide an appropriate selection, for example atoms.select('protein').
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 A in 0.16s.
   @> Channel pathfinding (graph Dijkstra) over 3 cavities completed in 2.28s.
   @> Detected 40 channels.
   @> No output path given.
   @> Channel calculation completed in 8.24s.
   @> Frame/model: 5
   @> WARNING The atoms supplied to calcChannels() contain non-protein components: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein structure, provide an appropriate selection, for example atoms.select('protein').
   @> Channel pathfinding (graph Dijkstra) over 2 cavities completed in 2.45s.
   @> Detected 32 channels.
   @> No output path given.
   @> Channel calculation completed in 8.38s.
   @> Frame/model: 3
   @> WARNING The atoms supplied to calcChannels() contain non-protein components: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein structure, provide an appropriate selection, for example atoms.select('protein').
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 A in 0.17s.
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 A in 0.17s.
   @> Delaunay tessellation of 48645 points constructed in 1.92s.
   @> Delaunay tessellation of 48645 points constructed in 2.04s.
   @> Delaunay tessellation of 48645 points constructed in 2.09s.
   @> Channel pathfinding (graph Dijkstra) over 2 cavities completed in 4.81s.
   @> Detected 39 channels.
   @> No output path given.
   @> Channel calculation completed in 10.84s.
   @> Frame/model: 1
   @> WARNING The atoms supplied to calcChannels() contain non-protein components: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein structure, provide an appropriate selection, for example atoms.select('protein').
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 A in 0.19s.
   @> Surface and inner simplices filtered in 2.79s.
   @> 4 surface cavities detected and filtered in 0.45s.
   @> Surface and inner simplices filtered in 2.68s.
   @> Delaunay tessellation of 48645 points constructed in 2.25s.
   @> Surface and inner simplices filtered in 2.75s.
   @> 4 surface cavities detected and filtered in 0.40s.
   @> 5 surface cavities detected and filtered in 0.39s.
   @> Channel pathfinding (graph Dijkstra) over 4 cavities completed in 1.86s.
   @> Detected 40 channels.
   @> No output path given.
   @> Channel calculation completed in 7.19s.
   @> Frame/model: 8
   @> WARNING The atoms supplied to calcChannels() contain non-protein components: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein structure, provide an appropriate selection, for example atoms.select('protein').
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 A in 0.19s.
   @> Surface and inner simplices filtered in 2.23s.
   @> Channel pathfinding (graph Dijkstra) over 4 cavities completed in 2.04s.
   @> Detected 41 channels.
   @> No output path given.
   @> Channel calculation completed in 7.35s.
   @> Frame/model: 10
   @> WARNING The atoms supplied to calcChannels() contain non-protein components: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein structure, provide an appropriate selection, for example atoms.select('protein').
   @> Channel pathfinding (graph Dijkstra) over 5 cavities completed in 1.98s.
   @> Detected 33 channels.
   @> No output path given.
   @> Channel calculation completed in 7.40s.
   @> Frame/model: 12
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 A in 0.18s.
   @> WARNING The atoms supplied to calcChannels() contain non-protein components: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein structure, provide an appropriate selection, for example atoms.select('protein').
   @> 4 surface cavities detected and filtered in 0.38s.
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 A in 0.18s.
   @> Delaunay tessellation of 48645 points constructed in 1.96s.
   @> Channel pathfinding (graph Dijkstra) over 4 cavities completed in 1.80s.
   @> Detected 27 channels.
   @> No output path given.
   @> Channel calculation completed in 6.87s.
   @> Frame/model: 14
   @> WARNING The atoms supplied to calcChannels() contain non-protein components: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein structure, provide an appropriate selection, for example atoms.select('protein').
   @> Delaunay tessellation of 48645 points constructed in 2.15s.
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 A in 0.24s.
   @> Delaunay tessellation of 48645 points constructed in 2.13s.
   @> Surface and inner simplices filtered in 2.73s.
   @> 2 surface cavities detected and filtered in 0.45s.
   @> Delaunay tessellation of 48645 points constructed in 2.29s.
   @> Surface and inner simplices filtered in 2.70s.
   @> Surface and inner simplices filtered in 2.57s.
   @> 1 surface cavities detected and filtered in 0.41s.
   @> 5 surface cavities detected and filtered in 0.39s.
   @> Surface and inner simplices filtered in 2.20s.
   @> 2 surface cavities detected and filtered in 0.35s.
   @> Channel pathfinding (graph Dijkstra) over 5 cavities completed in 1.89s.
   @> Detected 31 channels.
   @> No output path given.
   @> Channel calculation completed in 7.17s.
   @> Frame/model: 13
   @> WARNING The atoms supplied to calcChannels() contain non-protein components: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein structure, provide an appropriate selection, for example atoms.select('protein').
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 A in 0.16s.
   @> Channel pathfinding (graph Dijkstra) over 2 cavities completed in 2.96s.
   @> Detected 28 channels.
   @> No output path given.
   @> Channel calculation completed in 8.31s.
   @> Frame/model: 9
   @> WARNING The atoms supplied to calcChannels() contain non-protein components: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein structure, provide an appropriate selection, for example atoms.select('protein').
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 A in 0.17s.
   @> Channel pathfinding (graph Dijkstra) over 1 cavities completed in 2.65s.
   @> Detected 34 channels.
   @> No output path given.
   @> Channel calculation completed in 8.10s.
   @> Frame/model: 11
   @> WARNING The atoms supplied to calcChannels() contain non-protein components: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein structure, provide an appropriate selection, for example atoms.select('protein').
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 A in 0.17s.
   @> Delaunay tessellation of 48645 points constructed in 1.98s.
   @> Delaunay tessellation of 48645 points constructed in 2.04s.
   @> Channel pathfinding (graph Dijkstra) over 2 cavities completed in 2.90s.
   @> Detected 30 channels.
   @> No output path given.
   @> Channel calculation completed in 8.00s.
   @> Frame/model: 15
   @> WARNING The atoms supplied to calcChannels() contain non-protein components: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein structure, provide an appropriate selection, for example atoms.select('protein').
   @> Delaunay tessellation of 48645 points constructed in 2.06s.
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 A in 0.20s.
   @> Surface and inner simplices filtered in 2.84s.
   @> Surface and inner simplices filtered in 2.70s.
   @> Delaunay tessellation of 48645 points constructed in 2.33s.
   @> 4 surface cavities detected and filtered in 0.44s.
   @> 3 surface cavities detected and filtered in 0.39s.
   @> Surface and inner simplices filtered in 2.90s.
   @> Channel pathfinding (graph Dijkstra) over 3 cavities completed in 0.65s.
   @> Detected 15 channels.
   @> No output path given.
   @> Channel calculation completed in 5.96s.
   @> Frame/model: 16
   @> WARNING The atoms supplied to calcChannels() contain non-protein components: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein structure, provide an appropriate selection, for example atoms.select('protein').
   @> 6 surface cavities detected and filtered in 0.46s.
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 A in 0.20s.
   @> Channel pathfinding (graph Dijkstra) over 4 cavities completed in 2.26s.
   @> Detected 30 channels.
   @> No output path given.
   @> Channel calculation completed in 7.69s.
   @> Surface and inner simplices filtered in 2.38s.
   @> Frame/model: 18
   @> WARNING The atoms supplied to calcChannels() contain non-protein components: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein structure, provide an appropriate selection, for example atoms.select('protein').
   @> Channel pathfinding (graph Dijkstra) over 6 cavities completed in 1.56s.
   @> Detected 28 channels.
   @> No output path given.
   @> Channel calculation completed in 7.16s.
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 A in 0.18s.
   @> 2 surface cavities detected and filtered in 0.40s.
   @> Delaunay tessellation of 48645 points constructed in 1.98s.
   @> Delaunay tessellation of 48645 points constructed in 2.19s.
   @> Channel pathfinding (graph Dijkstra) over 2 cavities completed in 2.46s.
   @> Detected 23 channels.
   @> No output path given.
   @> Channel calculation completed in 7.78s.
   @> Surface and inner simplices filtered in 2.21s.
   @> 5 surface cavities detected and filtered in 0.44s.
   @> Surface and inner simplices filtered in 2.19s.
   @> Channel pathfinding (graph Dijkstra) over 5 cavities completed in 1.24s.
   @> Detected 27 channels.
   @> No output path given.
   @> Channel calculation completed in 6.09s.
   @> Frame/model: 17
   @> WARNING The atoms supplied to calcChannels() contain non-protein components: non-water hetero components: 48 atoms (resnames: RET). All supplied atoms except waters will be used for channel analysis. To analyze only the protein structure, provide an appropriate selection, for example atoms.select('protein').
   @> Substituted 3669 atoms with 48645 homogeneous balls of radius 1.20 A in 0.16s.
   @> 3 surface cavities detected and filtered in 0.38s.
   @> Delaunay tessellation of 48645 points constructed in 1.78s.
   @> Channel pathfinding (graph Dijkstra) over 3 cavities completed in 3.04s.
   @> Detected 36 channels.
   @> No output path given.
   @> Channel calculation completed in 7.99s.
   @> Surface and inner simplices filtered in 2.04s.
   @> 2 surface cavities detected and filtered in 0.32s.
   @> Channel pathfinding (graph Dijkstra) over 2 cavities completed in 2.64s.
   @> Detected 31 channels.
   @> No output path given.
   @> Channel calculation completed in 6.95s.


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

   @> Channel ID:        Volume [Å³]     Length [Å]      Bottleneck [Å]
   @> Frame 0
   @> channel 0:   50.35           7.87            0.94
   @> channel 1:   204.12          23.27           0.91
   @> channel 2:   219.95          30.39           0.75
   @> channel 3:   207.9           31.91           0.75
   @> channel 4:   229.21          33.58           0.75
   @> channel 5:   274.42          34.34           0.82
   @> channel 6:   231.8           35.53           0.82
   @> channel 7:   266.15          39.95           0.87
   @> channel 8:   229.22          40.47           0.84
   @> channel 9:   248.58          42.53           0.84
   @> channel 10:  235.72          41.62           0.84
   @> channel 11:  290.94          44.33           0.84
   @> channel 12:  240.1           42.24           0.84
   @> channel 13:  240.31          41.62           0.84
   @> channel 14:  233.43          42.56           0.84
   @> channel 15:  323.34          50.23           0.87
   @> channel 16:  237.42          43.93           0.84
   @> channel 17:  435.23          61.28           0.87
   @> channel 18:  478.35          64.0            0.87
   @> channel 19:  487.84          64.71           0.87
   @> channel 20:  328.13          54.29           0.87
   @> channel 21:  453.03          63.56           0.87
   @> channel 22:  448.05          64.02           0.87
   @> channel 23:  459.39          65.35           0.87
   @> channel 24:  434.39          64.32           0.87
   @> channel 25:  300.19          52.58           0.86
   @> channel 26:  333.64          55.79           0.87
   @> channel 27:  437.23          65.98           0.87
   @> channel 28:  421.99          64.59           0.87
   @> channel 29:  388.67          65.8            0.87
   @> channel 30:  508.69          79.85           0.87
   @> channel 31:  601.28          82.55           0.87
   @> channel 32:  563.9           87.34           0.87
   @> channel 33:  498.51          82.35           0.87
   @> channel 34:  588.61          89.38           0.87
   @> channel 35:  604.06          95.48           0.87
   @> channel 36:  645.49          98.07           0.87
   @> channel 37:  640.81          101.0           0.87
   @> channel 38:  651.76          101.68          0.87
   @> Frame 1
   @> channel 0:   31.02           6.98            0.88
   @> channel 1:   172.96          18.89           0.89
   @> channel 2:   113.37          14.76           0.95
   @> channel 3:   165.52          20.29           0.89
   @> channel 4:   170.4           22.51           0.89
   @> channel 5:   185.08          24.48           0.89
   @> channel 6:   197.74          25.2            0.89
   @> channel 7:   186.87          25.39           0.89
   @> channel 8:   169.1           24.42           0.76
   @> channel 9:   133.98          23.12           0.76
   @> channel 10:  174.4           25.22           0.76
   @> channel 11:  174.0           26.02           0.91
   @> channel 12:  268.87          31.87           0.91
   @> channel 13:  216.81          32.83           0.89
   @> channel 14:  189.25          28.8            0.91
   @> channel 15:  246.48          32.71           0.91
   @> channel 16:  209.96          33.44           0.89
   @> channel 17:  257.43          33.15           0.88
   @> channel 18:  191.3           31.78           0.76
   @> channel 19:  264.89          41.98           0.91
   @> channel 20:  289.9           41.66           0.76
   @> channel 21:  314.05          45.16           0.76
   @> channel 22:  338.35          46.35           0.76
   @> channel 23:  305.38          43.66           0.76
   @> channel 24:  296.0           46.45           0.76
   @> channel 25:  299.9           50.01           0.76
   @> channel 26:  293.71          46.72           0.26
   @> Frame 2
   @> channel 0:   27.89           5.51            0.91
   @> channel 1:   156.17          23.55           0.93
   @> channel 2:   219.99          29.07           0.93
   @> channel 3:   260.85          33.49           0.93
   @> channel 4:   225.6           31.97           0.93
   @> channel 5:   239.09          33.77           0.93
   @> channel 6:   244.01          34.13           0.92
   @> channel 7:   171.47          29.72           0.85
   @> channel 8:   244.8           36.09           0.89
   @> channel 9:   342.59          38.41           0.92
   @> channel 10:  271.43          37.81           0.92
   @> channel 11:  376.41          44.38           0.93
   @> channel 12:  372.05          45.48           0.93
   @> channel 13:  406.95          45.25           0.92
   @> channel 14:  138.58          26.43           0.68
   @> channel 15:  379.57          44.74           0.92
   @> channel 16:  143.53          27.32           0.68
   @> channel 17:  254.82          39.11           0.93
   @> channel 18:  260.97          35.88           0.83
   @> channel 19:  322.4           38.45           0.83
   @> channel 20:  377.61          46.64           0.85
   @> channel 21:  267.95          42.48           0.88
   @> channel 22:  275.63          42.79           0.93
   @> channel 23:  257.42          42.88           0.83
   @> channel 24:  271.99          45.25           0.88
   @> channel 25:  272.41          45.09           0.88
   @> channel 26:  215.81          36.99           0.5
   @> channel 27:  265.37          40.4            0.79
   @> channel 28:  225.24          41.15           0.79
   @> channel 29:  289.1           50.86           0.83
   @> channel 30:  331.49          55.31           0.83
   @> channel 31:  316.45          56.15           0.83
   @> Frame 3
   @> channel 0:   35.5            6.24            0.96
   @> channel 1:   38.15           7.22            0.95
   @> channel 2:   25.65           7.31            0.8
   @> channel 3:   106.64          15.12           0.73
   @> channel 4:   245.65          26.28           0.94
   @> channel 5:   226.2           21.73           0.84
   @> channel 6:   108.27          17.37           0.73
   @> channel 7:   181.35          20.12           0.73
   @> channel 8:   170.96          24.42           0.94
   @> channel 9:   184.72          26.3            0.94
   @> channel 10:  186.56          24.67           0.84
   @> channel 11:  160.6           24.18           0.94
   @> channel 12:  183.7           26.71           0.94
   @> channel 13:  174.53          26.44           0.84
   @> channel 14:  232.98          28.63           0.84
   @> channel 15:  282.1           31.15           0.94
   @> channel 16:  201.64          27.79           0.84
   @> channel 17:  235.83          31.74           0.94
   @> channel 18:  245.98          32.23           0.94
   @> channel 19:  235.05          31.27           0.93
   @> channel 20:  139.77          23.6            0.77
   @> channel 21:  231.52          32.12           0.93
   @> channel 22:  221.47          29.9            0.84
   @> channel 23:  242.5           34.87           0.84
   @> channel 24:  265.86          36.53           0.84
   @> channel 25:  150.62          27.24           0.77
   @> channel 26:  204.6           30.34           0.84
   @> channel 27:  244.46          36.09           0.84
   @> channel 28:  194.5           35.68           0.94
   @> channel 29:  238.23          35.95           0.78
   @> channel 30:  261.45          38.99           0.84
   @> channel 31:  225.28          36.67           0.77
   @> channel 32:  228.62          37.77           0.77
   @> Frame 4
   @> channel 0:   57.12           8.39            0.8
   @> channel 1:   151.42          14.99           0.94
   @> channel 2:   173.4           15.5            0.94
   @> channel 3:   116.65          16.04           0.94
   @> channel 4:   138.63          17.67           0.94
   @> channel 5:   104.44          16.17           0.94
   @> channel 6:   71.39           14.61           0.77
   @> channel 7:   203.58          22.03           0.77
   @> channel 8:   153.88          22.93           0.77
   @> channel 9:   150.84          24.62           0.77
   @> channel 10:  125.75          22.85           0.77
   @> channel 11:  148.25          26.91           0.77
   @> channel 12:  192.71          29.94           0.74
   @> channel 13:  238.41          33.1            0.74
   @> channel 14:  196.76          31.22           0.74
   @> channel 15:  237.12          35.08           0.74
   @> channel 16:  152.99          29.22           0.77
   @> channel 17:  148.89          28.11           0.72
   @> channel 18:  163.18          30.41           0.77
   @> channel 19:  151.82          29.99           0.77
   @> channel 20:  256.97          37.45           0.74
   @> channel 21:  236.11          36.4            0.74
   @> channel 22:  341.36          42.59           0.74
   @> channel 23:  266.02          40.09           0.74
   @> channel 24:  281.06          40.57           0.74
   @> channel 25:  305.79          43.46           0.74
   @> channel 26:  266.11          42.06           0.74
   @> channel 27:  286.53          44.55           0.74
   @> channel 28:  281.49          44.55           0.74
   @> channel 29:  303.82          47.72           0.74
   @> channel 30:  396.43          51.51           0.74
   @> channel 31:  352.95          51.11           0.74
   @> channel 32:  376.95          53.05           0.74
   @> channel 33:  375.11          54.74           0.74
   @> channel 34:  377.24          55.32           0.74
   @> channel 35:  403.1           57.97           0.74
   @> channel 36:  400.48          59.87           0.74
   @> channel 37:  388.38          57.6            0.74
   @> channel 38:  415.38          61.58           0.74
   @> channel 39:  438.08          60.59           0.74
   @> Frame 5
   @> channel 0:   33.06           6.93            0.82
   @> channel 1:   56.06           8.06            0.82
   @> channel 2:   59.28           9.92            0.84
   @> channel 3:   50.72           9.03            0.88
   @> channel 4:   75.81           11.42           0.84
   @> channel 5:   72.02           13.21           0.84
   @> channel 6:   47.94           9.85            0.68
   @> channel 7:   214.77          22.14           0.92
   @> channel 8:   97.84           16.78           0.84
   @> channel 9:   167.42          22.61           0.92
   @> channel 10:  141.26          23.02           0.92
   @> channel 11:  261.98          30.01           0.85
   @> channel 12:  199.41          26.86           0.85
   @> channel 13:  206.86          29.0            0.9
   @> channel 14:  252.05          32.97           0.85
   @> channel 15:  182.52          27.39           0.85
   @> channel 16:  285.0           35.72           0.85
   @> channel 17:  282.88          36.53           0.85
   @> channel 18:  362.2           40.89           0.85
   @> channel 19:  278.26          37.23           0.85
   @> channel 20:  190.95          30.51           0.85
   @> channel 21:  320.11          40.89           0.85
   @> channel 22:  266.25          36.87           0.85
   @> channel 23:  168.79          30.49           0.85
   @> channel 24:  167.28          30.41           0.85
   @> channel 25:  265.64          36.94           0.85
   @> channel 26:  203.79          33.58           0.85
   @> channel 27:  353.34          42.62           0.85
   @> channel 28:  162.78          30.73           0.85
   @> channel 29:  83.83           22.42           0.76
   @> channel 30:  613.27          56.96           0.85
   @> channel 31:  579.08          56.52           0.85
   @> channel 32:  364.24          46.29           0.85
   @> channel 33:  578.13          57.63           0.85
   @> channel 34:  173.53          33.66           0.85
   @> channel 35:  264.15          37.98           0.67
   @> channel 36:  178.04          34.24           0.85
   @> channel 37:  349.07          46.24           0.85
   @> channel 38:  352.96          46.88           0.85
   @> channel 39:  350.85          46.81           0.85
   @> channel 40:  381.65          48.84           0.85
   @> Frame 6
   @> channel 0:   60.29           8.2             1.06
   @> channel 1:   75.02           9.15            1.06
   @> channel 2:   73.07           8.72            0.9
   @> channel 3:   93.02           15.79           0.95
   @> channel 4:   89.96           16.57           0.88
   @> channel 5:   108.9           18.6            0.9
   @> channel 6:   113.12          19.49           0.91
   @> channel 7:   154.53          23.37           0.88
   @> channel 8:   99.23           21.01           0.91
   @> channel 9:   102.99          22.55           0.91
   @> channel 10:  127.57          24.32           0.79
   @> channel 11:  230.93          32.93           0.9
   @> channel 12:  117.22          25.59           0.8
   @> channel 13:  298.65          39.73           0.9
   @> channel 14:  360.0           41.8            0.9
   @> channel 15:  351.27          45.31           0.9
   @> channel 16:  310.19          44.64           0.9
   @> channel 17:  345.8           52.04           0.9
   @> channel 18:  339.05          53.52           0.9
   @> channel 19:  339.07          54.6            0.9
   @> channel 20:  346.32          55.14           0.9
   @> Frame 7
   @> channel 0:   51.7            7.81            0.93
   @> channel 1:   43.53           8.56            0.93
   @> channel 2:   55.43           10.19           0.93
   @> channel 3:   54.88           9.71            0.9
   @> channel 4:   179.89          16.38           0.92
   @> channel 5:   46.34           9.79            0.9
   @> channel 6:   317.97          28.22           0.91
   @> channel 7:   148.37          17.74           0.91
   @> channel 8:   316.77          29.28           0.91
   @> channel 9:   149.28          20.28           0.92
   @> channel 10:  168.69          21.67           0.92
   @> channel 11:  258.05          28.94           0.91
   @> channel 12:  113.13          18.56           0.92
   @> channel 13:  180.19          24.46           0.91
   @> channel 14:  324.49          32.0            0.91
   @> channel 15:  204.96          24.78           0.84
   @> channel 16:  188.8           24.85           0.84
   @> channel 17:  328.4           33.02           0.91
   @> channel 18:  389.09          38.61           0.91
   @> channel 19:  261.24          31.56           0.91
   @> channel 20:  229.59          29.46           0.91
   @> channel 21:  166.12          22.9            0.66
   @> channel 22:  452.23          42.72           0.91
   @> channel 23:  272.16          33.26           0.91
   @> channel 24:  497.59          44.38           0.91
   @> channel 25:  383.66          39.7            0.91
   @> channel 26:  141.4           24.72           0.66
   @> channel 27:  139.39          24.7            0.66
   @> channel 28:  137.51          25.17           0.66
   @> channel 29:  147.52          26.3            0.66
   @> channel 30:  142.04          23.21           0.54
   @> channel 31:  201.62          29.34           0.66
   @> channel 32:  147.11          27.22           0.66
   @> channel 33:  395.77          45.85           0.91
   @> channel 34:  376.65          45.37           0.91
   @> channel 35:  358.16          45.06           0.91
   @> channel 36:  338.59          43.81           0.85
   @> channel 37:  355.91          45.61           0.91
   @> channel 38:  145.92          27.85           0.54
   @> channel 39:  178.12          31.04           0.54
   @> Frame 8
   @> channel 0:   48.41           6.88            0.92
   @> channel 1:   51.17           8.75            0.92
   @> channel 2:   62.54           9.41            0.92
   @> channel 3:   169.68          26.06           0.92
   @> channel 4:   187.25          28.68           0.92
   @> channel 5:   300.98          30.74           0.84
   @> channel 6:   168.12          28.04           0.92
   @> channel 7:   181.47          29.94           0.92
   @> channel 8:   284.34          34.87           0.84
   @> channel 9:   293.7           35.51           0.84
   @> channel 10:  241.84          33.48           0.84
   @> channel 11:  271.64          35.24           0.84
   @> channel 12:  294.59          38.43           0.84
   @> channel 13:  227.73          37.08           0.92
   @> channel 14:  169.32          36.07           0.82
   @> channel 15:  304.84          49.38           0.92
   @> channel 16:  282.68          41.42           0.8
   @> channel 17:  283.72          47.51           0.92
   @> channel 18:  337.24          48.55           0.92
   @> channel 19:  177.95          38.12           0.82
   @> channel 20:  169.03          37.21           0.82
   @> channel 21:  327.41          51.7            0.92
   @> channel 22:  305.44          50.34           0.9
   @> channel 23:  320.9           52.71           0.92
   @> channel 24:  378.18          56.62           0.88
   @> channel 25:  378.59          56.69           0.88
   @> channel 26:  326.3           57.68           0.36
   @> channel 27:  374.74          68.81           0.36
   @> Frame 9
   @> channel 0:   35.5            6.87            0.92
   @> channel 1:   50.67           10.43           0.91
   @> channel 2:   75.65           14.35           0.88
   @> channel 3:   177.85          22.96           0.94
   @> channel 4:   253.86          32.51           0.94
   @> channel 5:   270.65          36.7            0.94
   @> channel 6:   214.63          33.03           0.76
   @> channel 7:   296.99          41.78           0.93
   @> channel 8:   297.0           41.66           0.93
   @> channel 9:   262.81          40.53           0.9
   @> channel 10:  269.34          41.09           0.93
   @> channel 11:  295.95          44.4            0.93
   @> channel 12:  415.1           57.63           0.77
   @> channel 13:  376.9           59.64           0.77
   @> channel 14:  406.08          62.93           0.77
   ..
   ..
   @> Frame 18
   @> channel 0:   24.52           5.65            0.91
   @> channel 1:   22.68           5.4             0.86
   @> channel 2:   38.07           9.51            0.91
   @> channel 3:   137.08          32.07           0.9
   @> channel 4:   139.86          33.05           0.86
   @> channel 5:   232.34          41.27           0.9
   @> channel 6:   343.38          47.78           0.9
   @> channel 7:   272.81          43.69           0.9
   @> channel 8:   246.6           42.78           0.9
   @> channel 9:   369.35          50.59           0.9
   @> channel 10:  246.61          38.6            0.83
   @> channel 11:  298.28          43.4            0.83
   @> channel 12:  334.6           45.16           0.83
   @> channel 13:  278.74          46.01           0.9
   @> channel 14:  366.33          47.07           0.83
   @> channel 15:  321.08          48.51           0.9
   @> channel 16:  359.88          47.66           0.83
   @> channel 17:  212.05          43.86           0.84
   @> channel 18:  239.35          45.28           0.84
   @> channel 19:  231.45          44.74           0.83
   @> channel 20:  412.58          56.24           0.83
   @> channel 21:  272.12          52.18           0.9
   @> channel 22:  410.73          56.77           0.83
   @> channel 23:  251.33          52.34           0.9
   @> channel 24:  261.48          51.3            0.87
   @> channel 25:  234.63          49.75           0.89
   @> channel 26:  246.15          52.26           0.9
   @> channel 27:  283.84          55.29           0.9
   @> channel 28:  406.37          57.48           0.83
   @> channel 29:  258.34          52.26           0.87
   @> channel 30:  243.38          52.79           0.9
   @> channel 31:  271.71          53.06           0.85
   @> channel 32:  223.31          50.92           0.89
   @> channel 33:  338.14          61.82           0.9
   @> channel 34:  354.53          62.5            0.9
   @> channel 35:  339.9           62.75           0.78

   [([7.865188922554236,
      23.268497762876585,
      30.39292266480067,
      31.906393020223632,
      33.57982277826004,
      34.3379231113981,
      35.53440706573264,
      39.94923607322629,
      40.472931196510665,
      42.53361013796159,
      41.62124888584993,
      44.32781411583137,
      42.23838649553053,
      41.61799131584387,
      42.55773739475321,
      50.22662719969935,
      43.93440615542735,
      61.2784657630814,
      64.00267058795325,
      64.70579348856032,
      54.289518800090875,
      63.557047844768014,
      64.01844082976504,
      65.34821900024978,
      64.31849421681108,
      52.583266494536346,
      55.785182496848414,
      65.9834289205483,
      64.58524819287172,
      65.80361633791189,
      79.85284824091534,
      82.55216710150766,
      87.34461993557275,
      82.34781023408898,
      89.3801156860001,
      95.48174882254847,
      98.07393309466514,
      100.99971784754959,
      101.6757430476978],
     [0.9375858145821543,
      0.9095124228459006,
      ..
      ..
      251.32842014632106,
      261.48169093480885,
      234.63406796892775,
      246.15089833842134,
      283.8371393582868,
      406.3725187354091,
      258.33865059227685,
      243.37849296530914,
      271.7094473114096,
      223.31221990418078,
      338.137435087287,
      354.53100792973135,
      339.902566902574])]

Access to the residues that are forming the channels is provided by
:func:`.getChannelResidueNames` function. Below, the example on how to
obtain information for frame #0.

.. ipython:: python
   :verbatim:

   p3.setACSIndex(0)
   getChannelResidueNames(p3, channels3[0])

.. parsed-literal::

   @> 3769 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3879 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3919 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3924 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3889 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3964 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4014 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3994 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3979 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4029 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4029 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4014 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4019 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4024 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4024 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4154 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4049 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4289 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4264 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4279 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4159 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4254 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4284 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4289 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4254 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4174 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4164 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4279 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4294 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4299 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4424 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4434 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4444 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4444 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4459 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4544 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4564 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4604 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4594 atoms and 1 coordinate set(s) were parsed in 0.05s.
   Out[18]: 
   ['channel0: ILE106, CYS107, TYR110, LEU111, ALA120, SER122, LEU123, LYS126',
    'channel1: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, TRP83, ARG94, PHE137, GLY141, TYR223, ASN224, ASP227, PHE228',
    'channel2: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, GLY196, TRP197, ILE199, TYR200, TYR223, ASN224, LEU225, ALA226, ASP227, PHE228, VAL229, ASN230',
    'channel3: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, ILE192, PHE195, GLY196, TRP197, ILE199, TYR200, TYR223, ASN224, ALA226, ASP227, PHE228, VAL229, ASN230',
    'channel4: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, GLY196, TRP197, ILE199, TYR200, GLY203, ILE222, TYR223, ASN224, LEU225, ALA226, ASP227, PHE228, ASN230',
    'channel5: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, PHE152, GLY155, CYS156, TRP159, TRP197, ALA198, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230',
    'channel6: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, GLY155, CYS156, TRP159, ILE194, TRP197, ALA198, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230',
    'channel7: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, ARG94, TYR95, TRP98, LEU99, VAL102, PRO103, LEU127, GLY130, SER131, MET134, PHE137, TYR223, ASN224, ALA226, ASP227, PHE228',
    'channel8: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, ARG94, TYR95, TRP98, LEU99, VAL102, PRO103, GLY130, SER131, MET134, PHE137, TYR223, ASN224, ALA226, ASP227, PHE228',
    'channel9: GLY21, VAL36, LEU40, PHE73, TYR76, MET77, ARG80, TYR200, TYR204, GLY207, TYR208, MET210, GLY211, ASP212, GLY214, SER215, ASN218, LEU219, ILE222, TYR223, ASN224, ALA226, ASP227, PHE228',
    'channel10: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, TYR200, TYR204, GLY207, TYR208, LEU209, MET210, GLY211, ASP212, GLY213, GLY214, SER215, ASN218, LEU219, ILE222, TYR223, ASN224, ALA226, ASP227, PHE228',
    'channel11: ASP22, LEU23, VAL36, LEU40, PHE73, TYR76, MET77, ARG80, TYR200, TYR204, GLY207, TYR208, MET210, GLY211, ASP212, GLY214, SER215, ASN218, LEU219, ILE222, TYR223, ASN224, ALA226, ASP227, PHE228',
    'channel12: GLY21, VAL36, LEU40, PHE73, TYR76, MET77, ARG80, TYR200, TYR204, THR206, GLY207, TYR208, MET210, GLY211, GLY214, SER215, ASN218, LEU219, ILE222, TYR223, ASN224, ALA226, ASP227, PHE228',
    'channel13: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, TYR200, TYR204, GLY207, TYR208, GLY211, GLY214, SER215, ALA216, ASN218, LEU219, ILE222, TYR223, ASN224, ALA226, ASP227, PHE228',
    'channel14: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, TYR200, TYR204, GLY207, TYR208, LEU209, MET210, GLY211, ASP212, GLY214, SER215, ASN218, LEU219, ILE222, TYR223, ASN224, ALA226, ASP227, PHE228',
    'channel15: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, LYS125, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, TYR161, MET162, ILE163, GLU165, LEU166, MET190, ILE193, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230',
    'channel16: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, TYR200, TYR204, GLY207, TYR208, MET210, GLY211, ASP212, GLY213, GLY214, SER215, ASN218, LEU219, ILE222, TYR223, ASN224, ALA226, ASP227, PHE228',
    'channel17: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, PHE109, TYR110, LEU113, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, LEU166, TRP167, GLY169, GLY171, LYS172, TYR186, MET189, MET190, ILE193, ILE194, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230',
    'channel18: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, PHE109, TYR110, LEU113, ALA114, ALA116, THR117, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, LEU166, GLY171, ALA174, ALA178, SER179, VAL182, TYR186, MET189, MET190, ILE193, ILE194, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230',
    'channel19: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, PHE109, TYR110, LEU113, ALA114, ALA116, THR117, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, LEU166, GLY171, ALA174, CYS175, ALA178, VAL182, TYR186, MET189, MET190, ILE193, ILE194, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230',
    'channel20: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, SER122, LYS125, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, GLU165, LEU166, GLY169, GLU170, MET190, ILE193, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230',
    'channel21: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, PHE109, TYR110, LEU113, ALA114, THR117, VAL119, ALA120, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, LEU166, GLU170, GLY171, ALA174, VAL182, TYR186, MET189, MET190, ILE193, ILE194, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230',
    'channel22: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, PHE109, TYR110, LEU113, ALA114, ALA120, GLY121, SER122, LYS125, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, LEU166, GLU170, GLY171, TYR186, MET189, MET190, ILE193, ILE194, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230',
    'channel23: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, PHE109, TYR110, LEU113, ALA114, VAL119, ALA120, GLY121, SER122, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, LEU166, GLU170, GLY171, TYR186, MET189, MET190, ILE193, ILE194, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230',
    'channel24: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, PHE109, TYR110, LEU113, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, LEU166, GLY171, LYS172, VAL182, GLN183, TYR186, MET189, MET190, ILE193, ILE194, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230',
    'channel25: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, LYS125, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, GLU165, LEU166, GLY169, MET190, ILE193, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230',
    'channel26: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, SER122, LYS125, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, GLU165, LEU166, ALA168, GLY169, MET190, ILE193, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230',
    'channel27: VAL36, LEU40, PHE73, TYR76, MET77, ARG80, PHE109, TYR110, LEU113, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, LEU166, GLY171, LYS172, THR177, ALA178, SER179, VAL182, GLN183, TYR186, MET189, MET190, ILE193, ILE194, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230',
    'channel28: VAL36, LEU40, LYS57, TRP58, PHE73, TYR76, MET77, ARG80, PHE109, TYR110, LEU113, ALA116, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, LEU166, ALA181, VAL182, ALA185, TYR186, MET189, MET190, ILE193, ILE194, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230',
    'channel29: VAL36, LEU40, TRP58, PHE73, TYR76, MET77, ARG80, PHE109, TYR110, LEU113, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, LEU166, SER179, ALA181, VAL182, ALA185, TYR186, MET189, MET190, ILE193, ILE194, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230',
    'channel30: VAL36, LEU40, PHE47, TRP58, SER61, LEU62, SER65, GLY66, THR69, GLY70, PHE73, TYR76, MET77, ARG80, PHE109, TYR110, LEU113, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, LEU166, VAL182, ALA185, TYR186, MET189, MET190, ILE193, ILE194, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230, PHE234, GLY235, ILE238',
    'channel31: VAL36, LEU40, THR44, PHE47, TRP58, SER61, LEU62, SER65, GLY66, THR69, GLY70, PHE73, TYR76, MET77, ARG80, PHE109, TYR110, LEU113, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, LEU166, VAL182, ALA185, TYR186, MET189, MET190, ILE193, ILE194, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230, ILE232, PHE234, GLY235, ILE238',
    'channel32: VAL36, LEU40, SER43, PHE46, PHE47, TRP58, SER61, LEU62, SER65, GLY66, THR69, PHE73, TYR76, MET77, ARG80, PHE109, TYR110, LEU113, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, LEU166, VAL182, ALA185, TYR186, MET189, MET190, ILE193, ILE194, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230, ILE232, PHE234, GLY235, LEU236, ILE238',
    'channel33: VAL36, LEU40, PHE47, TRP58, SER61, LEU62, SER65, GLY66, LEU67, THR69, GLY70, PHE73, TYR76, MET77, ARG80, PHE109, TYR110, LEU113, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, LEU166, VAL182, ALA185, TYR186, MET189, MET190, ILE193, ILE194, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230, PHE234, GLY235, ILE238',
    'channel34: VAL36, LEU40, SER43, PHE46, PHE47, GLU50, TRP58, SER61, LEU62, SER65, GLY66, THR69, PHE73, TYR76, MET77, ARG80, PHE109, TYR110, LEU113, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, LEU166, VAL182, ALA185, TYR186, MET189, MET190, ILE193, ILE194, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230, ILE232, PHE234, GLY235, LEU236, ILE238, TRP239',
    'channel35: VAL36, LEU40, SER43, PHE46, PHE47, PHE48, GLU50, ARG51, TRP58, SER61, LEU62, SER65, GLY66, THR69, PHE73, TYR76, MET77, ARG80, PHE109, TYR110, LEU113, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, LEU166, VAL182, ALA185, TYR186, MET189, MET190, ILE193, ILE194, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230, ILE232, PHE234, GLY235, LEU236, ILE238, TRP239',
    'channel36: VAL36, LEU40, SER43, PHE46, PHE47, GLU50, ARG51, TRP58, SER61, LEU62, THR63, SER65, GLY66, THR69, PHE73, TYR76, MET77, ARG80, PHE109, TYR110, LEU113, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, LEU166, VAL182, ALA185, TYR186, MET189, MET190, ILE193, ILE194, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230, ILE232, PHE234, GLY235, LEU236, ILE238, TRP239',
    'channel37: VAL36, LEU40, SER43, PHE46, PHE47, GLU50, ARG51, VAL54, TRP58, LYS59, SER61, LEU62, SER65, GLY66, THR69, PHE73, TYR76, MET77, ARG80, PHE109, TYR110, LEU113, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, LEU166, VAL182, ALA185, TYR186, MET189, MET190, ILE193, ILE194, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230, ILE232, PHE234, GLY235, LEU236, ILE238, TRP239',
    'channel38: VAL36, LEU40, SER43, PHE46, PHE47, GLU50, ARG51, ASP52, VAL54, TRP58, LYS59, SER61, LEU62, SER65, GLY66, THR69, PHE73, TYR76, MET77, ARG80, PHE109, TYR110, LEU113, LYS126, VAL129, VAL133, GLY155, ALA158, TRP159, MET162, ILE163, LEU166, VAL182, ALA185, TYR186, MET189, MET190, ILE193, ILE194, TRP197, TYR200, PRO201, TYR223, ASN224, ALA226, ASP227, PHE228, ASN230, ILE232, PHE234, GLY235, LEU236, ILE238, TRP239']

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

   @> 3739 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3884 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3794 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3899 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3914 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3929 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3939 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3944 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3869 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3854 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3859 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3879 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3939 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3974 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3914 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3944 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4004 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3934 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3934 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4009 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4034 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4009 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4039 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4034 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4094 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4114 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4059 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> Channel residues were saved to: results_frame1_Residues_All_channels.txt

   ['channel0: W159, M162, I163, L166, W167, M190, Y191, I194',
    'channel1: S43, F47, V54, W58, L62, S65, G66, T69, G70, I232, F234, G235, I238, W239, A242, V243',
    'channel2: F33, V36, L40, F73, Y76, M77, Y78, R80, G81, N224, D227, F228',
    'channel3: S43, F47, V54, W58, L62, S65, G66, T69, I232, F234, G235, L236, I238, W239, A242, V243',
    'channel4: F47, V54, K57, W58, S61, L62, S65, F109, I112, L113, G235, I238, W239, A242, V243',
    'channel5: F47, V54, K57, W58, S61, L62, S65, F109, I112, L113, A116, A181, V182, G235, I238, W239, A242, V243',
    'channel6: F47, V54, K57, W58, S61, L62, S65, F109, I112, L113, A116, V182, G235, I238, W239, A242, V243',
    'channel7: S43, F46, F47, E50, V54, W58, L62, S65, G66, T69, I232, F234, G235, L236, I238, W239, N240, A242, V243',
    'channel8: V64, S65, V68, T69, W98, L105, E108, F109, M189, I192, I193, F195, G196, I199, V229, N230, L233, F234',
    'channel9: V64, S65, V68, T69, W98, L105, E108, F109, M189, I192, I193, V229, N230, L233, F234, I237',
    'channel10: V64, S65, V68, T69, W98, L105, E108, F109, M189, I192, I193, G196, V229, N230, L233, F234, I237',
    'channel11: V64, S65, V68, L105, I106, E108, F109, Y110, L113, A120, G121, S122, K126, G169, E170, G171, Y186, M189, F234',
    'channel12: V64, S65, V68, L105, I106, E108, F109, Y110, L113, A114, T117, A120, S122, K126, E170, A174, T177, A178, V182, Y186, M189, F234',
    'channel13: F47, V54, S55, K57, W58, S61, L62, S65, F109, I112, L113, A181, V182, A185, G235, I238, W239, V241, A242, V243, E245, S246',
    'channel14: V64, S65, V68, L105, I106, E108, F109, Y110, L113, A114, V119, A120, G121, S122, K126, E170, G171, Y186, M189, F234',
    'channel15: V64, S65, V68, L105, I106, E108, F109, Y110, L113, A114, A116, T117, A120, S122, K126, E170, A174, A178, S179, V182, Y186, M189, F234',
    'channel16: F47, V54, W58, S61, L62, S65, F109, I112, L113, A181, V182, S184, A185, G235, I238, W239, V241, A242, V243, K244, E245',
    'channel17: V64, S65, V68, L105, I106, E108, F109, Y110, L113, A114, T117, A120, S122, K126, E170, A174, T177, A178, V182, Y186, M189, F234',
    'channel18: V64, S65, V68, T69, W98, L99, V102, P103, L105, I106, E108, F109, K126, L127, G130, S131, V133, M134, M189, I193, W197, N230, F234',
    'channel19: V64, S65, V68, L105, I106, E108, F109, Y110, L113, A120, S122, K125, K126, L128, V129, Y161, M162, E165, L166, G169, E170, G171, Y186, M189, F234',
    'channel20: V64, S65, V68, T69, W98, V102, L105, I106, E108, F109, K126, G130, V133, M134, F152, G155, C156, A158, W159, M189, I193, W197, A198, P201, N230, F234',
    'channel21: V64, S65, V68, T69, Y76, T91, R94, Y95, W98, V102, L105, I106, E108, F109, K126, G130, V133, M134, F137, G138, E142, M189, I193, W197, N230, F234',
    'channel22: V64, S65, V68, T69, Y76, T91, R94, Y95, W98, V102, L105, I106, E108, F109, K126, G130, V133, M134, F137, G138, G141, M189, I193, W197, N230, F234',
    'channel23: V64, S65, V68, T69, Y76, R94, W98, V102, L105, I106, E108, F109, K126, G130, V133, M134, F137, M189, I193, W197, Y200, Y223, N230, F234',
    'channel24: V64, S65, V68, T69, Y76, R80, W83, R94, W98, V102, L105, I106, E108, F109, K126, G130, V133, M134, F137, M189, I193, W197, Y200, L219, Y223, N230, F234',
    'channel25: V64, S65, V68, T69, Y76, R80, W83, I84, R94, W98, V102, L105, I106, E108, F109, K126, G130, V133, M134, F137, M189, I193, W197, Y200, A216, L219, N220, Y223, N230, F234',
    'channel26: V64, S65, V68, L105, I106, E108, F109, Y110, L113, A120, S122, K125, K126, V129, L157, A158, Y161, M162, E165, L166, G169, E170, G171, Y186, M189, F234']


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
:func:`.calcSurfaceCavitiesMultipleFrames`. We used ``r2=1.5`` controlling
the detection of accessible surface cavities. 
The results are also saved as separate PQR files for individual cavities
when ``separate`` parameter is set.


.. ipython:: python
   :verbatim:

   cavities, surface = calcSurfaceCavitiesMultipleFrames(atoms, r2=1.5, 
			output_path=PDB_ID+'_CAV_', separate=True)   

.. parsed-literal::

   @> Model: 0
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 2360 atoms with 30586 homogeneous balls of radius 1.20 A in 0.18s.
   @> Delaunay tessellation of 30586 points constructed in 1.26s.
   @> Surface and inner simplices filtered in 0.38s.
   @> 26 surface cavities detected and filtered in 0.27s.
   @> Returning surface cavities
   @> Saving multiple surface cavities to directory ..
   @> Surface cavity calculation completed in 2.21s.
   @> Model: 1
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 2360 atoms with 30586 homogeneous balls of radius 1.20 A in 0.10s.
   @> Delaunay tessellation of 30586 points constructed in 1.21s.
   @> Surface and inner simplices filtered in 0.35s.
   @> 32 surface cavities detected and filtered in 0.26s.
   @> Returning surface cavities
   @> Saving multiple surface cavities to directory ..
   @> Surface cavity calculation completed in 2.03s.
   @> Model: 2
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 2360 atoms with 30586 homogeneous balls of radius 1.20 A in 0.10s.
   @> Delaunay tessellation of 30586 points constructed in 1.16s.
   @> Surface and inner simplices filtered in 0.33s.
   @> 29 surface cavities detected and filtered in 0.26s.
   @> Returning surface cavities
   @> Saving multiple surface cavities to directory ..
   @> Surface cavity calculation completed in 1.97s.
   ..
   ..
   @> Model: 18
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 2360 atoms with 30586 homogeneous balls of radius 1.20 A in 0.11s.
   @> Delaunay tessellation of 30586 points constructed in 1.15s.
   @> Surface and inner simplices filtered in 0.38s.
   @> 25 surface cavities detected and filtered in 0.27s.
   @> Returning surface cavities
   @> Saving multiple surface cavities to directory ..
   @> Surface cavity calculation completed in 2.03s.
   @> Model: 19
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 2360 atoms with 30586 homogeneous balls of radius 1.20 A in 0.10s.
   @> Delaunay tessellation of 30586 points constructed in 1.16s.
   @> Surface and inner simplices filtered in 0.38s.
   @> 25 surface cavities detected and filtered in 0.26s.
   @> Returning surface cavities
   @> Saving multiple surface cavities to directory ..
   @> Surface cavity calculation completed in 2.02s.

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
   @> cavity 0: 	521.74 		4.69 		308
   @> cavity 1: 	599.86 		6.67 		360
   @> cavity 2: 	80.0 		2.24 		95
   @> cavity 3: 	132.38 		3.98 		127
   @> cavity 4: 	54.74 		1.9 		43
   @> cavity 5: 	484.91 		6.86 		156
   @> cavity 6: 	369.42 		8.6 		181
   @> cavity 7: 	343.3 		6.07 		160
   @> cavity 8: 	67.58 		3.14 		66
   @> cavity 9: 	389.91 		4.73 		297
   @> cavity 10: 	77.69 		1.86 		39
   @> cavity 11: 	81.99 		3.14 		52
   @> cavity 12: 	65.45 		3.93 		58
   @> cavity 13: 	96.49 		3.23 		70
   @> Model/frame: 1
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	1584.91 		7.96 		980
   @> cavity 1: 	78.4 		2.06 		75
   @> cavity 2: 	164.29 		2.28 		157
   @> cavity 3: 	174.84 		6.11 		71
   @> cavity 4: 	60.96 		1.59 		60
   @> cavity 5: 	211.44 		6.01 		110
   @> cavity 6: 	118.85 		3.8 		63
   @> cavity 7: 	66.94 		2.9 		48
   @> cavity 8: 	64.48 		4.45 		78
   @> cavity 9: 	123.68 		3.5 		59
   @> cavity 10: 	303.16 		4.19 		215
   @> cavity 11: 	180.99 		4.45 		105
   @> cavity 12: 	85.18 		5.62 		45
   @> Model/frame: 2
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	408.56 		8.82 		150
   @> cavity 1: 	409.74 		6.44 		222
   @> cavity 2: 	135.95 		2.56 		110
   @> cavity 3: 	704.89 		5.96 		451
   @> cavity 4: 	1041.57 	7.92 		504
   @> cavity 5: 	282.52 		6.58 		161
   @> cavity 6: 	242.03 		6.5 		151
   @> cavity 7: 	85.45 		2.02 		79
   @> cavity 8: 	54.49 		1.64 		40
   @> cavity 9: 	74.31 		3.18 		74
   @> Model/frame: 3
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	53.98 		2.47 		37
   @> cavity 1: 	186.54 		2.92 		120
   @> cavity 2: 	645.63 		4.06 		422
   @> cavity 3: 	104.42 		3.87 		91
   @> cavity 4: 	533.04 		7.76 		267
   @> cavity 5: 	316.88 		7.25 		154
   @> cavity 6: 	120.94 		5.19 		80
   @> cavity 7: 	498.7 		5.82 		313
   @> cavity 8: 	743.6 		11.04 		488
   @> cavity 9: 	113.18 		6.39 		61
   @> cavity 10: 	67.45 		1.69 		78
   @> cavity 11: 	65.77 		2.11 		41
   ..
   ..
   @> Model/frame: 18
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	103.13 		3.3 		89
   @> cavity 1: 	942.27 		7.62 		585
   @> cavity 2: 	536.1 		6.4 		378
   @> cavity 3: 	153.08 		2.18 		118
   @> cavity 4: 	1329.65 	7.74 		620
   @> cavity 5: 	400.23 		5.67 		248
   @> cavity 6: 	85.65 		1.84 		68
   @> cavity 7: 	206.96 		12.19 		114
   @> cavity 8: 	153.19 		2.51 		107
   @> cavity 9: 	74.67 		2.51 		68
   @> cavity 10: 	104.02 		3.82 		43
   @> Model/frame: 19
   @> Cavity ID: 	Volume [Å³] 	Depth [Å] 	Tetrahedra count
   @> cavity 0: 	546.92 		6.01 		302
   @> cavity 1: 	1243.19 	4.59 		614
   @> cavity 2: 	273.83 		6.53 		105
   @> cavity 3: 	1462.62 	7.46 		873
   @> cavity 4: 	58.85 		3.6 		52
   @> cavity 5: 	151.99 		2.73 		135
   @> cavity 6: 	154.38 		5.37 		79
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
   @> 1979 atoms and 1 coordinate sets were parsed in 0.03s.
   @> 1942 atoms and 1 coordinate sets were parsed in 0.03s.
   @> 1988 atoms and 1 coordinate sets were parsed in 0.05s.
   @> 2148 atoms and 1 coordinate sets were parsed in 0.05s.
   @> 2147 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1751 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2034 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2152 atoms and 1 coordinate sets were parsed in 0.02s.
   @> 2066 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2012 atoms and 1 coordinate sets were parsed in 0.02s.
   @> 2033 atoms and 1 coordinate sets were parsed in 0.02s.
   @> 2207 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1958 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1934 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2091 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2093 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2438 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2378 atoms and 1 coordinate sets were parsed in 0.03s.
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

   channels, surface, details = calcChannelsMultipleFrames(pdb_multi, r2=0.85,
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
   @> Model: 6
   @> 6636 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6651 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6621 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6636 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> Model: 7
   @> 6481 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6561 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6466 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6481 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6561 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6556 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6486 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6601 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6566 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6606 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6571 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> Model: 8
   @> 6531 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> Model: 9
   @> Model: 10
   @> 6461 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6471 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6516 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6526 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> Model: 11
   @> Model: 12
   @> 6476 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6511 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6561 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6526 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6556 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6521 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6531 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6491 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6551 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6531 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> Model: 13
   @> 6616 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6676 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6676 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6676 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6676 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6656 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6656 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> Model: 14
   @> Model: 15
   @> 6501 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6556 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6546 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6576 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6541 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6581 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6526 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6516 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6571 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6561 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6591 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6556 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6596 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6541 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6611 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6661 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6501 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6556 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6546 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6576 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6541 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6581 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6526 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6571 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6561 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6591 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6556 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6596 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6541 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> Model: 16
   @> Model: 17
   @> 6676 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6661 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6676 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6711 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6736 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6701 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6651 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6691 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6676 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6691 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6726 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6751 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6716 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6666 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6716 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6701 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6716 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6751 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6776 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6741 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6691 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6676 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6661 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6676 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6711 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6736 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6701 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6651 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6671 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6656 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6671 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6706 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6731 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6696 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6646 atoms and 1 coordinate set(s) were parsed in 0.06s.
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
    ['pore0: V26, L30, V41, P42, I43, I44, P45, S46, S200, M204, L207, E215, G218, M221, G222, A224, L225, L228, Y243, K248, E312, I317, W318, M319, T322, L330, F348, G349, H353, R357, D399, S400, S401, M403, P404, D426, F429',
     'pore1: V26, L30, V41, P42, I43, I44, P45, S46, S200, M204, L207, G218, M221, G222, A224, L225, L228, Y243, K248, E312, I317, W318, M319, T322, L330, F348, G349, A352, H353, M355, G356, R357, D399, S400, S401, M403, P404, D426, F429',
     'pore2: V26, L30, V41, P42, I43, I44, P45, S46, S200, M204, L207, E215, G218, M221, G222, A224, L225, L228, Y243, E312, I317, W318, M319, T322, M323, L330, F348, G349, H353, R357, D399, S400, S401, M403, P404, D426, F429',
     'pore3: V26, L30, V41, P42, I43, I44, P45, S46, S200, M204, L207, G218, M221, G222, A224, L225, L228, Y243, E312, I317, W318, M319, T322, M323, L330, F348, G349, A352, H353, M355, G356, R357, D399, S400, S401, M403, P404, D426, F429'],
    ['pore0: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, S46, S119, E120, L124, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L228, L270, E312, W318, M319, F334, D426, F429',
     'pore1: L30, D33, N34, L37, V41, I44, S46, S119, E120, L124, N128, V131, F135, K138, Q142, S200, M204, L207, A208, Y211, T212, R217, V220, M221, A224, L225, L228, E312, W318, M319, F334, M403, P404, G407, D411, Y418, Y422, D426, F429',
     'pore2: I22, I25, V26, A29, L30, D33, N34, L37, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L228, L270, L311, E312, L315, P316, I317, Q329, V332, A333, A337, I381, L384, I385, N388, D426, F429',
     'pore3: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, S46, L124, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L228, L270, E312, W318, M319, F334, D426, F429',
     'pore4: L30, D33, N34, L37, V41, I44, S46, L124, N128, V131, F135, K138, Q142, S200, M204, L207, A208, Y211, T212, R217, V220, M221, A224, L225, L228, E312, W318, M319, F334, M403, P404, G407, D411, Y418, Y422, D426, F429',
     'pore5: L30, D33, N34, L37, V41, P42, I43, I44, P45, S46, F135, K138, Q142, S200, M204, L207, A208, Y211, T212, R217, V220, M221, A224, L225, L228, Y243, E312, I317, W318, M319, T322, M323, F334, M403, P404, G407, D411, Y418, Y422, D426, F429',
     'pore6: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, S46, S119, E120, L124, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L228, L270, E312, W318, F334, D426, F429',
     'pore7: L30, D33, N34, L37, V41, I44, S46, S119, E120, L124, N128, V131, F135, K138, Q142, S200, M204, L207, A208, Y211, D214, E215, R217, G218, V220, M221, A224, L225, L228, E312, W318, F334, H353, R357, M403, P404, Y422, D426, F429',
     'pore8: L30, D33, N34, L37, V41, I44, S46, S119, E120, L124, N128, V131, F135, K138, Q142, S200, M204, L207, A208, Y211, T212, R217, V220, M221, A224, L225, L228, E312, W318, F334, M403, P404, G407, D411, Y418, Y422, D426, F429',
     'pore9: L30, D33, N34, L37, V41, P42, I44, P45, S46, F135, K138, Q142, S200, M204, L207, A208, Y211, D214, E215, R217, G218, V220, M221, A224, L225, L228, Y243, E312, I317, W318, M319, T322, M323, S325, F334, H353, R357, M403, P404, Y422, D426, F429',
     'pore10: L30, D33, N34, L37, V41, P42, I44, P45, S46, F135, K138, Q142, S200, M204, L207, A208, Y211, T212, R217, V220, M221, A224, L225, L228, Y243, E312, I317, W318, M319, T322, M323, S325, F334, M403, P404, G407, D411, Y418, Y422, D426, F429'],
    ['pore0: I22, I25, V26, A29, L30, D33, L37, V41, L124, F135, K138, Q142, P159, I162, S196, S199, S200, G203, M206, L207, L228, V269, L270, I308, L311, E312, P313, L315, P316, I317, W318, F334, G392'],
    [],
    ['pore0: I22, I25, V26, A29, L30, D33, N34, L37, L124, K138, I162, S196, S199, S200, G203, M204, M206, L225, L228, V232, V269, L270, I308, E312, P313, L315, P316, I317, W318, K379, M403, D426, F429',
     'pore1: I22, I25, V26, A29, L30, D33, N34, L37, L124, K138, I162, S196, S199, S200, G203, M204, M206, L225, L228, V232, V269, L270, I308, E312, P313, L315, P316, I317, W318, M403, D426, F429',
     'pore2: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, K138, I162, S196, S199, S200, G203, M204, M206, L225, L228, V232, V269, L270, I308, W318, M319, E321, T322, R326, K327, W328, Q329, L330, F334, M403, D426, F429',
     'pore3: I22, I25, V26, A29, L30, D33, N34, L37, V41, I44, K138, I162, S196, S199, S200, G203, M204, M206, L225, L228, V232, V269, L270, I308, P316, I317, W318, M319, Q329, L330, A333, F334, L384, M403, D426, F429'],
    [],
    ['pore0: R17, I22, I25, V26, A29, L30, N34, L37, K122, L124, F135, P159, I162, S199, S200, G203, M204, M206, L225, L228, V232, V269, L270, E312, P313, L315, P316, I317, W318, F334, M403, D426, F429, Y433',
     'pore1: R17, I22, I25, V26, A29, L30, N34, L37, D121, K122, L124, F135, P159, I162, S199, S200, G203, M204, M206, L225, L228, V232, V269, L270, E312, P313, L315, P316, I317, W318, M319, F334, M403, D426, F429, Y433',
     'pore2: L21, F24, I25, L28, A29, L30, D33, N34, L37, S119, E120, D121, K122, L124, F135, K138, Q142, F166, M169, S196, S199, S200, L228, V232, D262, I265, Q266, V269, E312, P313, L315, P316, I317, W318, F334, D426, F429, Y433',
     'pore3: R17, I22, I25, V26, A29, L30, N34, L37, S119, E120, D121, K122, L124, F135, P159, I162, S199, S200, G203, M204, M206, L225, L228, V232, V269, L270, E312, P313, L315, P316, I317, W318, F334, M403, D426, F429, Y433',
     'pore4: L21, F24, I25, L28, A29, L30, D33, N34, L37, D121, K122, L124, F135, K138, Q142, F166, M169, S196, S199, S200, L228, V232, D262, I265, Q266, V269, E312, P313, L315, P316, I317, W318, F334, D426, F429, Y433',
     'pore5: R17, I22, I25, V26, A29, L30, N34, L37, D121, K122, L124, F135, P159, I162, S199, S200, G203, M204, M206, L225, L228, V232, V269, L270, E312, P313, L315, P316, I317, W318, F334, M403, D426, F429, Y433',
     'pore6: R17, I22, I25, V26, A29, L30, N34, L37, S119, D121, K122, L124, F135, P159, I162, S199, S200, G203, M204, M206, L225, L228, V232, V269, L270, E312, P313, L315, P316, I317, W318, M319, M320, F334, M403, D426, F429, Y433',
     'pore7: R17, I22, I25, V26, A29, L30, N34, L37, V41, P42, I44, F135, P159, I162, S199, S200, G203, M204, M206, L225, L228, V232, V269, L270, E312, W318, M319, E321, T322, M323, R326, K327, W328, Q329, L330, F334, M403, D426, F429, Y433',
     'pore8: R17, I22, I25, V26, A29, L30, N34, L37, P159, I162, S199, S200, G203, M204, M206, L225, L228, V232, V269, L270, G307, I308, L311, A337, Y341, V368, S371, I372, C374, I375, F377, P387, N388, V391, G392, M403, D426, F429, Y433',
     'pore9: R17, I22, I25, V26, A29, L30, N34, L37, P159, I162, S199, S200, G203, M204, M206, L225, L228, V232, V269, L270, G307, I308, L311, A337, Y341, V368, V370, S371, I372, C374, I375, P387, N388, V391, G392, M403, D426, F429, Y433'],
    ['pore0: I22, I25, V26, A29, L30, D33, N34, L37, I44, P45, Y47, N128, V131, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L225, L228, V232, L270, Q276, P277, W318, Y341, I395, D399, M403, D426, F429, Y433',
     'pore1: I22, I25, V26, A29, L30, D33, N34, L37, V41, P42, I44, S46, S119, E120, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L225, L228, V232, L270, Q276, I317, W318, M319, M320, T322, L330, F334, Y341, I395, D399, M403, D426, F429, Y433',
     'pore2: I22, I25, V26, A29, L30, D33, N34, L37, V41, P42, I44, S46, S119, E120, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L225, L228, V232, L270, Q276, P277, I317, W318, M319, M320, T322, L330, F334, Y341, I395, D399, M403, D426, F429, Y433',
     'pore3: I22, I25, V26, A29, L30, D33, N34, L37, S119, E120, D121, K122, L124, E127, N128, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L225, L228, V232, L270, Q276, P313, W318, Y341, I395, D399, M403, D426, F429, Y433',
     'pore4: I22, I25, V26, A29, L30, D33, N34, L37, S119, E120, D121, K122, L124, E127, N128, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L225, L228, V232, L270, Q276, P277, P313, W318, Y341, I395, D399, M403, D426, F429, Y433',
     'pore5: I22, I25, V26, A29, L30, D33, N34, L37, V41, P42, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L225, L228, V232, L270, Q276, I317, W318, M319, E321, T322, R326, K327, W328, Q329, L330, F334, Y341, I395, D399, M403, D426, F429, Y433',
     'pore6: I22, I25, V26, A29, L30, D33, N34, L37, V41, P42, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L225, L228, V232, L270, Q276, P277, I317, W318, M319, E321, T322, R326, K327, W328, Q329, L330, F334, Y341, I395, D399, M403, D426, F429, Y433'],
    [],
    ['pore0: L30, L37, V41, K122, L124, E127, N128, V131, F135, Q142, I149, G150, T153, N154, Y158, S200, V201, A202, M204, G205, L225, L228, V232, G282, T283, L285, E312, P313, A314, L315, P316, I317, W318, F334, M403, V417, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore1: L30, L37, V41, K122, L124, E127, N128, V131, F135, Q142, I149, G150, T153, N154, Y158, S200, V201, A202, M204, G205, S209, L225, L228, V232, Q280, G282, L285, E312, P313, A314, L315, P316, I317, W318, F334, M403, S416, V417, Y418, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore2: L30, L37, V41, K122, L124, E127, N128, V131, F135, Q142, I149, G150, T153, N154, G157, Y158, S200, V201, A202, M204, G205, S209, L225, L228, V232, Q276, S279, Q280, G282, L285, E312, P313, A314, L315, P316, I317, W318, F334, M403, V417, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore3: L30, L37, V41, K122, L124, E127, N128, V131, F135, Q142, I149, G150, T153, N154, G157, Y158, S200, V201, A202, M204, G205, L225, L228, V232, Q276, S279, Q280, G282, L285, E312, P313, A314, L315, P316, I317, W318, F334, M403, V417, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore4: L30, L37, V41, K122, L124, E127, N128, V131, F135, Q142, I149, G150, T153, N154, Y158, S200, V201, A202, M204, G205, L225, L228, V232, S279, Q280, K281, G282, L285, E312, P313, A314, L315, P316, I317, W318, F334, M403, V417, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore5: L30, L37, V41, K122, L124, E127, N128, V131, F135, Q142, I149, G150, T153, N154, I156, G157, Y158, S200, V201, A202, M204, G205, L225, L228, V232, Q276, E278, S279, Q280, G282, L285, E312, P313, A314, L315, P316, I317, W318, F334, M403, V417, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore6: I22, I25, V26, A29, L30, D33, N34, L37, V41, K122, L124, E127, N128, V131, F135, Q142, P159, I162, F163, S196, S199, S200, G203, M204, M206, L225, L228, V232, Q266, L270, Q271, P272, R274, E312, P313, A314, L315, P316, I317, W318, F334, D426, F429, Y433',
     'pore7: L30, L37, V41, S119, E120, D121, K122, L124, E127, N128, V131, F135, Q142, I149, G150, T153, N154, Y158, S200, V201, A202, M204, G205, L225, L228, V232, G282, T283, L285, E312, P313, I317, W318, M319, F334, M403, V417, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore8: L30, L37, V41, S119, E120, D121, K122, L124, E127, N128, V131, F135, Q142, I149, G150, T153, N154, Y158, S200, V201, A202, M204, G205, S209, L225, L228, V232, Q280, G282, L285, E312, P313, I317, W318, M319, F334, M403, S416, V417, Y418, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore9: L30, L37, V41, S119, E120, D121, K122, L124, E127, N128, V131, F135, Q142, I149, G150, T153, N154, G157, Y158, S200, V201, A202, M204, G205, S209, L225, L228, V232, Q276, S279, Q280, G282, L285, E312, P313, I317, W318, M319, F334, M403, V417, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore10: L30, L37, V41, S119, E120, D121, K122, L124, E127, N128, V131, F135, Q142, I149, G150, T153, N154, G157, Y158, S200, V201, A202, M204, G205, L225, L228, V232, Q276, S279, Q280, G282, L285, E312, P313, I317, W318, M319, F334, M403, V417, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore11: L30, L37, V41, S119, E120, D121, K122, L124, E127, N128, V131, F135, Q142, I149, G150, T153, N154, Y158, S200, V201, A202, M204, G205, L225, L228, V232, S279, Q280, K281, G282, L285, E312, P313, I317, W318, M319, F334, M403, V417, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore12: L30, L37, V41, S119, E120, D121, K122, L124, E127, N128, V131, F135, Q142, I149, G150, T153, N154, I156, G157, Y158, S200, V201, A202, M204, G205, L225, L228, V232, Q276, E278, S279, Q280, G282, L285, E312, P313, I317, W318, M319, F334, M403, V417, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore13: I22, I25, V26, A29, L30, D33, N34, L37, V41, S119, E120, D121, K122, L124, E127, N128, V131, F135, Q142, P159, I162, F163, S196, S199, S200, G203, M204, M206, L225, L228, V232, Q266, L270, Q271, P272, R274, E312, P313, I317, W318, M319, F334, D426, F429, Y433',
     'pore14: L30, L37, V41, S119, E120, D121, K122, L124, E127, N128, V131, F135, Q142, I149, G150, T153, N154, Y158, S200, V201, A202, M204, G205, L225, L228, V232, Q280, G282, T283, L285, T287, L288, D291, I294, E312, P313, I317, W318, M319, F334, M403, V410, H414, V417, Y418, G419, S420, V421, Y422, A423, A425, D426, F429, Y433',
     'pore15: L30, L37, V41, S119, E120, D121, K122, L124, E127, N128, V131, F135, Q142, I149, G150, T153, N154, Y158, S200, V201, A202, M204, G205, L225, L228, V232, Q280, G282, T283, L285, T287, L288, L289, K290, D291, I294, E312, P313, I317, W318, M319, F334, M403, V410, H414, V417, Y418, G419, S420, V421, Y422, A423, A425, D426, F429, Y433',
     'pore16: L30, L37, V41, K122, L124, E127, N128, V131, F135, Q142, I149, G150, T153, N154, Y158, S200, V201, A202, M204, G205, L225, L228, V232, G282, T283, L285, E312, P313, P316, I317, W318, M319, F334, M403, V417, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore17: L30, L37, V41, K122, L124, E127, N128, V131, F135, Q142, I149, G150, T153, N154, Y158, S200, V201, A202, M204, G205, S209, L225, L228, V232, Q280, G282, L285, E312, P313, P316, I317, W318, M319, F334, M403, S416, V417, Y418, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore18: L30, L37, V41, K122, L124, E127, N128, V131, F135, Q142, I149, G150, T153, N154, G157, Y158, S200, V201, A202, M204, G205, S209, L225, L228, V232, Q276, S279, Q280, G282, L285, E312, P313, P316, I317, W318, M319, F334, M403, V417, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore19: L30, L37, V41, K122, L124, E127, N128, V131, F135, Q142, I149, G150, T153, N154, G157, Y158, S200, V201, A202, M204, G205, L225, L228, V232, Q276, S279, Q280, G282, L285, E312, P313, P316, I317, W318, M319, F334, M403, V417, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore20: L30, L37, V41, K122, L124, E127, N128, V131, F135, Q142, I149, G150, T153, N154, Y158, S200, V201, A202, M204, G205, L225, L228, V232, S279, Q280, K281, G282, L285, E312, P313, P316, I317, W318, M319, F334, M403, V417, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore21: L30, L37, V41, K122, L124, E127, N128, V131, F135, Q142, I149, G150, T153, N154, I156, G157, Y158, S200, V201, A202, M204, G205, L225, L228, V232, Q276, E278, S279, Q280, G282, L285, E312, P313, P316, I317, W318, M319, F334, M403, V417, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore22: I22, I25, V26, A29, L30, D33, N34, L37, V41, K122, L124, E127, N128, V131, F135, Q142, P159, I162, F163, S196, S199, S200, G203, M204, M206, L225, L228, V232, Q266, L270, Q271, P272, R274, E312, P313, P316, I317, W318, M319, F334, D426, F429, Y433',
     'pore23: L30, L37, V41, E120, D121, E127, N128, V131, F135, Q142, I149, G150, T153, N154, Y158, S200, V201, A202, M204, G205, S209, L225, L228, V232, Q280, G282, L285, E312, P313, I317, W318, F334, M403, S416, V417, Y418, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore24: L30, L37, V41, E120, D121, E127, N128, V131, F135, Q142, I149, G150, T153, N154, G157, Y158, S200, V201, A202, M204, G205, S209, L225, L228, V232, Q276, S279, Q280, G282, L285, E312, P313, I317, W318, F334, M403, V417, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore25: L30, L37, V41, E120, D121, E127, N128, V131, F135, Q142, I149, G150, T153, N154, G157, Y158, S200, V201, A202, M204, G205, L225, L228, V232, Q276, S279, Q280, G282, L285, E312, P313, I317, W318, F334, M403, V417, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore26: L30, L37, V41, E120, D121, E127, N128, V131, F135, Q142, I149, G150, T153, N154, Y158, S200, V201, A202, M204, G205, L225, L228, V232, S279, Q280, K281, G282, L285, E312, P313, I317, W318, F334, M403, V417, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore27: L30, L37, V41, E120, D121, E127, N128, V131, F135, Q142, I149, G150, T153, N154, I156, G157, Y158, S200, V201, A202, M204, G205, L225, L228, V232, Q276, E278, S279, Q280, G282, L285, E312, P313, I317, W318, F334, M403, V417, G419, S420, Y422, A423, A425, D426, F429, Y433',
     'pore28: I22, I25, V26, A29, L30, D33, N34, L37, V41, E120, D121, E127, N128, V131, F135, Q142, P159, I162, F163, S196, S199, S200, G203, M204, M206, L225, L228, V232, Q266, L270, Q271, P272, R274, E312, P313, I317, W318, F334, D426, F429, Y433'],
    [],
    ['pore0: R17, I22, I25, V26, A29, L30, N34, L37, T38, V41, K122, D123, L124, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, M221, A224, L225, L228, V232, L270, Q276, E312, P313, A314, L315, P316, I317, W318, M403, P404, D426, F429, Y433',
     'pore1: I22, I25, V26, A29, L30, N34, L37, T38, V41, K122, D123, L124, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, V210, M221, A224, L225, L228, V232, L270, Q276, E312, P313, A314, L315, P316, I317, W318, M403, P404, D426, F429, Y433',
     'pore2: I22, I25, V26, A29, L30, N34, L37, T38, V41, K122, D123, L124, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, S209, V210, M221, A224, L225, L228, V232, L270, Q276, P277, E312, P313, A314, L315, P316, I317, W318, M403, P404, D426, F429, Y433',
     'pore3: I22, I25, V26, A29, L30, N34, L37, T38, V41, K122, D123, L124, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, M221, A224, L225, L228, V232, L270, Q271, R274, V275, Q276, E312, P313, A314, L315, P316, I317, W318, M403, P404, D426, F429, Y433',
     'pore4: I22, I25, V26, A29, L30, N34, L37, T38, V41, K122, D123, L124, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, M221, A224, L225, L228, V232, L270, Q271, P272, R274, Q276, E312, P313, A314, L315, P316, I317, W318, M403, P404, D426, F429, Y433',
     'pore5: I22, I25, V26, A29, L30, N34, L37, T38, V41, K122, D123, L124, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, M221, A224, L225, L228, V232, L270, Q271, P272, R274, Q276, E312, P313, A314, L315, P316, I317, W318, M403, P404, D426, F429, Y433',
     'pore6: L30, N34, L37, T38, V41, K122, D123, L124, N128, V131, F135, K138, Q142, S196, S200, M204, A208, Y211, T212, D214, E215, R217, G218, M221, G222, L225, L228, V232, E312, P313, A314, L315, P316, I317, W318, Y341, T345, G349, R357, I395, G396, D399, S400, M403, P404, G407, Y408, D411, Y418, Y422, D426, F429, Y433',
     'pore7: R17, I22, I25, V26, A29, L30, N34, L37, T38, V41, S119, E120, D121, K122, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, M221, A224, L225, L228, V232, L270, Q276, E312, P313, L315, P316, I317, W318, M403, P404, D426, F429, Y433',
     'pore8: I22, I25, V26, A29, L30, N34, L37, T38, V41, S119, E120, D121, K122, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, V210, M221, A224, L225, L228, V232, L270, Q276, E312, P313, L315, P316, I317, W318, M403, P404, D426, F429, Y433',
     'pore9: I22, I25, V26, A29, L30, N34, L37, T38, V41, S119, E120, D121, K122, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, S209, V210, M221, A224, L225, L228, V232, L270, Q276, P277, E312, P313, L315, P316, I317, W318, M403, P404, D426, F429, Y433',
     'pore10: I22, I25, V26, A29, L30, N34, L37, T38, V41, S119, E120, D121, K122, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, M221, A224, L225, L228, V232, L270, Q271, R274, V275, Q276, E312, P313, L315, P316, I317, W318, M403, P404, D426, F429, Y433',
     'pore11: I22, I25, V26, A29, L30, N34, L37, T38, V41, S119, E120, D121, K122, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, M221, A224, L225, L228, V232, L270, Q271, P272, R274, Q276, E312, P313, L315, P316, I317, W318, M403, P404, D426, F429, Y433',
     'pore12: I22, I25, V26, A29, L30, N34, L37, T38, V41, S119, E120, D121, K122, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, M221, A224, L225, L228, V232, L270, Q271, P272, R274, Q276, E312, P313, L315, P316, I317, W318, M403, P404, D426, F429, Y433',
     'pore13: L30, N34, L37, T38, V41, S119, E120, D121, K122, N128, V131, F135, K138, Q142, S196, S200, M204, A208, Y211, T212, D214, E215, R217, G218, M221, G222, L225, L228, V232, E312, P313, L315, P316, I317, W318, Y341, T345, G349, R357, I395, G396, D399, S400, M403, P404, G407, Y408, D411, Y418, Y422, D426, F429, Y433',
     'pore14: R17, I22, I25, V26, A29, L30, N34, L37, T38, V41, K122, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, M221, A224, L225, L228, V232, L270, Q276, E312, P313, L315, P316, I317, W318, M319, M403, P404, D426, F429, Y433',
     'pore15: I22, I25, V26, A29, L30, N34, L37, T38, V41, K122, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, V210, M221, A224, L225, L228, V232, L270, Q276, E312, P313, L315, P316, I317, W318, M319, M403, P404, D426, F429, Y433',
     'pore16: I22, I25, V26, A29, L30, N34, L37, T38, V41, K122, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, S209, V210, M221, A224, L225, L228, V232, L270, Q276, P277, E312, P313, L315, P316, I317, W318, M319, M403, P404, D426, F429, Y433',
     'pore17: I22, I25, V26, A29, L30, N34, L37, T38, V41, K122, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, M221, A224, L225, L228, V232, L270, Q271, R274, V275, Q276, E312, P313, L315, P316, I317, W318, M319, M403, P404, D426, F429, Y433',
     'pore18: I22, I25, V26, A29, L30, N34, L37, T38, V41, K122, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, M221, A224, L225, L228, V232, L270, Q271, P272, R274, Q276, E312, P313, L315, P316, I317, W318, M319, M403, P404, D426, F429, Y433',
     'pore19: I22, I25, V26, A29, L30, N34, L37, T38, V41, K122, N128, V131, F135, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, M221, A224, L225, L228, V232, L270, Q271, P272, R274, Q276, E312, P313, L315, P316, I317, W318, M319, M403, P404, D426, F429, Y433',
     'pore20: L30, N34, L37, T38, V41, K122, N128, V131, F135, K138, Q142, S196, S200, M204, A208, Y211, T212, D214, E215, R217, G218, M221, G222, L225, L228, V232, E312, P313, L315, P316, I317, W318, M319, Y341, T345, G349, R357, I395, G396, D399, S400, M403, P404, G407, Y408, D411, Y418, Y422, D426, F429, Y433',
     'pore21: R17, I22, I25, V26, A29, L30, N34, L37, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, M221, A224, L225, L228, V232, L270, Q276, I308, E312, L315, P316, I317, M319, Q329, A333, F334, A337, S338, Y341, I381, L384, I385, N388, M403, P404, D426, F429, Y433',
     'pore22: I22, I25, V26, A29, L30, N34, L37, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, V210, M221, A224, L225, L228, V232, L270, Q276, I308, E312, L315, P316, I317, M319, Q329, A333, F334, A337, S338, Y341, I381, L384, I385, N388, M403, P404, D426, F429, Y433',
     'pore23: I22, I25, V26, A29, L30, N34, L37, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, S209, V210, M221, A224, L225, L228, V232, L270, Q276, P277, I308, E312, L315, P316, I317, M319, Q329, A333, F334, A337, S338, Y341, I381, L384, I385, N388, M403, P404, D426, F429, Y433',
     'pore24: I22, I25, V26, A29, L30, N34, L37, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, M221, A224, L225, L228, V232, L270, Q271, R274, V275, Q276, I308, E312, L315, P316, I317, M319, Q329, A333, F334, A337, S338, Y341, I381, L384, I385, N388, M403, P404, D426, F429, Y433',
     'pore25: I22, I25, V26, A29, L30, N34, L37, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, M221, A224, L225, L228, V232, L270, Q271, P272, R274, Q276, I308, E312, L315, P316, I317, M319, Q329, A333, F334, A337, S338, Y341, I381, L384, I385, N388, M403, P404, D426, F429, Y433',
     'pore26: I22, I25, V26, A29, L30, N34, L37, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, M221, A224, L225, L228, V232, L270, Q271, P272, R274, Q276, I308, E312, L315, P316, I317, M319, Q329, A333, F334, A337, S338, Y341, I381, L384, I385, N388, M403, P404, D426, F429, Y433',
     'pore27: L30, N34, L37, K138, Q142, S196, S200, M204, A208, Y211, T212, D214, E215, R217, G218, M221, G222, L225, L228, V232, I308, E312, L315, P316, I317, M319, Q329, A333, F334, A337, S338, Y341, T345, G349, R357, I381, L384, I385, N388, I395, G396, D399, S400, M403, P404, G407, Y408, D411, Y418, Y422, D426, F429, Y433',
     'pore28: R17, I22, I25, V26, A29, L30, N34, L37, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, M221, A224, L225, L228, V232, L270, Q276, I308, E312, L315, P316, I317, A333, F334, A337, S338, Y341, I381, L384, I385, N388, M403, P404, D426, F429, Y433',
     'pore29: I22, I25, V26, A29, L30, N34, L37, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, V210, M221, A224, L225, L228, V232, L270, Q276, I308, E312, L315, P316, I317, A333, F334, A337, S338, Y341, I381, L384, I385, N388, M403, P404, D426, F429, Y433',
     'pore30: I22, I25, V26, A29, L30, N34, L37, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, S209, V210, M221, A224, L225, L228, V232, L270, Q276, P277, I308, E312, L315, P316, I317, A333, F334, A337, S338, Y341, I381, L384, I385, N388, M403, P404, D426, F429, Y433',
     'pore31: I22, I25, V26, A29, L30, N34, L37, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, M221, A224, L225, L228, V232, L270, Q271, R274, V275, Q276, I308, E312, L315, P316, I317, A333, F334, A337, S338, Y341, I381, L384, I385, N388, M403, P404, D426, F429, Y433',
     'pore32: I22, I25, V26, A29, L30, N34, L37, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, M221, A224, L225, L228, V232, L270, Q271, P272, R274, Q276, I308, E312, L315, P316, I317, A333, F334, A337, S338, Y341, I381, L384, I385, N388, M403, P404, D426, F429, Y433',
     'pore33: I22, I25, V26, A29, L30, N34, L37, K138, Q142, Y158, P159, I162, S196, S199, S200, G203, M204, M206, L207, M221, A224, L225, L228, V232, L270, Q271, P272, R274, Q276, I308, E312, L315, P316, I317, A333, F334, A337, S338, Y341, I381, L384, I385, N388, M403, P404, D426, F429, Y433',
     'pore34: L30, N34, L37, K138, Q142, S196, S200, M204, A208, Y211, T212, D214, E215, R217, G218, M221, G222, L225, L228, V232, I308, E312, L315, P316, I317, A333, F334, A337, S338, Y341, T345, G349, R357, I381, L384, I385, N388, I395, G396, D399, S400, M403, P404, G407, Y408, D411, Y418, Y422, D426, F429, Y433'],
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
   @> pore 1: 	1094.62 		79.18 		0.83
   @> pore 2: 	1071.63 		76.13 		0.83
   @> pore 3: 	1085.25 		78.72 		0.83
   @> Frame/model: 7
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	959.56 		65.08 		0.91
   @> pore 1: 	872.05 		69.54 		0.81
   @> pore 2: 	917.19 		65.01 		0.91
   @> pore 3: 	969.25 		67.51 		0.91
   @> pore 4: 	881.74 		71.97 		0.81
   @> pore 5: 	821.42 		71.65 		0.81
   @> pore 6: 	950.45 		67.11 		0.89
   @> pore 7: 	914.65 		75.02 		0.81
   @> pore 8: 	862.94 		71.57 		0.81
   @> pore 9: 	891.37 		77.4 		0.81
   @> pore 10: 	839.66 		73.96 		0.81
   @> Frame/model: 8
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	1285.2 		71.85 		1.01
   @> Frame/model: 9
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 10
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	1272.43 		72.75 		0.99
   @> pore 1: 	1247.9 		72.69 		1.05
   @> pore 2: 	1102.47 		73.66 		0.9
   @> pore 3: 	1131.07 		74.98 		0.92
   @> Frame/model: 11
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 12
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	885.93 		68.25 		0.8
   @> pore 1: 	911.8 		72.55 		0.8
   @> pore 2: 	936.28 		74.49 		0.98
   @> pore 3: 	931.71 		74.82 		0.8
   @> pore 4: 	912.02 		72.71 		0.85
   @> pore 5: 	907.45 		73.05 		0.8
   @> pore 6: 	948.39 		76.4 		0.8
   @> pore 7: 	848.04 		74.38 		0.71
   @> pore 8: 	685.75 		74.95 		0.7
   @> pore 9: 	697.44 		75.87 		0.7
   @> Frame/model: 13
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	915.86 		73.24 		0.86
   @> pore 1: 	973.99 		80.09 		0.86
   @> pore 2: 	993.16 		81.56 		0.86
   @> pore 3: 	934.36 		75.9 		0.83
   @> pore 4: 	953.54 		77.37 		0.83
   @> pore 5: 	974.12 		82.27 		0.81
   @> pore 6: 	993.3 		83.74 		0.81
   @> Frame/model: 14
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 15
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	867.57 		65.62 		0.89
   @> pore 1: 	915.06 		71.53 		0.89
   @> pore 2: 	919.31 		71.94 		0.89
   @> pore 3: 	935.51 		73.14 		0.89
   @> pore 4: 	872.24 		69.15 		0.88
   @> pore 5: 	933.88 		75.21 		0.85
   @> pore 6: 	886.71 		75.42 		0.81
   @> pore 7: 	860.05 		67.69 		0.89
   @> pore 8: 	907.53 		73.59 		0.89
   @> pore 9: 	911.78 		74.01 		0.89
   @> pore 10: 	927.99 		75.2 		0.89
   @> pore 11: 	864.72 		71.21 		0.88
   @> pore 12: 	926.36 		77.27 		0.85
   @> pore 13: 	879.19 		77.48 		0.81
   @> pore 14: 	895.72 		80.73 		0.89
   @> pore 15: 	906.85 		86.35 		0.73
   @> pore 16: 	840.03 		67.91 		0.89
   @> pore 17: 	887.52 		73.82 		0.89
   @> pore 18: 	891.77 		74.24 		0.89
   @> pore 19: 	907.98 		75.43 		0.89
   @> pore 20: 	844.7 		71.44 		0.88
   @> pore 21: 	906.35 		77.5 		0.85
   @> pore 22: 	859.17 		77.71 		0.81
   @> pore 23: 	869.96 		74.3 		0.89
   @> pore 24: 	874.21 		74.71 		0.89
   @> pore 25: 	890.41 		75.9 		0.89
   @> pore 26: 	827.14 		71.92 		0.88
   @> pore 27: 	888.78 		77.98 		0.85
   @> pore 28: 	841.61 		78.19 		0.81
   @> Frame/model: 16
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> Frame/model: 17
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	854.46 		84.52 		0.7
   @> pore 1: 	832.03 		83.19 		0.7
   @> pore 2: 	894.84 		85.61 		0.7
   @> pore 3: 	851.62 		87.64 		0.7
   @> pore 4: 	866.6 		89.58 		0.7
   @> pore 5: 	873.57 		90.01 		0.7
   @> pore 6: 	878.32 		86.75 		0.84
   @> pore 7: 	843.84 		87.66 		0.7
   @> pore 8: 	821.41 		86.32 		0.7
   @> pore 9: 	884.22 		88.75 		0.7
   @> pore 10: 	841.0 		90.78 		0.7
   @> pore 11: 	855.98 		92.72 		0.7
   @> pore 12: 	862.94 		93.15 		0.7
   @> pore 13: 	867.69 		89.89 		0.84
   @> pore 14: 	831.7 		86.94 		0.7
   @> pore 15: 	809.27 		85.61 		0.7
   @> pore 16: 	872.08 		88.03 		0.7
   @> pore 17: 	828.86 		90.06 		0.7
   @> pore 18: 	843.84 		92.01 		0.7
   @> pore 19: 	850.81 		92.43 		0.7
   @> pore 20: 	855.56 		89.17 		0.84
   @> pore 21: 	722.32 		83.19 		0.7
   @> pore 22: 	699.89 		81.85 		0.7
   @> pore 23: 	762.7 		84.28 		0.7
   @> pore 24: 	719.48 		86.3 		0.7
   @> pore 25: 	734.46 		88.25 		0.7
   @> pore 26: 	741.42 		88.68 		0.7
   @> pore 27: 	746.17 		85.42 		0.84
   @> pore 28: 	701.63 		83.13 		0.7
   @> pore 29: 	679.2 		81.79 		0.7
   @> pore 30: 	742.01 		84.22 		0.7
   @> pore 31: 	698.79 		86.24 		0.7
   @> pore 32: 	713.77 		88.19 		0.7
   @> pore 33: 	720.74 		88.62 		0.7
   @> pore 34: 	725.49 		85.35 		0.84
   @> Frame/model: 18
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	880.75 		71.11 		0.79
   @> pore 1: 	883.86 		71.9 		0.79
   @> pore 2: 	903.47 		84.71 		0.79
   @> pore 3: 	811.74 		72.54 		0.79
   @> pore 4: 	704.54 		64.47 		0.74
   @> pore 5: 	719.6 		65.67 		0.74
   @> pore 6: 	820.11 		86.68 		0.79
   @> pore 7: 	781.35 		65.67 		0.74
   @> pore 8: 	796.42 		66.87 		0.74
   @> pore 9: 	896.92 		87.87 		0.79
   @> pore 10: 	774.87 		65.48 		0.74
   @> pore 11: 	789.93 		66.68 		0.74
   @> Frame/model: 19
   @> Pore ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> pore 0: 	917.75 		66.68 		0.73
   @> pore 1: 	903.64 		68.75 		0.73
   @> pore 2: 	1005.38 		79.53 		0.75
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
    ([65.49187390762492,
      69.74914632506119,
      71.52262812400917,
      72.7733284937729,
      70.8007142570642,
      75.063257094092,
      76.83213948274864,
      78.08606293309248,
      72.68178362741891,
      73.07368296067034,
      76.94336840869165,
      78.71204702710276,
      79.96616523324155,
      73.64851869588655,
      74.0404187950495,
      77.91010620767416,
      79.67878339869065,
      80.93290273764383,
      79.22759984948198,
      79.61750031295861,
      83.48296968080439,
      85.25186248855925,
      86.50517537270744],
     [0.748030901015996,
      0.748030901015996,
      0.748030901015996,
      0.6926669440694282,
      0.748030901015996,
      0.748030901015996,
      0.748030901015996,
      0.6926669440694282,
      0.748030901015996,
      0.748030901015996,
      0.748030901015996,
      0.748030901015996,
      0.6926669440694282,
      0.7458901659707238,
      0.7458901659707238,
      0.7458901659707238,
      0.7458901659707238,
      0.6926669440694282,
      0.748030901015996,
      0.748030901015996,
      0.748030901015996,
      0.748030901015996,
      0.6926669440694282],
     [825.7290019086264,
      834.7235555986805,
      879.7191876465255,
      843.1221951617829,
      823.6319404487533,
      832.6264941388073,
      877.6221261866525,
      841.0251337019097,
      805.4632114742075,
      800.6942671386256,
      814.4577651642614,
      859.4533972121067,
      822.856404727364,
      809.2981737540832,
      804.5292294185012,
      818.2927274441372,
      863.2883594919823,
      826.6913670072396,
      839.0404648719796,
      834.2715205363977,
      848.0350185620338,
      893.0306506098788,
      856.4336581251359]),
    ([76.59035345440319, 79.18162358293792, 76.1291573504211, 78.72016935667617],
     [0.8298276094079913,
      0.8298276094079913,
      0.8298276094079913,
      0.8298276094079913],
     [1081.0009523072954,
      1094.62026833373,
      1071.6286171659622,
      1085.2479331923969]),
    ([65.07866715066874,
      69.53808736123777,
      65.01461543950278,
      67.51081282801462,
      71.96971734263597,
      71.64593491645554,
      67.10991973600689,
      75.01547272653391,
      71.5698914487432,
      77.40246995375526,
      73.95688529731291],
     [0.9118922696158223,
      0.813576645807432,
      0.9118922696158223,
      0.9118922696158223,
      0.813576645807432,
      0.813576645807432,
      0.8935285215508568,
      0.813576645807432,
      0.813576645807432,
      0.813576645807432,
      0.813576645807432],
     [959.5591530967635,
      872.0499867930597,
      917.1902789209371,
      969.2489193570057,
      881.7397530533019,
      821.4237084754667,
      950.4500753030867,
      914.6500818581253,
      862.9409089993828,
      891.3678498747569,
      839.6586770160143]),
    ([71.84677180411046], [1.0104558297332396], [1285.2009039552493]),
    ([], [], []),
    ([72.74712788967892, 72.68719290556788, 73.66058180257426, 74.98192073087397],
     [0.9935284974640439,
      1.052007326607179,
      0.8999179901412739,
      0.9195967311006099],
     [1272.4266446468707,
      1247.902393132334,
      1102.4735112692838,
      1131.065926597962]),
    ([], [], []),
    ([68.25269580315992,
      72.55315633184057,
      74.48644504215622,
      74.81917459304985,
      72.71372399926511,
      73.04640285902778,
      76.40333224095926,
      74.37986132343786,
      74.95009729418638,
      75.87021679201985],
     [0.7980325956849226,
      0.7980325956849226,
      0.9771393588924104,
      0.7980325956849226,
      0.8506848374545151,
      0.7980325956849226,
      0.7980325956849226,
      0.7145541112891076,
      0.701643210667199,
      0.701643210667199],
     [885.9301447013927,
      911.797377583824,
      936.2797810091117,
      931.70828572543,
      912.0183345668365,
      907.446839283155,
      948.3943343326728,
      848.0380047674669,
      685.7489814869723,
      697.4409164999504]),
    ([73.23540381968851,
      80.08939091267527,
      81.55516468299065,
      75.90001305295283,
      77.36578682326822,
      82.27318347367724,
      83.73884127245586],
     [0.8590086239058461,
      0.8590086239058461,
      0.8590086239058461,
      0.8306002979099667,
      0.8306002979099667,
      0.808705875839397,
      0.808705875839397],
     [915.8631574849971,
      973.9858988745997,
      993.1616551551109,
      934.3646830086334,
      953.5404392891448,
      974.1215949224061,
      993.2973512029174]),
    ([], [], []),
    ([65.62115037402914,
      71.53042500169187,
      71.94495993587996,
      73.13596764052856,
      69.14935644570515,
      75.21057053015521,
      75.41852108658443,
      67.68749183474073,
      73.59341610395667,
      74.00807257054497,
      75.19831628499098,
      71.21260520228009,
      77.27285604019785,
      77.48354012077928,
      80.73116878743019,
      86.35071464854181,
      67.91311324028494,
      73.82237230230253,
      74.23690968717791,
      75.42791045208884,
      71.44130747880759,
      77.50251229790055,
      77.71047620998692,
      74.29603275648617,
      74.71066327149282,
      75.90098579161851,
      71.91518306799668,
      77.97553892502322,
      78.18608012550663],
     [0.8880679059215846,
      0.8880679059215846,
      0.8880679059215846,
      0.8880679059215846,
      0.8785444091924826,
      0.8490506238172444,
      0.8102864597349899,
      0.8880679059215846,
      0.8880679059215846,
      0.8880679059215846,
      0.8880679059215846,
      0.8785444091924826,
      0.8490506238172444,
      0.8102864597349899,
      0.8880679059215846,
      0.734349938926306,
      0.8880679059215846,
      0.8880679059215846,
      0.8880679059215846,
      0.8880679059215846,
      0.8785444091924826,
      0.8490506238172444,
      0.8102864597349899,
      0.8880679059215846,
      0.8880679059215846,
      0.8880679059215846,
      0.8785444091924826,
      0.8490506238172444,
      0.8102864597349899],
     [867.5691070465471,
      915.0571082621244,
      919.3080573275604,
      935.5118864403646,
      872.2410322484818,
      933.8841442078436,
      886.7097725815132,
      860.045382840555,
      907.5333840561325,
      911.7843331215684,
      927.9881622343727,
      864.7173080424898,
      926.3604200018518,
      879.1860483755213,
      895.7204387925915,
      906.851825264617,
      840.0326274337124,
      887.5206286492897,
      891.7715777147259,
      907.9754068275302,
      844.7045526356471,
      906.3476645950092,
      859.1732929686787,
      869.9566884402092,
      874.2076375056453,
      890.4114666184494,
      827.1406124265665,
      888.7837243859285,
      841.6093527595981]),
    ([], [], []),
    ([84.52472492120339,
      83.18517166215759,
      85.6116243800833,
      87.63585253941989,
      89.58165950651599,
      90.01291582330632,
      86.74710026335319,
      87.66488852869145,
      86.32497048349092,
      88.75179157508275,
      90.777482492146,
      92.72387104704652,
      93.15411195442789,
      89.88776506127502,
      86.9448941641825,
      85.60546913897271,
      88.03179353420344,
      90.05867204300597,
      92.00606004750155,
      92.43477506853809,
      89.16966073325125,
      83.19346479682866,
      81.85396442345473,
      84.28036425570858,
      86.30452142098932,
      88.25033684761284,
      88.68159627492663,
      85.4159351041474,
      83.1305995056015,
      81.79126114672023,
      84.21749669246464,
      86.24113006422118,
      88.18683574527624,
      88.61834382141099,
      85.35294774682225],
     [0.6975364266842847,
      0.6975364266842847,
      0.6975364266842847,
      0.6975364266842847,
      0.6975364266842847,
      0.6975364266842847,
      0.8418977703001918,
      0.6975364266842847,
      0.6975364266842847,
      0.6975364266842847,
      0.6975364266842847,
      0.6975364266842847,
      0.6975364266842847,
      0.8418977703001918,
      0.6975364266842847,
      0.6975364266842847,
      0.6975364266842847,
      0.6975364266842847,
      0.6975364266842847,
      0.6975364266842847,
      0.8418977703001918,
      0.6975364266842847,
      0.6975364266842847,
      0.6975364266842847,
      0.6975364266842847,
      0.6975364266842847,
      0.6975364266842847,
      0.8418977703001918,
      0.6975364266842847,
      0.6975364266842847,
      0.6975364266842847,
      0.6975364266842847,
      0.6975364266842847,
      0.6975364266842847,
      0.8418977703001918],
     [854.4598782279375,
      832.0290655726196,
      894.8389005466951,
      851.621390583104,
      866.6017409708695,
      873.566863623062,
      878.3154754074413,
      843.8364798624227,
      821.4056672071049,
      884.2155021811803,
      840.9979922175892,
      855.9783426053549,
      862.9434652575472,
      867.6920770419265,
      831.7002354906767,
      809.2694228353588,
      872.0792578094344,
      828.8617478458432,
      843.8420982336088,
      850.8072208858013,
      855.5558326701806,
      722.3169456508326,
      699.8861329955148,
      762.6959679695902,
      719.4784580059992,
      734.4588083937647,
      741.4239310459571,
      746.1725428300083,
      701.6311882116455,
      679.2003755563277,
      742.0102105304032,
      698.7927005668121,
      713.7730509545777,
      720.7381736067701,
      725.4867853908213]),
    ([71.10879068607957,
      71.89909986983707,
      84.71038931655131,
      72.53780220393708,
      64.472820243521,
      65.67411620356634,
      86.67805207852878,
      65.6662630451571,
      66.8675625557789,
      87.87179570521266,
      65.48160153745087,
      66.68289049054799],
     [0.7853632391543444,
      0.7853632391543444,
      0.7853632391543444,
      0.7853632391543444,
      0.7391908821133051,
      0.7391908821133051,
      0.7853632391543444,
      0.7391908821133051,
      0.7391908821133051,
      0.7853632391543444,
      0.7391908821133051,
      0.7391908821133051],
     [880.7471463815876,
      883.864224406766,
      903.4710369413804,
      811.7427929294366,
      704.53841218915,
      719.6040544317701,
      820.1088738097437,
      781.3494384396813,
      796.4150806823014,
      896.9199000605231,
      774.8687654631555,
      789.9344077057756]),
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
