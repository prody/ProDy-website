.. _cavitracer_single:

Exploring CaviTracer Parameters for Channels and Surface Cavities
===============================================================================

Choosing suitable parameters for channel and surface-cavity calculations may
depend on the analyzed system, its size, structural resolution, and the type of
pathway or pocket of interest. This tutorial shows how to explore different
CaviTracer parameter combinations in ProDy and compare the resulting channels
and surface cavities.

The examples demonstrate how to run parameter scans, save the detected objects
as PQR files, and generate spatial occupancy maps. These maps help identify
regions that are consistently detected across parameter sets, as well as
channels or surface cavities that appear only for selected parameter values.
This can guide the choice of parameters for further analysis and visualization.


Channels
-------------------------------------------------------------------------------

As an example for this tutorial, we will analyze the structure of cytochrome
P450 which contains 486 residues. To analyze the structure, we need to parse
a structure :file:`1tqn` using :func:`.parsePDB` and select protein
structure:

.. ipython:: python
   :verbatim:

   atoms = parsePDB('1tqn')
   protein = atoms.select('protein')

.. parsed-literal::

   @> Connecting wwPDB FTP server RCSB PDB (USA).
   @> Downloading PDB files via FTP failed, trying HTTP.
   @> 1tqn downloaded (1tqn.pdb.gz)
   @> PDB download via HTTP completed (1 downloaded, 0 failed).
   @> 3999 atoms and 1 coordinate set(s) were parsed in 0.14s.


Next, :func:`.scanChannelParameters` can be use to explore various parameters 
and how they affects channel detection. This is useful when the optimal values are
not known for a given protein or when the user wants to identify channels that
are robustly detected across several parameter combinations.

By default, the function tests all combinations of the following parameters:

- ``inner_radius``: ``1.2``, ``1.4``, and ``1.6`` Å,
- ``sparsity``: ``1.0``, ``3.0``, and ``5.0`` Å,
- ``min_depth``: ``3.0``, ``5.0``, and ``10.0`` Å.

This gives 27 parameter combinations in total. For each combination, channels
are calculated and saved as a PQR file. The function also generates an
occupancy map showing which regions are detected repeatedly across the parameter
scan. Regions with higher occupancy correspond to channels that are less
sensitive to parameter choice, whereas low-occupancy regions indicate channels
detected only under selected parameter settings.

The returned ``channels_all`` object contains the channels detected for each
parameter combination, ``parameter_sets`` stores the corresponding parameter
values, and ``occupancy_file`` is the path to the PDB file containing the spatial
occupancy map.


.. ipython:: python
   :verbatim:

   channels_all, parameter_sets, occupancy_file = scanChannelParameters(protein, 
						output_path='channel_param_grid')


.. parsed-literal::

    @> Calculating channels for 27 parameter combinations.
    @> Grid run 1/27: inner_radius=1.2, sparsity=2, min_depth=3
    @> The atoms supplied to calcChannels contain protein atoms only.
    @> Substituted 3766 atoms with 23638 homogeneous balls of radius 1.52 Å in 0.17s.
    @> Delaunay tessellation of 23638 points constructed in 0.78s.
    @> Surface and inner simplices filtered in 1.64s.
    @> Cavities: 129 found, 16 deeper than min_depth=3.0 Å, 15 of them at least seed_volume=50 Å³ and searched for channels; the 1 smaller ones are tessellation debris and are left unsearched, in 0.31s.
    @> Chambers (probe 1.40 Å): 11 of the 15 searched cavities have them; the other 9 are searched whole.
    @>     cavity 0: 13 chambers, 4 of them seeded.
    @>     cavity 1: 10 chambers, 1 of them seeded.
    @>     cavity 2: 1 chamber, seeded.
    @>     cavity 3: 3 chambers, 1 of them seeded.
    @>     cavity 4: 1 chamber, none of them deep and large enough to seed; searched whole.
    @>     cavity 5: 4 chambers, 1 of them seeded.
    @>     cavity 7: 1 chamber, seeded.
    @>     cavity 8: 1 chamber, none of them deep and large enough to seed; searched whole.
    @>     cavity 9: 1 chamber, none of them deep and large enough to seed; searched whole.
    @>     cavity 10: 2 chambers, none of them deep and large enough to seed; searched whole.
    @>     cavity 11: 1 chamber, none of them deep and large enough to seed; searched whole.
    @> 18 search sites (sp) in 0.08s: one per seeded chamber, one per cavity searched whole.
    @> Channel search (Dijkstra) over 18 search sites in 15 cavities completed in 0.28s.
    @> Found 20 channels and 2 links (a link joins a deep chamber to a shallower one and never reaches the surface).
    @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
    @>     site  start_point [Å]              void                   volume [Å³]  depth [Å]  channels  links
    @>     sp0   [-17.336, -19.734, -11.982]  cavity 0, chamber 1/4         4919       14.9         3      -
    @>     sp1   [-27.873, -18.703, 2.491]    cavity 4, whole                311        3.6         1      -
    @>     sp2   [-11.141, -24.081, 7.425]    cavity 6, whole                275        3.1         1      -
    @>     sp3   [-27.878, -30.525, -18.960]  cavity 0, chamber 2/4          255        3.7         1      -
    @>     sp4   [-13.078, -37.175, -1.924]   cavity 8, whole                247        3.3         1      -
    @>     sp5   [-25.439, -42.211, -10.996]  cavity 2, chamber 1/1          216        6.9         -      -  sealed
    @>     sp6   [-24.989, 1.816, -23.778]    cavity 9, whole                214        3.1         1      -
    @>     sp7   [-30.276, -32.754, -28.375]  cavity 10, whole               184        3.6         2      -
    @>     sp8   [-28.112, -17.893, -23.847]  cavity 3, chamber 1/1          173        3.6         1      -
    @>     sp9   [-15.595, -34.280, 11.278]   cavity 11, whole               148        3.1         1      -
    @>     sp10  [-13.854, -37.671, 5.928]    cavity 12, whole               139        3.4         1      -
    @>     sp11  [-9.276, -31.119, -5.477]    cavity 1, chamber 1/1          120        4.8         1      -
    @>     sp12  [-30.654, -44.562, -24.180]  cavity 13, whole               109        3.0         1      -
    @>     sp13  [-13.624, -43.494, -13.155]  cavity 7, chamber 1/1           97        3.1         1      -
    @>     sp14  [-7.797, -9.475, -4.093]     cavity 14, whole                87        3.8         1      -
    @>     sp15  [-17.760, -11.110, -5.738]   cavity 0, chamber 3/4           69        3.9         1      1  -> sp0
    @>     sp16  [-28.252, -30.840, -12.953]  cavity 0, chamber 4/4           51        3.5         1      1  -> sp3
    @>     sp17  [-27.498, -24.585, 0.741]    cavity 5, chamber 1/1           50        9.3         1      -
    @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
    @> The 1 site marked sealed above report neither a channel nor a link: no route out of them survived - either narrower than bottleneck=1.20 Å, or dropped as a duplicate of a shallower site's, or the void is its own mouth and has nowhere to path to. Lower bottleneck to see how the first kind connect.
    @> Saving 20 channels to channel_param_grid/channels_run000_inner_radius_1p2_sparsity_2_depth_3.pqr and 2 links to channel_param_grid/channels_run000_inner_radius_1p2_sparsity_2_depth_3_links.pqr.
    @> Channel calculation completed in 3.19s.
    @> Grid run 2/27: inner_radius=1.2, sparsity=2, min_depth=5
    @> The atoms supplied to calcChannels contain protein atoms only.
    @> Substituted 3766 atoms with 23638 homogeneous balls of radius 1.52 Å in 0.17s.
    @> Delaunay tessellation of 23638 points constructed in 0.75s.
    @> Surface and inner simplices filtered in 1.65s.
    @> Cavities: 129 found, 7 deeper than min_depth=5.0 Å and searched for channels, in 0.29s.
    @> Chambers (probe 1.40 Å): 6 of the 7 searched cavities have them; the other 4 are searched whole.
    @>     cavity 0: 10 chambers, 2 of them seeded.
    @>     cavity 1: 1 chamber, none of them deep and large enough to seed; searched whole.
    @>     cavity 2: 1 chamber, seeded.
    @>     cavity 3: 1 chamber, none of them deep and large enough to seed; searched whole.
    @>     cavity 4: 3 chambers, 1 of them seeded.
    @>     cavity 5: 1 chamber, none of them deep and large enough to seed; searched whole.
    @> 8 search sites (sp) in 0.06s: one per seeded chamber, one per cavity searched whole.
    @> Channel search (Dijkstra) over 8 search sites in 7 cavities completed in 0.20s.
    @> Found 9 channels.
    @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
    @>     site  start_point [Å]              void                   volume [Å³]  depth [Å]  channels  links
    @>     sp0   [-17.336, -19.734, -11.982]  cavity 0, chamber 1/2         4701       14.9         3      -
    @>     sp1   [-9.260, -31.344, -5.479]    cavity 1, whole               1410        5.0         1      -
    @>     sp2   [-27.710, -17.308, -22.652]  cavity 3, whole                359        5.0         1      -
    @>     sp3   [-15.652, -41.491, -12.372]  cavity 5, whole                256        5.6         1      -
    @>     sp4   [-13.208, -37.795, -3.550]   cavity 6, whole                247        5.0         1      -
    @>     sp5   [-25.439, -42.211, -10.996]  cavity 2, chamber 1/1          199        6.9         -      -  sealed
    @>     sp6   [-26.621, -29.554, -20.717]  cavity 0, chamber 2/2          116        5.1         1      -
    @>     sp7   [-27.498, -24.585, 0.741]    cavity 4, chamber 1/1           50        9.3         1      -
    @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
    @> The 1 site marked sealed above report neither a channel nor a link: no route out of them survived - either narrower than bottleneck=1.20 Å, or dropped as a duplicate of a shallower site's, or the void is its own mouth and has nowhere to path to. Lower bottleneck to see how the first kind connect.
    @> Saving 9 channels to channel_param_grid/channels_run001_inner_radius_1p2_sparsity_2_depth_5.pqr.
    ..
    ..
    @> The 3 sites marked sealed above report neither a channel nor a link: no route out of them survived - either narrower than bottleneck=1.60 Å, or dropped as a duplicate of a shallower site's, or the void is its own mouth and has nowhere to path to. Lower bottleneck to see how the first kind connect.
    @> Saving 5 channels to channel_param_grid/channels_run024_inner_radius_1p6_sparsity_10_depth_3.pqr.
    @> Channel calculation completed in 2.65s.
    @> Grid run 26/27: inner_radius=1.6, sparsity=10, min_depth=5
    @> The atoms supplied to calcChannels contain protein atoms only.
    @> Substituted 3766 atoms with 23638 homogeneous balls of radius 1.52 Å in 0.17s.
    @> Delaunay tessellation of 23638 points constructed in 0.74s.
    @> Surface and inner simplices filtered in 1.51s.
    @> Cavities: 88 found, 4 deeper than min_depth=5.0 Å and searched for channels, in 0.13s.
    @> Channel search (Dijkstra) over 4 search sites in 4 cavities completed in 0.06s.
    @> Found 3 channels.
    @> Search sites (sp), the void each search ran from, largest first; sp<n> tags every channel, link and output file:
    @>     site  start_point [Å]              void             volume [Å³]  depth [Å]  channels  links
    @>     sp0   [-22.642, -26.761, 3.943]    cavity 0, whole         4930       41.1         2      -
    @>     sp1   [-26.197, -29.882, -22.220]  cavity 1, whole          660        6.9         1      -
    @>     sp2   [-28.410, -21.246, -23.836]  cavity 2, whole          250        5.6         -      -  sealed
    @>     sp3   [-15.652, -41.491, -12.372]  cavity 3, whole          142        6.3         -      -  sealed
    @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
    @> The 2 sites marked sealed above report neither a channel nor a link: no route out of them survived - either narrower than bottleneck=1.60 Å, or dropped as a duplicate of a shallower site's, or the void is its own mouth and has nowhere to path to. Lower bottleneck to see how the first kind connect.
    @> Saving 3 channels to channel_param_grid/channels_run025_inner_radius_1p6_sparsity_10_depth_5.pqr.
    @> Channel calculation completed in 2.62s.
    @> Grid run 27/27: inner_radius=1.6, sparsity=10, min_depth=10
    @> The atoms supplied to calcChannels contain protein atoms only.
    @> Substituted 3766 atoms with 23638 homogeneous balls of radius 1.52 Å in 0.17s.
    @> Delaunay tessellation of 23638 points constructed in 0.74s.
    @> Surface and inner simplices filtered in 1.49s.
    @> Cavities: 88 found, 1 deeper than min_depth=10.0 Å and searched for channels, in 0.12s.
    @> Channel search (Dijkstra) over 1 search sites in 1 cavities completed in 0.05s.
    @> Found 2 channels.
    @> The void the search ran from:
    @>     start_point [Å]            void             volume [Å³]  depth [Å]  channels
    @>     [-22.642, -26.761, 3.943]  cavity 0, whole         4930       41.1         2
    @>     (site volumes measure the void itself and are not on the swept-sphere scale of the channel volumes)
    @> Saving 2 channels to channel_param_grid/channels_run026_inner_radius_1p6_sparsity_10_depth_10.pqr.
    @> Channel calculation completed in 2.58s.
    @> Number of PQR files: 27
    @> Resolution: 0.5
    @> max_proc: 2
    @> Calculating overlaps using 2 processes.
    @> 975 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1010 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 720 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 720 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 445 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 445 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 975 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1475 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 720 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 445 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 1250 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1475 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1045 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1250 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1475 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1045 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 1250 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 940 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1045 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 760 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 695 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 855 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 760 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 855 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 695 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 760 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 695 atoms and 1 coordinate sets were parsed in 0.00s.
    @> Overlap written to: channel_param_grid/channel_parameter_occupancy.pdb
    @> Number of occupied overlap voxels: 29093
    @> Channel parameters scan completed in 85.37s.


The results can be visualized using VMD_ as shown below. Red regions 
correspond to the channel's regions present in the largest number of 
results (across various parameters), and blue regions (displayed as 
dots) correspond only to certain parameters.

.. figure:: images/cavitracer_figure34.jpg
   :scale: 50 %


If other values are required, the scan can be customized as follows:

.. ipython:: python
   :verbatim:

   channels_all, parameter_sets, occupancy_file = scanChannelParameters(
       protein,
       inner_radius_values=[1.2, 1.5, 2.0],
       sparsity_values=[1.0, 3.0],
       min_depth_values=[5.0, 10.0, 15.0],
       output_path='channel_param_grid_custom')



Surface Cavities
-------------------------------------------------------------------------------

For surface cavity analysis, we will use the S. aureus Sortase A, the NMR 
structure (PDB ID: 2KID). 

We first load the protein structure and select protein structure and model 0
in the NMR ensemble. This step removes the bound substrate peptide from the 
analysis, allowing the cavity detection procedure to identify the surface 
groove that accommodates the substrate.

.. ipython:: python
   :verbatim:

   protein = parsePDB('2KID').select('protein')


.. parsed-literal::

   @> Connecting wwPDB FTP server RCSB PDB (USA).
   @> Downloading PDB files via FTP failed, trying HTTP.
   @> 2kid downloaded (2kid.pdb.gz)
   @> PDB download via HTTP completed (1 downloaded, 0 failed).
   @> 2437 atoms and 20 coordinate set(s) were parsed in 0.33s.


.. ipython:: python
   :verbatim:

   protein.setACSIndex(0)

Now, we use model 1 for the analysis using :func:`.scanSurfaceCavityParameters`.

By default, the function tests all combinations of the following parameters:

- ``surf_radius``: ``4.0``, ``4.5``, and ``5.0`` Å,
- ``inner_radius``: ``1.5`` and ``2.0`` Å,
- ``min_depth``: ``1.5`` and ``2.0`` Å,
- ``max_depth``: ``2.5`` and ``3.0`` Å,
- ``min_volume``: ``None`` and ``50`` Å³.

This gives 48 parameter combinations in total. For each combination, surface
cavities are calculated and saved as a PQR file. The function also generates an
occupancy map showing which surface regions are detected repeatedly across the
parameter scan. Regions with higher occupancy correspond to cavities that are
less sensitive to parameter choice, whereas low-occupancy regions indicate
cavities detected only under selected parameter settings.


.. ipython:: python
   :verbatim:

   cavities_all, parameter_sets, occupancy_file = scanSurfaceCavityParameters(protein)


.. parsed-literal::

    @> Calculating surface cavities for 48 parameter combinations.
    @> Grid run 1/48: surf_radius=4, inner_radius=1.5, min_depth=1.5, max_depth=2.5, min_volume=None
    @> The atoms supplied to calcChannels contain protein atoms only.
    @> The structure carries its hydrogens (101% of what a complete protein would hold), so inner_radius=1.50 is more conservative than it needs to be: the 1.2 Å floor exists only to keep a sub-water probe out of the space that missing hydrogens leave open, and 
    @> Substituted 2403 atoms with 31123 homogeneous balls of radius 1.20 Å in 0.11s.
    @> Delaunay tessellation of 31123 points constructed in 1.17s.
    @> Surface and inner simplices filtered in 0.40s.
    @> Cavities: 281 found, 30 deeper than min_depth=1.5 Å and kept, in 0.26s.
    @> Returning surface cavities
    @> Saving surface cavities to surface_cavity_parameter_grid/surface_cavities_run000_surf_radius_4_inner_radius_1p5_mindepth_1p5_maxdepth_2p5_minvol_none.pqr.
    @> Surface cavity calculation completed in 2.04s.
    @> Grid run 2/48: surf_radius=4, inner_radius=1.5, min_depth=1.5, max_depth=2.5, min_volume=50.0
    @> The atoms supplied to calcChannels contain protein atoms only.
    @> The structure carries its hydrogens (101% of what a complete protein would hold), so inner_radius=1.50 is more conservative than it needs to be: the 1.2 Å floor exists only to keep a sub-water probe out of the space that missing hydrogens leave open, and 
    @> Substituted 2403 atoms with 31123 homogeneous balls of radius 1.20 Å in 0.10s.
    @> Delaunay tessellation of 31123 points constructed in 1.19s.
    @> Surface and inner simplices filtered in 0.41s.
    @> Cavities: 281 found, 30 deeper than min_depth=1.5 Å and kept, in 0.26s.
    @> Returning surface cavities
    @> Saving surface cavities to surface_cavity_parameter_grid/surface_cavities_run001_surf_radius_4_inner_radius_1p5_mindepth_1p5_maxdepth_2p5_minvol_50.pqr.
    @> Surface cavity calculation completed in 2.06s.
    @> Grid run 3/48: surf_radius=4, inner_radius=1.5, min_depth=1.5, max_depth=3, min_volume=None
    @> The atoms supplied to calcChannels contain protein atoms only.
    @> The structure carries its hydrogens (101% of what a complete protein would hold), so inner_radius=1.50 is more conservative than it needs to be: the 1.2 Å floor exists only to keep a sub-water probe out of the space that missing hydrogens leave open, and 
    @> Substituted 2403 atoms with 31123 homogeneous balls of radius 1.20 Å in 0.10s.
    @> Delaunay tessellation of 31123 points constructed in 1.17s.
    @> Surface and inner simplices filtered in 0.40s.
    @> Cavities: 281 found, 30 deeper than min_depth=1.5 Å and kept, in 0.26s.
    @> Returning surface cavities
    @> Saving surface cavities to surface_cavity_parameter_grid/surface_cavities_run002_surf_radius_4_inner_radius_1p5_mindepth_1p5_maxdepth_3_minvol_none.pqr.
    @> Surface cavity calculation completed in 2.04s.
    @> Grid run 4/48: surf_radius=4, inner_radius=1.5, min_depth=1.5, max_depth=3, min_volume=50.0
    @> The atoms supplied to calcChannels contain protein atoms only.
    @> The structure carries its hydrogens (101% of what a complete protein would hold), so inner_radius=1.50 is more conservative than it needs to be: the 1.2 Å floor exists only to keep a sub-water probe out of the space that missing hydrogens leave open, and 
    @> Substituted 2403 atoms with 31123 homogeneous balls of radius 1.20 Å in 0.11s.
    @> Delaunay tessellation of 31123 points constructed in 1.15s.
    @> Surface and inner simplices filtered in 0.40s.
    @> Cavities: 281 found, 30 deeper than min_depth=1.5 Å and kept, in 0.27s.
    @> Returning surface cavities
    @> Saving surface cavities to surface_cavity_parameter_grid/surface_cavities_run003_surf_radius_4_inner_radius_1p5_mindepth_1p5_maxdepth_3_minvol_50.pqr.
    @> Surface cavity calculation completed in 2.03s.
    @> Grid run 5/48: surf_radius=4, inner_radius=1.5, min_depth=2, max_depth=2.5, min_volume=None
    @> The atoms supplied to calcChannels contain protein atoms only.
    @> The structure carries its hydrogens (101% of what a complete protein would hold), so inner_radius=1.50 is more conservative than it needs to be: the 1.2 Å floor exists only to keep a sub-water probe out of the space that missing hydrogens leave open, and 
    @> Substituted 2403 atoms with 31123 homogeneous balls of radius 1.20 Å in 0.11s.
    @> Delaunay tessellation of 31123 points constructed in 1.15s.
    @> Surface and inner simplices filtered in 0.35s.
    @> Cavities: 281 found, 17 deeper than min_depth=2.0 Å and kept, in 0.23s.
    @> Returning surface cavities
    @> Saving surface cavities to surface_cavity_parameter_grid/surface_cavities_run004_surf_radius_4_inner_radius_1p5_mindepth_2_maxdepth_2p5_minvol_none.pqr.
    ..
    ..
    @> Surface cavity calculation completed in 1.74s.
    @> Grid run 48/48: surf_radius=5, inner_radius=2, min_depth=2, max_depth=3, min_volume=50.0
    @> The atoms supplied to calcChannels contain protein atoms only.
    @> The structure carries its hydrogens (101% of what a complete protein would hold), so inner_radius=2.00 is more conservative than it needs to be: the 1.2 Å floor exists only to keep a sub-water probe out of the space that missing hydrogens leave open, and
    @> Substituted 2403 atoms with 31123 homogeneous balls of radius 1.20 Å in 0.09s.
    @> Delaunay tessellation of 31123 points constructed in 1.07s.
    @> Surface and inner simplices filtered in 0.33s.
    @> Cavities: 244 found, 19 deeper than min_depth=2.0 Å and kept, in 0.18s.
    @> Returning surface cavities
    @> Saving surface cavities to surface_cavity_parameter_grid/surface_cavities_run047_surf_radius_5_inner_radius_2_mindepth_2_maxdepth_3_minvol_50.pqr.
    @> Surface cavity calculation completed in 1.73s.
    @> Number of PQR files: 48
    @> Resolution: 0.5
    @> max_proc: 2
    @> Calculating overlaps using 2 processes.
    @> 1961 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1794 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1734 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1730 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 2076 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 654 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 591 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 685 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 622 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 1847 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1679 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1617 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 512 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 512 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 543 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 543 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 2214 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1979 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 2359 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 2120 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 850 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 1979 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 728 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 902 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 779 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 1897 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 688 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 673 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 2124 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 2038 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 740 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 724 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 2472 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 2236 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 2629 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 2387 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 2263 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 2154 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 2420 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 2305 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1087 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 932 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1164 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1005 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 1015 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 911 atoms and 1 coordinate sets were parsed in 0.00s.
    @> 1092 atoms and 1 coordinate sets were parsed in 0.01s.
    @> 984 atoms and 1 coordinate sets were parsed in 0.01s.
    @> Overlap written to: surface_cavity_parameter_grid/surface_cavity_param_occupancy.pdb
    @> Number of occupied overlap voxels: 33322
    @> Surface cavity parameters scan completed in 97.08s.


The results can be visualized using VMD_. As we can see, the red
region (present in at least 50% of results; ``occupancy > 0.5``)
is also shown at the ligand binding spot.

.. figure:: images/cavitracer_figure35.jpg
   :scale: 50 %

