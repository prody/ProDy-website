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


Next, :func:`scanChannelParameters` can be use to explore various parameters 
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
   @> Grid run 1/27: inner_radius=1.2, sparsity=1, min_depth=3
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 3766 atoms with 23638 homogeneous balls of radius 1.52 A in 0.27s.
   @> Delaunay tessellation of 23638 points constructed in 0.88s.
   @> Surface and inner simplices filtered in 1.18s.
   @> 17 surface cavities detected and filtered in 0.26s.
   @> Channel pathfinding (graph Dijkstra) over 17 cavities completed in 0.46s.
   @> Detected 26 channels.
   @> Saving results to channel_param_grid/channels_run000_inner_radius_1p2_sparsity_1_depth_3.pqr.
   @> Channel calculation completed in 3.06s.
   @> Grid run 2/27: inner_radius=1.2, sparsity=1, min_depth=5
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 3766 atoms with 23638 homogeneous balls of radius 1.52 A in 0.20s.
   @> Delaunay tessellation of 23638 points constructed in 0.83s.
   @> Surface and inner simplices filtered in 1.16s.
   @> 7 surface cavities detected and filtered in 0.20s.
   @> Channel pathfinding (graph Dijkstra) over 7 cavities completed in 0.31s.
   @> Detected 13 channels.
   @> Saving results to channel_param_grid/channels_run001_inner_radius_1p2_sparsity_1_depth_5.pqr.
   @> Channel calculation completed in 2.70s.
   @> Grid run 3/27: inner_radius=1.2, sparsity=1, min_depth=10
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 3766 atoms with 23638 homogeneous balls of radius 1.52 A in 0.18s.
   @> Delaunay tessellation of 23638 points constructed in 0.83s.
   @> Surface and inner simplices filtered in 1.15s.
   @> 3 surface cavities detected and filtered in 0.20s.
   @> Channel pathfinding (graph Dijkstra) over 3 cavities completed in 0.27s.
   @> Detected 7 channels.
   ..
   ..
   @> Channel calculation completed in 2.23s.
   @> Grid run 27/27: inner_radius=1.6, sparsity=5, min_depth=10
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 3766 atoms with 23638 homogeneous balls of radius 1.52 A in 0.18s.
   @> Delaunay tessellation of 23638 points constructed in 0.80s.
   @> Surface and inner simplices filtered in 1.03s.
   @> 1 surface cavities detected and filtered in 0.10s.
   @> Channel pathfinding (graph Dijkstra) over 1 cavities completed in 0.07s.
   @> Detected 2 channels.
   @> Saving results to channel_param_grid/channels_run026_inner_radius_1p6_sparsity_5_depth_10.pqr.
   @> Channel calculation completed in 2.19s.
   @> Number of PQR files: 27
   @> Resolution: 0.5
   @> max_proc: 2
   @> Calculating overlaps using 2 processes.
   @> 2345 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2985 atoms and 1 coordinate sets were parsed in 0.03s.
   @> 1975 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2495 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1650 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1990 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1955 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 2025 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1585 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1260 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1625 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1420 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1910 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1620 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1525 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1415 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1235 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1030 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1455 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1385 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1235 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1230 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1060 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1055 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1010 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 855 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 680 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Overlap written to: channel_param_grid/channel_parameter_occupancy.pdb
   @> Number of occupied overlap voxels: 30028
   @> Channel parameters scan completed in 80.07s.


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

Now, we use model 1 for the analysis using :func:`scanSurfaceCavityParameters`.

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
   @> Substituted 2403 atoms with 31123 homogeneous balls of radius 1.20 A in 0.16s.
   @> Delaunay tessellation of 31123 points constructed in 1.30s.
   @> Surface and inner simplices filtered in 0.39s.
   @> 30 surface cavities detected and filtered in 0.31s.
   @> Returning surface cavities
   @> Saving surface cavities to surface_cavity_parameter_grid/surface_cavities_run000_surf_radius_4_inner_radius_1p5_mindepth_1p5_maxdepth_2p5_minvol_none.pqr.
   @> Surface cavity calculation completed in 2.27s.
   @> Grid run 2/48: surf_radius=4, inner_radius=1.5, min_depth=1.5, max_depth=2.5, min_volume=50.0
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 2403 atoms with 31123 homogeneous balls of radius 1.20 A in 0.11s.
   @> Delaunay tessellation of 31123 points constructed in 1.20s.
   @> Surface and inner simplices filtered in 0.38s.
   @> 30 surface cavities detected and filtered in 0.30s.
   @> Returning surface cavities
   @> Saving surface cavities to surface_cavity_parameter_grid/surface_cavities_run001_surf_radius_4_inner_radius_1p5_mindepth_1p5_maxdepth_2p5_minvol_50.pqr.
   @> Surface cavity calculation completed in 2.07s.
   @> Grid run 3/48: surf_radius=4, inner_radius=1.5, min_depth=1.5, max_depth=3, min_volume=None
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 2403 atoms with 31123 homogeneous balls of radius 1.20 A in 0.10s.
   @> Delaunay tessellation of 31123 points constructed in 1.19s.
   @> Surface and inner simplices filtered in 0.37s.
   @> 30 surface cavities detected and filtered in 0.30s.
   @> Returning surface cavities
   @> Saving surface cavities to surface_cavity_parameter_grid/surface_cavities_run002_surf_radius_4_inner_radius_1p5_mindepth_1p5_maxdepth_3_minvol_none.pqr.
   @> Surface cavity calculation completed in 2.07s.
   @> Grid run 4/48: surf_radius=4, inner_radius=1.5, min_depth=1.5, max_depth=3, min_volume=50.0
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 2403 atoms with 31123 homogeneous balls of radius 1.20 A in 0.11s.
   @> Delaunay tessellation of 31123 points constructed in 1.21s.
   @> Surface and inner simplices filtered in 0.38s.
   @> 30 surface cavities detected and filtered in 0.30s.
   @> Returning surface cavities
   @> Saving surface cavities to surface_cavity_parameter_grid/surface_cavities_run003_surf_radius_4_inner_radius_1p5_mindepth_1p5_maxdepth_3_minvol_50.pqr.
   @> Surface cavity calculation completed in 2.10s.
   @> Grid run 5/48: surf_radius=4, inner_radius=1.5, min_depth=2, max_depth=2.5, min_volume=None
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 2403 atoms with 31123 homogeneous balls of radius 1.20 A in 0.10s.
   @> Delaunay tessellation of 31123 points constructed in 1.20s.
   @> Surface and inner simplices filtered in 0.38s.
   @> 17 surface cavities detected and filtered in 0.30s.
   ..
   ..
   @> Saving surface cavities to surface_cavity_parameter_grid/surface_cavities_run046_surf_radius_5_inner_radius_2_mindepth_2_maxdepth_3_minvol_none.pqr.
   @> Surface cavity calculation completed in 1.85s.
   @> Grid run 48/48: surf_radius=5, inner_radius=2, min_depth=2, max_depth=3, min_volume=50.0
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 2403 atoms with 31123 homogeneous balls of radius 1.20 A in 0.10s.
   @> Delaunay tessellation of 31123 points constructed in 1.11s.
   @> Surface and inner simplices filtered in 0.35s.
   @> 19 surface cavities detected and filtered in 0.23s.
   @> Returning surface cavities
   @> Saving surface cavities to surface_cavity_parameter_grid/surface_cavities_run047_surf_radius_5_inner_radius_2_mindepth_2_maxdepth_3_minvol_50.pqr.
   @> Surface cavity calculation completed in 1.84s.
   @> Number of PQR files: 48
   @> Resolution: 0.5
   @> max_proc: 2
   @> Calculating overlaps using 2 processes.
   @> 1794 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1961 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1734 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1730 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 654 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 2076 atoms and 1 coordinate sets were parsed in 0.01s.
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
   @> 902 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1897 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 779 atoms and 1 coordinate sets were parsed in 0.00s.
   ..
   ..
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
   @> 911 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 1092 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 984 atoms and 1 coordinate sets were parsed in 0.01s.
   @> Overlap written to: surface_cavity_parameter_grid/surface_cavity_param_occupancy.pdb
   @> Number of occupied overlap voxels: 33322
   @> Surface cavity parameters scan completed in 100.97s.

The results can be visualized using VMD_. As we can see, the red
region (present in at least 50% of results; ``occupancy > 0.5``)
is also shown at the ligand binding spot.

.. figure:: images/cavitracer_figure35.jpg
   :scale: 50 %

