.. _setup:

Simulation Setup
================

System that contains protein structures with water and counter ions in collective 
molecular dynamics simulations can be prepared using the following
interface:

.. figure:: ../../_static/figures/comd_setup_gui.png
   :scale: 80%

Input Files
-----------

:file:`.pdb` files for both initial and final structure of the protein are required from the user. You can learn how to prepare these files from `NAMD tutorials`_.


Protein Structures
------------------

1. First, select initial :file:`.pdb` (structure) files and final :file:`.pdb` 
   (structure) files. You can use the adelynate kinase files provided for 
   this tutorial. (1ake.pdb, 4ake.pdb) Alternatively, structure and 
   coordinate files for a protein of interest can also be used. These 
   files should contain all atoms required for protein stability and
   function, these may include cofactors and metal atoms. Crystallographic 
   water molecules may also be retained in the structure. 

2. Select the corresponding chain ID's for initial and final structure. 
   For adelynate kinase both chain ID's should be 'A'. 

Ionization Parameters
--------------------- 

1. To build the topology of the structure, a topology file needs to be
   provided. This file will parametrize the topology of structures e.g.
   bond lengths and bond angles. If not given, it will use CHARMM22 
   all-hydrogen topology file for proteins and lipids. For more information
   visit `Topology`_ pages. 

2. Enter simulation box padding (distance from the protein to box surfaces)
   and select whether you would like to add counter ions.
   Collective Molecular Dynamics GUI uses `Solvate`_ and `Autoionize`_ plugins to add
   water and ions.

Minimization Parameters
-----------------------

1. For every cycle, the structures (protein and water) will be minimized by using NAMD. 
   For parametrization of molecular dynamics simulation in NAMD, an input file is required. 
   If not given, it will use Combined CHARMM22 all-hydrogen parameter file for CHARMM22 
   proteins and CHARMM27 lipids. For more information visit `Parameter`_ pages. 

2. For minimization, two additional parameters are required. The first one is the temperature
   that the minimization will be performed and the second one is the total number of 
   steps to perform minimization. Default parameters are 310 K and 50 steps. 

Conformational Sampling with ANM Modes
--------------------------------------

For sampling of conformations with ANM modes, two parameters are required. ANM cutoff
is the parameter to determine the contacts on the proteins for ANM calculations. 
The second parameter is the deviation parameter. For each step in conformational sampling
there is a scale factor for perturbing structure on the direction of ANM mode. 

Targeted Molecular Dynamics Parameters
--------------------------------------

For every cycle, after conformational sampling with ANM modes toward to final structure, 
the side chains of the protein needs to move by performing a fast targeted molecular dynamics
simulation (TMD). In TMD simulations, two parameters are required. The first one is the spring 
constant (kcal/mol/A^2) and the second one is the length of TMD simulation (ps). For further 
information visit `TMD`_ pages. 

Simulation Options
------------------

coMD simulations will run in cycles and the number of cycles is the global parameter for simulations. 
The simulation will run until the number of cycles achieved or the RMSD between initial and final
structure is less than 1.5 A or the reduction in the RMSD is less than 0.15 A. The other parameter
in this block is the number of physical cores on your machine for parallel NAMD simulations.  

Output Options
--------------

Select output folder, file prefix to distinguish different runs. Performing multiple simulations 
to see whether results are reproducible is always a good idea.

Output Files
------------

There will be an output folder with a python script (py), a tcl script (tcl), two :file:`pdb` and 
:file:`psf` files and two parameter files will be prepared. The tcl script should run from vmd. 

For a summary of contents of the final system, see :file:`prefix.log` file.

Simulation
----------

Now you need to run coMD simulations. To perform those simulations it is possible 
to use the following command::

   vmd –dispdev text –e your_output_prefix.tcl
￼￼

If you want to run this simulation on cluster, copy this directory on cluster and put following command on file queue system.

When simulations are complete, you can continue with following analysis steps.



.. _NAMD tutorials: http://www.ks.uiuc.edu/Training/Tutorials/
.. _Autoionize: http://www.ks.uiuc.edu/Research/vmd/plugins/autoionize/
.. _Solvate: http://www.ks.uiuc.edu/Research/vmd/plugins/solvate/
.. _Topology: http://www.ks.uiuc.edu/Research/vmd/plugins/membrane/top_all27_prot_lipid.inp
.. _Parameter: http://www.ks.uiuc.edu/Training/Tutorials/science/nanotubes/files/par_all27_prot_lipid.prm
.. _TMD: http://www.ks.uiuc.edu/Research/namd/2.10b1/ug/node47.html
