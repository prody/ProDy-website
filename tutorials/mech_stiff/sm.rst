.. _mech_stiff:

Mechanical Stiffness Calculations
===============================================================================

This example shows how to perform mechanical resistance calculation for GFP
protein (**1gfl**) and visualize the results using Matplotlib_ library and VMD_
program.

An :class:`.ANM` instance that stores Hessian matrix and normal mode data 
describing the intrinsic dynamics of the protein structure will be used as 
an input (*model*) as well as cooridinates of protein structure (*coords*, *pdb*).

See [EB08]_ for more information about the theory of mechanical resistance 
calculations and more examples.

.. [EB08] Eyal E., Bahar I. Toward a Molecular Understanding of 
   the Anisotropic Response of Proteins to External Forces: Insights from 
   Elastic Network Models. *Biophys. J.* **2008** 94:3424-34355.


Parse structure
-------------------------------------------------------------------------------

We start by importing everything from the ProDy package:

.. ipython:: python

   from prody import *
   from pylab import *
   ion()   # turn interactive mode on

We start by parsing chain A of PDB structure **1gfl**.

.. ipython:: python

   gfp, header = parsePDB('1gfl', header=True, chain='A')
   gfp

We want to use only Cα atoms for our calculations, so we select them into 
a new object *calphas*:

.. ipython:: python

   calphas = gfp.ca
   calphas


Build Hessian
-------------------------------------------------------------------------------

In the next step we instantiate an :class:`.ANM` instance:

.. ipython:: python

   anm = ANM('GFP ANM analysis')

Then, build the Hessian matrix by passing selected atoms (230 Cα's) to
:meth:`.ANM.buildHessian` method:

.. ipython:: python

   anm.buildHessian(calphas, cutoff=13.0)


Those two actions are required and sufficient to perform mechanical stiffness 
calculations without calculating modes. The Hessian defines the elastic network model, 
which will be used to calculate mechanical striffness.


Stiffness Matrix Calculations
-------------------------------------------------------------------------------

Mechanical stiffness calculations for the selected group of atoms can be 
performed using :meth:`.ANM.buildMechStiff` method:

.. ipython:: python

   anm.buildMechStiff(calphas)
   anm.getStiffness()

Mechanical stiffness matrix is available using the :meth:`.ANM.getStiffness` 
method. To show the stiffness matrix as an image map use the following function:

.. ipython:: python
	
   showMechStiff(anm, calphas, 'jet_r')


Note that 'jet_r' will reverse the colormap of image map which will be 
similar to coloring method of VMD_ program. 

The mean values of the mechanical stiffness matrix for each residue 
can be calculated using :meth:`showMeanMechStiff` function where 
the secoundary structure of the protein is drawn using header information.

.. ipython:: python

   showMeanMechStiff(anm, calphas, header, 'A', 'jet_r')

 
Mechanical Stiffness in VMD
-------------------------------------------------------------------------------

We can generate tcl files for visualizing mechanical stiffness with VMD_ 
using the :func:`.writeVMDstiffness` function. Select one residue in *indices* (**[3]**) 
or series of residues (**[3, 7]**, means from 3 aa to 7 aa including) and 
a range of effective spring constant *k_range* (**[0, 7.5]**). 

We provide *gfp* as well as *calphas* so VMD_ has information about the complete protein structure,
which it can use for graphical representations.


.. ipython:: python
   :verbatim:

   writeVMDstiffness(anm, gfp, [3,7], [0,7.5], filename='1gfl_3-7aa', loadToVMD=False)
   writeVMDstiffness(anm, gfp, [3], [0,7], filename='1gfl_3')

Results will be loaded automatically to VMD_ by default. Use ``loadToVMD=False`` to 
change it. The TCL file will be saved automatically and can be used later by using 
linux command line: 

::  vmd -e 1gfl_3aa.tcl

or in VMD_ *TKConsole* (*VMD Main*) for Linux, Windows and Mac users: 
::  play 1gfl_3aa.tcl


The tcl file contains a method for drawing lines between selected pairs of 
residues, which are highlighted as spheres. The color of the line can be modified 
by changing the ``draw color red`` line in the output file. Only colors from VMD_ 
Coloring Method will work. Other changes can be done within VMD_ in the
*Graphical Representations* menu.

.. figure:: images/1gfl_chA.png
   :scale: 60 %

The figure shows GFP results from :meth:`.vmdfile.writeVMDstiffness` method opened in VMD_. 
Pairs of found residues LYS3-GLY116, LYS3-PRO211 and PRO211-ASN212 are shown as VDW 
spheres connected with red lines.

Additionally, :file:`1gfl_3aa.txt` file is created. It contains a list 
of residue pairs with the value of effective spring constant (in a.u. because 
*kBT=1*) obtained from :meth:`.ANM.buildMechStiff` method.
::

     LYS3    GLY116  6.91650667766
     LYS3    PRO211  6.85989128668
     LYS3    ASN212  6.69507284967
     ...


The range of spring constant for *k_range* can be checked as follows:  

.. ipython:: python

   anm.getStiffnessRange()

See also :meth:`.ANM.getMechStiffStatistic` and :meth:`.ANM.getStiffnessRangeSel`
functions for more detailed analysis of the stiffness matrix.

The results of the mean value of mechanical stiffness calculation can be seen 
in VMD_ using:

.. ipython:: python
   :verbatim:
	
   writeDeformProfile(anm, gfp, selstr='chain A and name CA', pdb_selstr='protein')


.. figure:: images/1gfl_defprofile_vmd.png
   :scale: 90 %



Calculate Distribution of Deformation 
-------------------------------------------------------------------------------

The distribution of the deformation in the distance contributed by each mode 
for a selected pair of residues has been described in [EB08]_, see *Eq. (10)*
and plots are shown on *Fig. (2)*. 

These results can be plotted using :meth:`.plotting.showPairDeformationDist` 
or a list can be obtained using :meth:`.analysis.calcPairDeformationDist`.

.. ipython:: python

   calcPairDeformationDist(anm, calphas, 3, 132)

   showPairDeformationDist(anm, calphas, 3, 132)


Figure shows the plotted distribution for deformations between 3-132 residue in each mode *k*.

To obtain results without saving any file type:

.. ipython:: python

   d1 = calcPairDeformationDist(anm, calphas, 3, 212)
   d2 = calcPairDeformationDist(anm, calphas, 132, 212)
   print d1[0], d1[1]

   plot(d1[0], d1[1], 'k-', d2[0], d2[1], 'r-')

