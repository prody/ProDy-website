.. _essa_profile:

Generation of ESSA profile 
===============================================================================

First, we will make necessary imports from ProDy_, NumPy_ and Matplotlib_
packages if you haven't already done it:

.. ipython:: python

    from prody import *
    from numpy import *
    from matplotlib.pyplot import *
    ion()


First, we parse a structure on which we want to perform ESSA. 
For this tutorial, we will fetch the TEM-1 beta-lactamase (1pzo) from the PDB.

.. ipython:: python

    ag = parsePDB('1pzo', compressed=False)

ESSA is implemented as a ProDy class, so we can instantiate an :class:`.ESSA` 
object.

.. ipython:: python

    essa = ESSA()

Before starting residue scanning, we first need to set the system. Even though 
the ligand(s) are not included in the scanning, they can be specified for further 
analysis. For this purpose, users can simply give the chain ids and residue numbers 
as a string shown in the following. The protein residues interacting with each ligand 
can be determined within a specific cutoff distance (default ``dist=4.0``). 
ESSA will use the title of the atomic object provided to it in the names of the 
output files.

ESSA, by default, generates a :class:`.ModeEnsemble` containing the modes 
resulting from the perturbation of each residue. However, if the user does not 
have enough memory resources for conducting computation on a large structure, 
then ``lowmem=True`` should be chosen, thus only storing eigenvalues/eigenvectors 
for each perturbed model. 

In the structure 1pzo, there are two identical allosteric inhibitors (residues 300 and 301 of chain A).
ESSA can be informed about them when setting the system as follows:

.. ipython:: python

    essa.setSystem(ag, lig='A 300 A 301')

Scanning can be done using :func:`.scanResidue`. Here one can specify the number 
of softest modes, the type of ENM (GNM or ANM), and the cutoff, which by default 
are ``n_modes=10``, ``enm='gnm'``, ``cutoff=None``, respectively. If ``cutoff`` 
is not specified, its value is adopted from the default value of the specified ENM. 
During the scanning, the progress will also be displayed.

.. ipython:: python

    essa.scanResidues()

The generated ESSA profile can be shown by :func:`.showESSAProfile`. On this profile, 
the residues interacting with the previously specified ligand(s) are automatically 
highlighted. The blue dashed baseline shows the q-th quantile of the profile, which is by 
default ``q=0.75``, representing the top quartile. Other residues of interest 
can also be indicated with their single-letter code and residue number on this plot 
if specified by a ProDy selection string provided as a parameter. Users can dynamically 
customize the properties of the plot using the matplotlib context manager as shown below.

.. ipython:: python

    with style.context({'figure.figsize': (10, 5), 'figure.dpi': 300}):
        essa.showESSAProfile()

ESSA z-scores can be obtained as a NumPy array using :func:`.getESSAZscores`, and saved with 
:func:`.saveESSAZscores`.

.. ipython:: python

    essa.getESSAZscores()[:10]

.. ipython:: python

    essa.saveESSAZscores()

In order to visualize the essential residues, a PDB file can be generated, in 
which the z-scores are written in the B-factor column. Later, this file can be 
opened in a molecular graphics program such as PyMOL or VMD, where the structure 
can be colored according to the B-factors. 

.. ipython:: python

    essa.writeESSAZscoresToPDB()

Please check the other getter and save methods and their docstrings, such as those 
for ligand binding residues.