.. _essa_profile:

Generation of ESSA profile 
===============================================================================

First, we parse a structure on which we want to perform ESSA. 
For this tutorial, we will fetch the TEM-1 beta-lactamase (1pzo) from the PDB.

.. ipython:: python

    pdb = parsePDB('1pzo', compressed=False)

ESSA is implemented as a ProDy :class:`.ESSA`, so we can instantiate an ESSA 
object.

.. ipython:: python

    essa = ESSA()

Before starting residue scanning, we first need to set the system. Even though 
the ligand(s) are not included in the scanning, they can be specified for further 
analysis. For this purpose, users can simply give the chain ids and residue numbers 
as a string shown below, and the protein residues interacting with each ligand 
can be determined within a specific cutoff distance (default ``dist=4.0``). 
ESSA will use the title of the parsed PDB file in the names of the output files.

ESSA, by default, generates a :class:`.ModeEnsemble` containing the modes 
resulting from the perturbation of each residue. However, if the user does not 
have enough memory resources for conducting computation on a large structure, 
then ``lowmem=True`` should be chosen, thus only storing eigenvalues/eigenvectors 
for each perturbed model. 

In the structure 1pzo, there are two identical allosteric inhibitors (A300 and A301).

.. ipython:: python

    essa.setSystem(pdb, lig='A 300 A 301')

Scanning can be done using :func:`.scanResidue`. Here one can specify the number 
of softest modes, the type of ENM (GNM or ANM), and the cutoff, which by default 
are ``n_modes=10``, ``enm='gnm'``, ``cutoff=None``, respectively. If ``cutoff`` 
is not specified, its value is adopted from the specified ENM. During the scanning, 
the progression will also be displayed.

.. ipython:: python

    essa.scanResidues()

The generated ESSA profile can be shown by :func:.showESSAProfile. On this profile, 
the residues interacting with the previously specified ligand(s) are automatically 
shown. The blue dashed baseline shows the q-th quantile of the profile, which is by 
default ``q=0.75``, representing the top quartile thereof. Other residues of interest 
can also be indicated with their single-letter code and residue number on this plot 
if specified by a ProDy selection algebra as a parameter. User can dynamically 
customize the properties of the plot by the matplotlib context manager as shown below.

.. ipython:: python

    with plt.style.context({'figure.figsize': (10, 5), 'figure.dpi': 300}):
        essa.showESSAProfile()

ESSA z-scores can be obtained by :func:`.getESSAZscores`, and saved by 
:func:`.saveESSAZscores` as a Numpy array.

.. ipython:: python

    essa.getESSAZscores()[:10]

.. ipython:: python

    essa.saveESSAZscores()

In order to visualize the essential residues, a pdb file can be generated, in 
which the z-scores are written in the B-factor column. Later, this file can be 
opened in PyMol or VMD, and the structure can be colored according to the B-factors. 

.. ipython:: python

    essa.writeESSAZscoresToPDB()

Please check the other getter/save-methods and their docstrings, such as those 
for ligand binding residues.