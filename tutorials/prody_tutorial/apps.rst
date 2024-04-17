.. _commands-tutorial:

Applications Tutorial
===============================================================================

You can use :ref:`prody-apps` to perform some automated tasks, such as
ANM/GNM/PCA calculations, fetching and aligning PDB files, making atom
selections, or identify contacts.  ProDy applications are handled by a script
that comes with all installation packages.  You can run the script
from a central location such as :file:`/usr/local/bin`::

  $ prody -h

or from the current working directory::

  $ ./prody -h

or::

  $ prody -h

or on Windows::

  $ C:\Python27\Scripts\prody -h

These lines will print available ProDy applications.  You can get more help
on a specific commands as follows::

  $ prody anm -h

Align PDB files
-------------------------------------------------------------------------------

:ref:`prody-align` command can be used to download and align structures for
given PDB identifiers::

  $ prody align 5uoj 1r39 1zz2

Structures will be automatically downloaded from wwPDB FTP servers and saved
in the current working directory.  Additionally, you can configure ProDy
to use a local mirror of PDB or to store downloaded files in a local folder.
See :ref:`prody-basics` part of the tutorial.


ANM calculations
-------------------------------------------------------------------------------

:ref:`prody-anm` can be used to perform ANM calculations::

  $ prody anm 5uoj -a -A

``-a`` and ``-A`` options will make ProDy output all data and figure files.


PCA calculations
-------------------------------------------------------------------------------

:ref:`prody-pca` can be used to perform PCA calculations.  The following
example will perform PCA calculations for Cα atoms of the p38 MAP kinase
using files:

  * `ProDy Tutorial files (ZIP) <prody_tutorial_files.zip>`_
  * `ProDy Tutorial files (TGZ) <prody_tutorial_files.tgz>`_

::

  $ tar -xzf p38_trajectory.tar.gz
  $ prody pca -a -A --select calpha --pdb p38.pdb p38_100frames.dcd
