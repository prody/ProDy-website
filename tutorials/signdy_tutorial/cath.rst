.. _signdy-cath:

Data Collection with CATH
===============================================================================

The first step in signature dynamics analysis is to collect a set of related 
protein structures and build a :class:`.PDBEnsemble`. This can be achieved by 
multiple routes: a query search of the PDB using :func:`.blastPDB` or :func:`.searchDali`, 
extraction of PDB IDs from the Pfam or CATH database, or input of a pre-defined list. 

Here, we demonstrate the usage of CATH for ensemble building.

First, make necessary imports from ProDy_ and Matplotlib_ packages if you haven't already.

.. ipython:: python

   from prody import *
   from pylab import *
   ion()

First, we initialise a :class:`.CATHDB` object. By default, this is done 
by downloading data from the CATH website.

.. ipython:: python

   cath = CATHDB()
   cath

We can also use this object to save this data to an :file:`.xml` file and 
load it later:

.. ipython:: python

   cath.save('cath.xml')

.. ipython:: python

   cath = CATHDB('cath.xml')

One way of using the :class:`.CATHDB` class is to navigate the CATH tree, using 
modified versions of methods and properties inherited from base classes in 
:module:`~xml.etree.ElementTree`. 

The root of the tree and all other elements in it are instances of the 
:class:`.CATHElement` class, which is based on :class:`~xml.etree.ElementTree.Element`, 
allowing us to easily navigate the CATH tree structure using parent/child relationships 
as follows:

.. ipython:: python

   root = cath.root
   root

.. ipython:: python

   node = root.getchildren()
   node

Any branching point node containing a collection of children is an instance of the 
:class:`.CATHCollection` class, which is based on the :class:`.CATHElement` class 
but has additional and modified properties and methods. 

For example, collections return a list of values for the properties *cath* (CATH ID) and 
*name*, while elements return single values:

.. ipython:: python

   node.name

.. ipython:: python

   node.cath

.. ipython:: python

   element = node[0]
   element.name

.. ipython:: python

   element.cath

We can also use the :class:`.CATHDB` class to find a particular part of the CATH hierarchy 
by CATH ID:

.. ipython:: python

   node = cath.find('1.10.8')
   node.name

We can also then examine its children:

.. ipython:: python

   node.getchildren().name

We can also use it to get PDB IDs associated with particular levels:

.. ipython:: python

   node = cath.find('1.10.8.40')
   node.getPDBs()

Another useful method is for seeing the associated CATH domains 
and the associated selection strings.

.. ipython:: python

   node.getDomains()

.. ipython:: python

   node.getSelStrs()

We can combine all of these together to fetch and parse structures from 
the PDB and make the appropriate selections at the same time:

.. ipython:: python

   proteins = node.parsePDBs(subset='ca')
   proteins

This then allows us to build a :class:`.PDBEnsemble` from them:

.. ipython:: python

   ens = buildPDBEnsemble(proteins, mapping='CE')
   ens

Lastly, the :class:`.CATHDB` object can be used to find different CATH domains within 
a particular PDB structure:

.. ipython:: python

   result = cath.search('3kg2A')
   result.name

.. ipython:: python

   result.getSelstrs()

This iGluR example also illustrates that CATH domains may also not correspond to biological domains 
identified by other methods as the N-terminal domain (NTD; residues 1 to 376), a type-I PBP domain, 
is split into CATH domains corresponding to the two lobes, which each belong to 'Superfamily 3.40.50.2300'. 

Likewise, the two lobes of the ligand-binding domain (LBD) are assigned as separate domains that both belong 
to 'Periplasmic binding protein-like II', which is usually the whole bi-lobed clamshell structure.
