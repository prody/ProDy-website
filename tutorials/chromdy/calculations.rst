Calculations
===============================================================================

Here are the required imports again. You do not need to repeat them if you are
still in the same python session.

.. ipython:: python

    from prody import *
    from pylab import *
    ion()


Loading a HiC file and preparing the data
-------------------------------------------------------------------------------
First, we parse structural data from HiC, which is essentially a contact matrix 
for constructing an elastic network model. The output is a :class:`.HiC` class instance.

.. ipython:: python
   :verbatim:

    hic = parseHiC(hic_file)

Then, we can normalize the data as follows::

.. ipython:: python
   :verbatim:

    hic.normalize()

By default ProDy uses the VCnorm, but you can choose different normalization schemes too.
If the normalization scheme that you desired is not supported yet, you can implement 
your own normalization function and supply it to hic.normalize.

Since Hi-C data contains regions that don't have data (e.g. telomere, centromere), 
these regions are "masked" by the "HiC" class. You can access the mask as follows::

.. ipython:: python
   :verbatim:

    print(hic.mask)

You can also turn the mask on and off::

.. ipython:: python
   :verbatim:

    hic.masked = True
    hic.masked = False

Gaussian Network Model (GNM) Analysis
-------------------------------------------------------------------------------

To apply the GNM to the HiC, it is a best practice to turn the mask on::

.. ipython:: python
   :verbatim:

    hic.masked = True
    gnm = hic.calcGNM()

The output object `gnm` will be a :class:`.MaskedGNM` instance. It behaves similarly to 
a normal :class:`.GNM` class, and the difference is that it comes with a `gnm.mask` property 
just like the HiC class. 

It also has a special :meth:`fixTail` function that is implemented specifically for Hi-C 
in case the tail is missing. For example, if your gnm has only 890 loci, whereas your 
ATAC-seq data has 900 loci (and you know it is the tail that is missing, not the head), 
then you can do the following::

.. ipython:: python
   :verbatim:

    gnm.fixTail(900)

There are other useful functions too, such as :function:`.chromatin.cluster` for domain 
identification. All these functions are in chromatin module, so feel free to check them out. 
