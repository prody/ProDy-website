Background
==========

Biological functions of biomolecules is guided by conformational transitions, which involve 
conformational states that are short-lived, and therefore difficult (or even impossible) to detect 
experimentally.  Here, computational tools are useful to predict plausible transition pathways, which can 
be probed experimentally.  We developed an efficient method to construct an energetically favorable 
pathway between two experimental structures, which represent stable states.

Methodology
===========

In CoMD, the protein is represented as coarse-grained elastic network model (ENM), accounting for the 
topology of inter-residue contacts.  The potential is two-state potential obtained by combining two 
ENMs, denoted here as end states.  This potential has a cusp hypersurface where energies from both the 
states are the same.  A minimum energy structure on the cusp hypersurface is identified, and is treated 
as a transition state.  Conformers collected from these two steepest descent paths, along with the 
transition state provide the pathway.   This pathway can be regarded as the minimum energy path 
between the end structures.  A detailed description of the algorithm can be found in [MG13]_.  In short:

1. Two end structures XA and XB are represented by the positions of their C-alpha, which are then aligned 
(structural superimposition).  XT  is a conformer for which the energy is equal from both XA and XB.   

2. Starting from XT (n), a transition state is identified by an iterative procedure.  At each step, a steepest 
descent minimization is performed for each of the two end-states XA and XB,  generating two new sets of 
coordinates XA(n+1) and XB(n+1).  

3.  A linear interpolation is performed between XA(n+1) and XB(n+1) to determine the conformer that 
residues on the cusp hpersurface.  This is the new transition state XT(n+1). The steps 1 and 2 are 
iterated until the energy difference between the two transition state XT(n) and XT(n+1) conformers 
are within a tolerance.

4.  Two separate steepest descent minimizations are then performed, starting from the final transition 
state XTf and conformers separated by a user-defined RMSD are collected.

5. These conformers are then indexed to construct a pathway, known here as the transition pathway.
