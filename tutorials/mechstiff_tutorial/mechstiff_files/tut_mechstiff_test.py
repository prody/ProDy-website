#!/usr/bin/env python

from prody import *

'''Example for the GFP as in:
Eyal E., Bahar I. Toward a Molecular Understanding of 
the Anisotropic Response of Proteins to External Forces:
Insights from Elastic Network Models. *Biophys J* **2008** 94:3424-34355.

Problems: karolami@pitt.edu
'''

gfp 	      = parsePDB('1gfl')
calphas       = gfp.select('protein and chain A and name CA')

anm = ANM('prot analysis')
anm.buildHessian(calphas, cutoff=13)

#anm.py: buildSM(self, coords, n_modes=None, kbt=1., saveMap=False, saveMatrix=False, filename='sm')
anm.buildSM(calphas, kbt=1, saveMap=True, saveMatrix=True, filename='1gfl_stiffmatrix')

# analysis.py: calcPairDeformationDist(model, coords, ind1, ind2, kbt=1., saveFile=False, filename='out', savePlot=False)
calcPairDeformationDist(anm, calphas, 132, 212, kbt=1, saveFile=True, savePlot=True, filename='1gfl_132_212')
calcPairDeformationDist(anm, calphas, 3, 132, kbt=1, saveFile=True, savePlot=True, filename='1gfl_3_132')

# vmdfile.py: writeVMDstiffness(model, pdb, indices, k_range, filename='vmd_out', selstr='protein and name CA', loadToVMD=True)
pdb = gfp.select('chain A')
writeVMDstiffness(anm, pdb, [3,7], [0,7.5], filename='1gfl_3-7aa', loadToVMD=False)
writeVMDstiffness(anm, pdb, [3], [0,7], filename='1gfl_3aa')

