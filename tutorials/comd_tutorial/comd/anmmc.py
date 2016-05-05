from prody import *
from numpy import linalg, sqrt, loadtxt, savetxt, argmax, floor, abs, max, exp, mod, zeros
from random import random
import os.path
import sys

ar =[]
for arg in sys.argv:
        ar.append(arg)

initial_pdb=ar[1]
final_pdbn=ar[2]
anm_cut=float(ar[3])
devi=float(ar[4])

pdb = parsePDB(initial_pdb)
final_pdb = parsePDB(final_pdbn)

# Current Structure 
pdb_ca = pdb.ca
stepcutoff = 0.5 * (len(pdb_ca) ** 0.5)
N=1000000

# ANM calculation based on current
pdb_anm = ANM('pdb ca')
# Build Hessian Matrix
pdb_anm.buildHessian(pdb_ca, cutoff=anm_cut)
pdb_anm.calcModes(n_modes=None)

# Cumulative sum vector preparation for metropolis sampling
eigs = 1/sqrt(pdb_anm.getEigvals())
eigs_nn = zeros(eigs.shape)
eigs_nn[:] = eigs[:];
eigs = eigs / sum(eigs)
eigscumsum = eigs.cumsum()
U = pdb_anm.getEigvecs()

# Target Structure
final_pdb_ca = final_pdb.ca
# Number of residues on protein structure 		
size=pdb_ca.getResnums().shape[0]	
# Difference between current and final structure 
deviation = final_pdb_ca.getCoords() - pdb_ca.getCoords()
# Cutoff to check the structure deviated a lot	
scale_devi = devi
	# Scale factor for estimation of energy
scale_factor = sqrt(abs(scale_devi*min(pdb_anm.getEigvals())))
	# counts for metropolis sampling
count1 = 0 # Up-hill moves
count2 = 0 # Accepted up-hill moves
count3 = 0 # Down-hill moves
# read MC parameter from file
if os.path.isfile(initial_pdb + '_ratio.dat') and os.stat(initial_pdb + '_ratio.dat').st_size != 0:
	MCpara = loadtxt(initial_pdb + '_ratio.dat')
	accept_para = MCpara[4]
	if MCpara[1]>0.95:
		accept_para*=1.5
	elif MCpara[1]<0.85:
		accept_para/=1.5
	else:
		savetxt(initial_pdb + '_status.dat',[1])
else:
	accept_para = 0.1
# the best parameter is around 0.9 so that the parameters below than 0.85 and higher than 0.95 are not preferred and adjusted to limits.  
	
# difference from the target structure is defined as the energy and the minimum is zero. 
native_dist = buildDistMatrix(final_pdb_ca)
Ep = 0
dist = buildDistMatrix(pdb_ca)
for i in range(size-1):
	for j in range(i+1, size):
		Ep += (native_dist[i][j]-dist[i][j])**2

pdb_ca_ini = pdb_ca.copy()
ensemble = Ensemble()
ensemble_final = Ensemble()

#exit
# MC Loop 
for k in range(N):
	pdb_ca_temp = pdb_ca.copy()
	rand = random()	
	ID = argmax(rand<eigscumsum)		
	direction = 2*(random()>0.5)-1

	coords_temp = pdb_ca_temp.getCoords()
	coords_temp[0:,0] = coords_temp[0:,0] + direction * U[range(0,len(U),3),ID] * eigs_nn[ID] * scale_factor
	coords_temp[0:,1] = coords_temp[0:,1] + direction * U[range(1,len(U),3),ID] * eigs_nn[ID] * scale_factor
	coords_temp[0:,2] = coords_temp[0:,2] + direction * U[range(2,len(U),3),ID] * eigs_nn[ID] * scale_factor
	pdb_ca_temp.setCoords(coords_temp)
	
	En = 0
	dist = buildDistMatrix(pdb_ca_temp)
	for i in range(size-1):
		for j in range(i+1, size):
			En += (native_dist[i][j]-dist[i][j])**2
	#print exp(-(En-Ep)*accept_para)
	#print k
	#print accept_para
	if Ep > En:
		count3 += 1
		pdb_ca = pdb_ca_temp.copy()
		Ep = En
	elif exp(-(En-Ep)*accept_para)>random():
		pdb_ca = pdb_ca_temp.copy() 
		count1 += 1
		count2 += 1
		Ep = En
	else:
		count1 += 1
	print En
	if (mod(k,25)==0 and not(k==0)):
	# 	print k
	# 	print count2
	# 	print count3
	# 	print count1
	# 	print count2*1.0/count1
	# 	print accept_para
	 	if count2*1.0/count1>0.95:
	 		accept_para*=1.5;
	 	elif count2*1.0/count1<0.85:
	 		accept_para/=1.5
	coord_diff = pdb_ca.getCoords() - pdb_ca_ini.getCoords()
	print linalg.norm(coord_diff.ravel())
	if linalg.norm(coord_diff.ravel())	> stepcutoff: 
		break
		
	ensemble.addCoordset(pdb_ca.getCoords())
	
ensemble_final.addCoordset(pdb_ca.getCoords())
	
writeDCD(initial_pdb + '_' + final_pdbn + '_final_structure.dcd', ensemble_final)
ratios = [count2*1.0/N, count2*1.0/count1 if count1 != 0 else 0, count2, k, accept_para ]
savetxt(initial_pdb + '_ratio.dat', ratios)

