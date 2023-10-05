from prody import *

PDBfile = '5kqm'
coords = parsePDB(PDBfile)

waterBridges_chain = calcWaterBridges(coords)

print(waterBridges_chain)