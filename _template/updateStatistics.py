import numpy as np
overall = np.loadtxt('statistics.txt', dtype=str)
overall = np.char.replace(overall, ',','').astype(np.int)

statistics_bak = '../../../statistics.txt.bak_dont_remove'
overall_bak = np.loadtxt(statistics_bak, dtype=str)
overall_bak = np.char.replace(overall_bak, ',','').astype(np.int)

if overall_bak > overall:
    overall = overall_bak

try:
	last6 = np.loadtxt('last-6hour.txt')
except ValueError:
	last6 = 0

overall = overall + last6.astype(int)
if overall_bak < overall:
    f = open(statistics_bak,'w')
    f.write("{:,}".format(int(overall)) + "\n")
    f.close()

f = open('statistics.txt','w')
f.write("{:,}".format(int(overall)) + "\n")
f.close()

f = open('statistics.json', 'w')
f.write('{"stat": "' + '{:,}'.format(int(overall)) + '"}\n')
f.close()
