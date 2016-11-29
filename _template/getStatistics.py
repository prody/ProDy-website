from vanity import count_downloads
from time import sleep
while True:
	f = open('statistics.txt','w')
	f.write(str(count_downloads('ProDy', verbose=False)) + "\n")
	f.close()
	sleep(60)


