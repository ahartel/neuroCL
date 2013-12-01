import sys
import numpy as np 
import matplotlib.pyplot as plt

f = open(sys.argv[1])
header = f.readline().rstrip()
header = header.split(",")[1:]
f.close()

data = np.loadtxt(sys.argv[1],skiprows=1,delimiter=",")

for x in range(1,len(data[0])):
	plt.plot(data[:,0],data[:,x])
plt.legend(header,loc=(0,0))

plt.show()

