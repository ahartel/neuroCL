import sys 
import matplotlib.pyplot as plt

for argument in range(1,len(sys.argv),1):
	plt.plotfile(sys.argv[argument],(0,1),linestyle='',marker='|',markersize=10.)

plt.show()


