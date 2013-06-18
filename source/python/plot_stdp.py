import sys 
import matplotlib.pyplot as plt


for argument in range(1,len(sys.argv),2):
    print (0,int(sys.argv[argument+1])*2+1,int(sys.argv[argument+1])*2+2)
    plt.plotfile(sys.argv[argument],(0,int(sys.argv[argument+1])*2+1,int(sys.argv[argument+1])*2+2))

plt.show()


