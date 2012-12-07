import sys 
import matplotlib.pyplot as plt

#plt.figure()
#plt.subplot(121)
plt.plotfile(sys.argv[1],(0,1))#,newfig=False)
plt.title('neuroCL')
#plt.subplot(122)
plt.plotfile(sys.argv[2],(0,1))
plt.title('compare')
plt.show()

