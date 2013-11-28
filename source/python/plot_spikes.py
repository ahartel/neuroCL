import sys 
import matplotlib.pyplot as plt


f = open(sys.argv[1])
spikes = {}
first = True
for line in f:
	if first:
		first = False
	else:
		line = line.rstrip()
		(time,number) = line.split(",")
		try:
			spikes[number].append(time)
		except KeyError:
			spikes[number] = []
			spikes[number].append(time)


for number,times in spikes.iteritems():
	plt.vlines(times,int(number)-0.45,int(number)+0.45)

plt.show()


