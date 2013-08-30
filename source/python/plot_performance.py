import sys
import matplotlib.pyplot as plt

average_by = 20
averages = []
cnt = 0
cur_average = 0
data = []

f = open(sys.argv[1])
for line in f:
    if cnt % average_by == 0 and cnt > 0:
        averages.append(cur_average/average_by)
        cur_average = 0

    line = line.rstrip('\n')
    data.append(float(line))
    cur_average += float(line)
    cnt += 1

f.close()

print averages

plt.subplot(211)
plt.plot(data)
plt.subplot(212)
plt.plot(averages)
plt.show()
