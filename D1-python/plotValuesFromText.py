import numpy as np
from constants import z0, zf, gamma
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import sys

data = np.loadtxt(sys.argv[1])
x = data[:,0]
y = data[:,1]
#plt.plot(x, y, linestyle="None", marker="o", color="r")

#plt.ylim(0, 0.004)
#plt.ylim(0, 0.002)
#plt.ylim(0, 0.0004)

plt.clf()
plt.xlabel("time")	

print(x[len(x)-1])
print("max")
print(np.max(x))
print("amax")
print(np.argmax(x))

print("max y = %e at x = %e" % (np.max(y),x[np.argmax(y)] ))
#plt.plot(x, y, 'bo', markersize=1)
plt.plot(x, y)



plt.draw()
plt.show()

