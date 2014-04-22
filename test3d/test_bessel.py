from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from constants import zf_0, z0_0, zf_1, z0_1

import numpy as np
from common import getZArray, getDz0, getDz1
from math import pi
import math



from scipy.special import hankel1

func = lambda z: np.sqrt(z[0] **2 +  z[1] **2)

wl =  math.sqrt((zf_0 - z0_0)**2 + (zf_1 - z0_1)**2)
k = 2 * math.pi / wl

def w(z):
	print("z=")
	print(z)	
	return hankel1(0,k*func(z))


def derivW(z):
	return -hankel1(1,k*func(z)) * k


def derivZ0(f):
	dx = getDz0()
	n = len(f)
	print("n=%d"%n)
	res = np.zeros((n, n))
	for i in range(0,n):	
		for j in range(1,n-1):	
			#centered
			res[i][j] = (f[i][j+1] - f[i][j-1])/2 * dx
		res[i][0] = res[i][1]
		res[i][n-1] = res[i][n-2]
	return res	

def derivZ1(f):
	dx = getDz1()
	n = len(f)
	res = np.zeros((n, n))
	for j in range(0,n):	
		for i in range(1,n-1):	
			#centered
			res[i][j] = (f[i+1][j] - f[i-1][j])/2 * dx
		res[0][j] = res[1][j]
		res[n-1][j] = res[n-2][j]
	return res	

bw = np.blackman(10)
print("blackman")
print(bw)


z = getZArray()

#vals = w(z)

#num deriv
r = func(z)
f = w(z)
vals = z[0] / r * derivZ0(f) + z[1] / r * derivZ1(f)

#vals = derivW(z)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_wireframe(z[0], z[1], vals)
plt.draw()
plt.show()
