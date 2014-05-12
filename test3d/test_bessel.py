from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from constants import zf_0, z0_0, zf_1, z0_1, nint

import numpy as np
from common import getZArray, getDz0, getDz1, getZIndex0, getZIndex1
from math import pi
import math



from scipy.special import hankel1

func = lambda z: np.sqrt(z[0] **2 +  z[1] **2)

wl =  math.sqrt((zf_0 - z0_0)**2 + (zf_1 - z0_1)**2)
k = 3*  2 * math.pi / wl

print("k=%E" % k)

np.set_printoptions(threshold='nan')

	

def derivW(z):
	return -hankel1(1,k*func(z))  * k


def gradNum(f):
	from common import getDz0, getDz1
	dz0 = getDz0()
	dz1 = getDz1()
	return np.gradient(f,dz0, dz1 )






#bw = np.blackman(10)
#print("blackman")
#print(bw)


z = getZArray()

#print("z=")
#print(z)	

#vals = w(z)


#print("f=")
#print(f)
#vals = w1(z,k)
#vals = f
#vals = k * z[0] / r * derivZ0(f) + k * z[1] / r * derivZ1(f)
#vals =  z[0] / r * derivZ0(f) + z[1] / r * derivZ1(f)

#vals = derivW(z)

#vals = np.outer(np.blackman(nint+2), np.blackman(nint+2))

from sound_wave_params import w
from sound_wave_params import symDerivW

f = w(z)
#vals = w(z)
vals = symDerivW(z)[0]
#vals = gradNum(f)[0]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_wireframe(z[0], z[1], vals)
plt.draw()
plt.show()
