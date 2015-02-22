import numpy as np
from constants import z0, zf
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from math import pi

#I may import it
from sound_wave_params import densFunc


sqrtDens = lambda z: np.sqrt(densFunc(z))


from common import getZArray

z = getZArray()
dz = z[1] - z[0]
numPoints = len(z)

plt.plot(z,sqrtDens(z), 'r-')
plt.draw()
plt.show()
