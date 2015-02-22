import numpy as np
from math import pi
phi0 = pi / 6.0 


def getSoundWaveFunction(wl, phi):
		def sinFunction(z):
			return np.sin(np.multiply((2.0 * pi/wl),z) + phi )
		return sinFunction
