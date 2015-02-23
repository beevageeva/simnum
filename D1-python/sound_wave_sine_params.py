import numpy as np
from math import pi
from constants import z0, zf

phi0 = pi / 6.0 
k0 = 5
phi = phi0 - 2.0 * pi * z0 / (zf - z0)


def getSoundWaveFunction(k0, phi):
		def sinFunction(z):
			return np.sin(np.multiply((2.0 * pi * k0/(zf - z0)),z) + phi )
		return sinFunction
