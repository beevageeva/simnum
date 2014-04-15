import numpy as np
from math import pi

functionType = "sine" 
#functionType = "gauss" 

if functionType == "sine":
	phi0 = pi / 6.0 
	#sine function
	def w(z, wl, z0):
		return np.sin(np.multiply((2.0 * pi/wl),(z-z0)) + phi0 )

if functionType == "gauss":
	R = 0.05
	def w(z, wl, z0):
		return np.exp(-(z -z0 - 0.5 * wl)**2 / R) 
