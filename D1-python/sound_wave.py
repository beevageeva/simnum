import numpy as np
from math import pi

	#u(x) = periodicFunction (x)
	#csSign = 1 if wave travelling right, -1 if wave travelling left
	#after time t
	#u(x) = periodicFunction(x - csSign * cs * t)

class SoundWave:

	def __init__(self, cs, csSign, periodicFunction):
		self.cs = cs
		self.csSign = csSign
		self.periodicFunction = periodicFunction	
	

	def getWaveShape(self, z, t=0):
		#TODO  if t = 0 there is no need to getPeriodicX for all the elements, but I have to apply bounday conditions for first and last element which are not in the domain because of the cells, but it's a good practice to assure that all elements of z are in the domain(why shouldn't be ??)
		from common import getPeriodicXArray
		newz = z - self.csSign * self.cs * t
		z = getPeriodicXArray(newz)
		return self.periodicFunction(z)



#monochromatic sound wave with amplitude 1

class SoundWaveSine(SoundWave):
	
	#u(x) = sin (2*pi*x/wl + phi)
	#csSign = 1 if wave travelling right, -1 if wave travelling left
	#after time t
	#u(x) =  sin (2 * pi * (x - csSign * cs * t)/ wl  + phi)
	def __init__(self, cs,csSign, wl, phi):
		def sinFunction(z):
			return np.sin(np.multiply((2.0 * pi/wl),z) + phi )
		SoundWave.__init__(self, cs, csSign, sinFunction)

	
	
