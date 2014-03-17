import numpy as np
from math import pi

	#u(x) = periodicFunction (x)
	#after time t
	#u(x) = periodicFunction(x -  cs * t)

class SoundWave:

	def __init__(self, A,  cs,  periodicFunction):
		self.cs = cs
		self.A = A
		self.periodicFunction = periodicFunction	
	

	def getWaveShape(self, z, t=0):
		#TODO  if t = 0 there is no need to getPeriodicX for all the elements, but I have to apply bounday conditions for first and last element which are not in the domain because of the cells, but it's a good practice to assure that all elements of z are in the domain(why shouldn't be ??)
		from common import getPeriodicXArray
		newz = z - self.cs * t
		z = getPeriodicXArray(newz)
		return self.A * self.periodicFunction(z)

	def printString(self):
		print("A=%E, cs=%E" %(self.A, self.cs) )


#monochromatic sound wave with amplitude 1

class SoundWaveSine(SoundWave):
	
	#u(x) = sin (2*pi*x/wl + phi)
	#after time t
	#u(x) =  sin (2 * pi * (x -  cs * t)/ wl  + phi)
	def __init__(self, A, cs, wl, phi):
		def sinFunction(z):
			return np.sin(np.multiply((2.0 * pi/wl),z) + phi )
		SoundWave.__init__(self, A,  cs, sinFunction)

class SuperpositionSoundWave(SoundWave):

	#array of simple waves
	def __init__(self, waveArray):
		self.waveArray = waveArray

	def getWaveShape(self, z, t=0):
		first = True
		for wave in self.waveArray:
			if first:
				result = wave.getWaveShape(z, t)
				first = False
			else:
				result = np.add(result, wave.getWaveShape(z, t))
		return result	

	def printString(self):
		print("SUPERPOSITION")
		for wave in self.waveArray:
			wave.printString()
		print("END")


	
