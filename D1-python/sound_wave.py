import numpy as np

	 


	#u(x) = periodicFunction (x)
	#after time t
	#u(x) = periodicFunction(x -  phaseVel * t)

class SoundWave:

	def __init__(self, A,  phaseVel,  periodicFunction):
		self.phaseVel = phaseVel
		self.A = A
		self.periodicFunction = periodicFunction	
	

	def getWaveShape(self, z, t=0):
		#TODO  if t = 0 there is no need to getPeriodicX for all the elements, but I have to apply bounday conditions for first and last element which are not in the domain because of the cells, but it's a good practice to assure that all elements of z are in the domain(why shouldn't be ??)
		from common import getPeriodicXArray
		newz = z - self.phaseVel * t
		z = getPeriodicXArray(newz)
		return self.A * self.periodicFunction(z)

	def printString(self):
		print("A=%E, phaseVel=%E" %(self.A, self.phaseVel) )


#monochromatic sound wave with amplitude A


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




	
