import numpy as np
from math import pi

def getSoundWaveGaussFunction(R):
	def w(z):
		return np.exp(-(z - 0.5 *(z0 + zf))**2 / R)
	return w;
	 
def getSoundWaveBesselFunction():
	from scipy.special import jv
	def w(z):
		return jv(0,  np.multiply(z, 2.2)) * (np.add(np.tanh((np.subtract(z,4.3))/0.3),1.0)) * (np.subtract(1.0,np.tanh((np.subtract(z,6.2))/0.3))) 
	return w

def getSoundWaveSineFunction(wl, phi):
		def sinFunction(z):
			return np.sin(np.multiply((2.0 * pi/wl),z) + phi )
		return sinFunction

def getSoundWavePacketFunction(k0, zc, W):
		from constants import z0, zf
		kf = 2.0 * pi/ (zf - z0)
		k = k0 * kf
		def gaussPacketFunction(z):
		  t2 = np.subtract(z,zc) ** 2 
		  return np.multiply(np.exp(-np.divide(t2, W**2)), np.cos(k * (z-z0)))
		return gaussPacketFunction

def getSoundWavePacketFFTAnalytical(k0, zc, W):
	from constants import z0, zf
	def analyticFFT(k):
		a = 	1/(2 * np.sqrt(2) * np.sqrt(1/W**2))
		t3 = k*k0	* W**2 + 1j * 2*k0*z0
		t1 = np.cos(2 * k0 * zc) +  np.cosh(t3)  +   np.sinh(t3) + 1j *  np.sin(2 * k0 *zc)	
		t2 = np.cosh(0.25*(k**2 * W**2) + 0.5* k*k0 * W**2 + 0.25*(k0**2 * W**2) + 1j *( k0 * z0 -  k* zc + k0 * zc))  - np.sinh(0.25*(k**2 * W**2) + 0.5* k*k0 * W**2 + 0.25*(k0**2 * W**2) + 1j * ( k0 * z0 -  k* zc + k0 * zc)) 
		return  a * t1 * t2  
	return analyticFFT

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




	
