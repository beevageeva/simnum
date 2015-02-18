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

def getSoundWavePacketFFTAnalyticalAbs(k0, zc, W):
	from constants import z0, zf
	def analyticFFT(k):
		from math import sqrt	
		real = (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(z0 - zf))**2))/4)*sqrt(pi)*np.cos(2.0*pi*zc*(k + k0/(z0 - zf)))*np.cos((2.0*k0*pi*z0)/(z0 - zf)))/(2.0*sqrt(W**(-2))) + (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(-z0 + zf))**2))/4)*sqrt(pi)*np.cos((2.0*k0*pi*z0)/(z0 - zf))*np.cos((4.0*k0*pi*z0)/(-z0 + zf) - 2.0*pi*zc*(k + k0/(-z0 + zf))))/(2.0*sqrt(W**(-2))) + (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(z0 - zf))**2))/4)*sqrt(pi)*np.sin(2.0*pi*zc*(k + k0/(z0 - zf)))*np.sin((2.0*k0*pi*z0)/(z0 - zf)))/(2.0*sqrt(W**(-2))) - (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(-z0 + zf))**2))/4)*sqrt(pi)*np.sin((2.0*k0*pi*z0)/(z0 - zf))*np.sin((4.0*k0*pi*z0)/(-z0 + zf) - 2.0*pi*zc*(k + k0/(-z0 + zf))))/(2.0*sqrt(W**(-2)))  
		imag = (-(np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(z0 - zf))**2))/4)*sqrt(pi)*np.cos((2.0*k0*pi*z0)/(z0 - zf))*np.sin(2.0*pi*zc*(k + k0/(z0 - zf))))/(2.0*sqrt(W**(-2))) + (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(z0 - zf))**2))/4)*sqrt(pi)*np.cos(2.0*pi*zc*(k + k0/(z0 - zf)))*np.sin((2.0*k0*pi*z0)/(z0 - zf)))/(2.0*sqrt(W**(-2))) + (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(-z0 + zf))**2))/4)*sqrt(pi)*np.cos((4.0*k0*pi*z0)/(-z0 + zf) - 2.0*pi*zc*(k + k0/(-z0 + zf)))*np.sin((2.0*k0*pi*z0)/(z0 - zf)))/(2.0*sqrt(W**(-2))) + (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(-z0 + zf))**2))/4)*sqrt(pi)*np.cos((2.0*k0*pi*z0)/(z0 - zf))*np.sin((4.0*k0*pi*z0)/(-z0 + zf) - 2.0*pi*zc*(k + k0/(-z0 + zf))))/(2.0*sqrt(W**(-2))))
		return np.sqrt(real**2+imag**2)
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




	
