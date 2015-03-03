import numpy as np
from constants import z0, zf
from math import pi,sqrt


k0 = 60.0
#k0 = 15.0#second exp of inhom
#zc = z0 + 0.2 * (zf - z0)
zc = z0 + (3.0/20.0)*(zf - z0)#second exp of inhom and first new
#W = 0.01
W = 0.05
#W = 0.25 #second exp of inhom

def getSoundWaveGaussFunction(zc, W):
	def gaussFunction(z):
		t2 = np.subtract(z,zc) ** 2
		return np.exp(-np.divide(t2, W**2))	 
	return gaussFunction


def getSoundWaveFunction(k0, zc, W):
		def gaussPacketFunction(z):
			return np.multiply(getSoundWaveGaussFunction(zc, W)(z),  np.cos(2.0 * pi * k0 * (z-z0)/ (zf - z0) ) )
		return gaussPacketFunction

def getSoundWaveFFTAnalytical(k0, zc, W):
	from constants import z0, zf
	def analyticFFT(k):
		return ((np.exp(-((pi*(k0**2*pi*W**2 - 2*k0*(k*pi*W**2 - 1j*z0 + 1j*zc)*(z0 - zf) + k*(k*pi*W**2 + (2*1j)*zc)*(z0 - zf)**2))/(z0 - zf)**2)) +  np.exp(-((pi*(k0**2*pi*W**2 + 2*k0*(k*pi*W**2 - 1j*z0 + 1j*zc)*(z0 - zf) + k*(k*pi*W**2 + (2*1j)*zc)*(z0 - zf)**2))/(z0 - zf)**2)))*sqrt(pi)*W)/2
	return analyticFFT
