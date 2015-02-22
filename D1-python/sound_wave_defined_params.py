import numpy as np

#here go all other definitions

#bessel
def getSoundWaveFunction():
	from scipy.special import jv
	def w(z):
		return jv(0,  np.multiply(z, 2.2)) * (np.add(np.tanh((np.subtract(z,4.3))/0.3),1.0)) * (np.subtract(1.0,np.tanh((np.subtract(z,6.2))/0.3))) 
	return w
