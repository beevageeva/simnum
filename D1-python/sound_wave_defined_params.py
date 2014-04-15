import numpy as np
from scipy.special import jv
from constants import z0, zf


#def w(z):
#	return jv(0,  np.multiply(z, 2.2)) * (np.add(np.tanh((np.subtract(z,4.3))/0.3),1.0)) * (np.subtract(1.0,np.tanh((np.subtract(z,6.2))/0.3))) 
#

#gauss
R = 0.05

def w(z):
	return np.exp(-(z - 0.5 *(z0 + zf))**2 / R) 
