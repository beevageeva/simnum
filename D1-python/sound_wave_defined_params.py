import numpy as np
from constants import z0, zf
from math import pi


#from scipy.special import jv
#def w(z):
#	return jv(0,  np.multiply(z, 2.2)) * (np.add(np.tanh((np.subtract(z,4.3))/0.3),1.0)) * (np.subtract(1.0,np.tanh((np.subtract(z,6.2))/0.3))) 
#

#gauss
#R = 0.05
#
#def w(z):
#	return np.exp(-(z - 0.5 *(z0 + zf))**2 / R) 

#wave packet

kf = 2.0 * pi/ (zf - z0)
k0 = 60.0
zc = z0 + 0.2 * (zf - z0)
W = 0.05
def w(z, nwav=k0):
  k = k0 * kf
  t2 = np.subtract(z,zc) ** 2 
  return np.multiply(np.exp(-np.divide(t2, W**2)), np.cos(k * (z-z0)))

def wFFTAn(k):
	t2 = np.subtract(k,zc) ** 2
	print("ft an ")
	print(t2)	 
	return np.exp(-np.divide(t2, W**2)) 


