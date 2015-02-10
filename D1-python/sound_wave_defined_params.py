import numpy as np
from constants import z0, zf
from math import pi


#ftype = "gauss"
ftype = "wavepacket"
#ftype = "bessel"

if ftype == "bessel":
	from scipy.special import jv
	def w(z):
		return jv(0,  np.multiply(z, 2.2)) * (np.add(np.tanh((np.subtract(z,4.3))/0.3),1.0)) * (np.subtract(1.0,np.tanh((np.subtract(z,6.2))/0.3))) 


elif ftype == "gauss":
	R = 0.05
	def w(z):
		return np.exp(-(z - 0.5 *(z0 + zf))**2 / R) 

elif ftype == "wavepacket":

	kf = 2.0 * pi/ (zf - z0)
	k0 = 60.0
	#k0 = 15.0#second exp of inhom
	#zc = z0 + 0.2 * (zf - z0)
	zc = z0 + (3.0/20.0)*(zf - z0)#second exp of inhom and first new
	W = 0.05
	#W = 0.25 #second exp of inhom
	def w(z, nwav=k0):
	  k = k0 * kf
	  t2 = np.subtract(z,zc) ** 2 
	  return np.multiply(np.exp(-np.divide(t2, W**2)), np.cos(k * (z-z0)))
	#
	def wFFTAn(k):
		a = 	1/(2 * np.sqrt(2) * np.sqrt(1/W**2))
		t3 = k*k0	* W**2 + 1j * 2*k0*z0
		t1 = np.cos(2 * k0 * zc) +  np.cosh(t3)  +   np.sinh(t3) + 1j *  np.sin(2 * k0 *zc)	
		t2 = np.cosh(0.25*(k**2 * W**2) + 0.5* k*k0 * W**2 + 0.25*(k0**2 * W**2) + 1j *( k0 * z0 -  k* zc + k0 * zc))  - np.sinh(0.25*(k**2 * W**2) + 0.5* k*k0 * W**2 + 0.25*(k0**2 * W**2) + 1j * ( k0 * z0 -  k* zc + k0 * zc)) 
		return a * t1 * t2 
	
	
	#	(1/(2 Sqrt[2] Sqrt[1/W^2]))(Cos[2 k0 zc] + 
	#   Cosh[k k0 W^2 + 2 I k0 z0] + I Sin[2 k0 zc] + 
	#   Sinh[k k0 W^2 + 2 I k0 z0]) (Cosh[(k^2 W^2)/4 + 1/2 k k0 W^2 + (
	#     k0^2 W^2)/4 + I k0 z0 - I k zc + I k0 zc] - 
	#   Sinh[(k^2 W^2)/4 + 1/2 k k0 W^2 + (k0^2 W^2)/4 + I k0 z0 - I k zc +
	#      I k0 zc])
	
