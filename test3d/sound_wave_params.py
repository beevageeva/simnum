import numpy as np
import math
from constants import gamma, z0_0, zf_0, z0_1, zf_1
from common import getArrayZShape


rho00 = 1.0
p00 = 1.0


v00 = 0.0

cs00 = math.sqrt(gamma * p00 / rho00)
print("cs00 = %4.3f" % cs00)

#v00 = - cs00 /  5.5
#v00 =  0.5 * cs00
#v00 = - cs00
#v00 =  cs00

A = 3.0 * 10.0 ** (-4)
#A = 5.0 * 10.0 ** (-2)

#periodicType = "repeat" # repeat| refl
#periodicType = "refl" 
periodicType = "diff" 
wType = "all"
#wType = "pot"

timesZArgW = 1 #1(sine, gauss) or 2(wave packet)

if(timesZArgW == 1):
	#functionType = "sine" 
	#functionType = "gauss" 
	functionType = "hankel" 
	#argType = "x"
	#argType = "y"
	argType = "r"
	#argType = "d1" 
	

	if(argType == "x"):
		wl = zf_0 - z0_0
		velFunc = lambda x,z: (x,np.zeros(x.shape))
		func = lambda z: z[0] 
		n = len(z[0])
		z0 = getArrayZShape(z0_0,z0_1,n)
	elif argType == "y":
		wl = zf_1 - z0_1
		velFunc = lambda x,z: (np.zeros(x.shape),x)
		func = lambda z: z[1] 
		n = len(z[0])
		z0 = getArrayZShape(z0_0,z0_1,n)
	elif argType == "d1":
		wl =  math.sqrt((zf_0 - z0_0)**2 + (zf_1 - z0_1)**2)
		k1 = wl/(zf_0 - z0_0)
		k2 = wl/(zf_1 - z0_1)
		print("k1=%d, k2=%d" % (k1,k2))
		velFunc = lambda x,z: (k1 * x ,k2 * x)
		func = lambda z: z[0] * k1 + z[1] * k2 
		n = len(z[0])
		z0 = getArrayZShape(z0_0,z0_1,n)
	elif argType == "r":
		wl =  math.sqrt((zf_0 - z0_0)**2 + (zf_1 - z0_1)**2)
		def velFunc(x, z):
			r = np.sqrt(z[0]**2 + z[1]**2)
			return x * z[0]/r , x * z[1] / r
		func = lambda z: np.sqrt(z[0] **2 +  z[1] **2)
		n = len(z[0])
		z0 = np.zeros((n,n))

	
	if functionType == "sine":
		phi0 = math.pi / 6.0 
		#sine function
		def w(z):
			return np.sin(np.multiply((2.0 * math.pi/wl),func(z-z0)) + phi0 )
	
	elif functionType == "gauss":
		R = 0.05
		def w(z):
			return np.exp(-(func(z -z0) - 0.5 * wl)**2 / R) 

	elif functionType == "hankel":
		#k = 1.0
		k = 2 * math.pi / wl
		from scipy.special import hankel1
		def w(z):
			return hankel1(0,k*func(z-z0))
		if wType == "pot":
			def derivW(z, z0):
				return -hankel1(1,k*func(z-z0)) * k
		

#wave packet
elif(timesZArgW == 2):
	kf = [2.0 * math.pi/ (zf_0 - z0_0), 2.0 * math.pi/ (zf_1 - z0_1)]
	k0 = [60.0,60.0]
	zc = [0.5 * (z0_0+ zf_0), 0.5 * (z0_1 + zf_1)] #in the middle
	W = 0.05
	def w(z, nwav=k0):
		k = [k0[0] * kf[0], k0[1] * kf[1]]
		n = len(z[0])
		#z0 = getArrayZShape(z0_0,z0_1,n)
		#direct operations
		#k0Reshaped = getArrayZShape(k0[0],k0[1],n)
		#zcReshaped = getArrayZShape(zc[0],zc[1],n)
		#t1 = np.subtract(z, zcReshaped)
		#t2 = t1[0,:,:] ** 2 + t1[1,:,:] ** 2
		#direct operations:
		t2 = np.subtract(z[0],zc[0]) ** 2 + np.subtract(z[1], zc[1]) ** 2
		#direct operations:
		#t3 = np.subtract(z, z0)
		#return np.multiply(np.exp(-np.divide(t2, W**2)), np.cos(t3[0] * k0[0] + t3[1] * k0[1]))
		#direct operations:
		sumt1 = np.multiply(np.subtract(z[0], z0_0), k[0])
		sumt2 = np.multiply(np.subtract(z[1], z0_1), k[1])
		return np.multiply(np.exp(-np.divide(t2, W**2)), np.cos(np.add(sumt1, sumt2)))

	modk = math.sqrt(k0[0] ** 2 + k0[1]**2)
	velFunc = lambda x, z: ((k0[0]/modk) * x , (k0[1]/modk) * x)
