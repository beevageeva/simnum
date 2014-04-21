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
periodicType = "refl" 

timesZArgW = 2 #1(sine, gauss) or 2(wave packet)

if(timesZArgW == 1):
	functionType = "sine" 
	#functionType = "gauss" 
	#these sholuld be defined only in case of r argType: see below
	#velDir = "x"
	#velDir = "y"
	velDir = "d1"
	#argType = "x"
	#argType = "y"
	#argType = "r"
	argType = "d1" 
	
	#be sure 
	if(argType == "d1"):
		velDir = "d1"
	
	elif(argType == "x"):
		velDir = "x"
	
	elif(argType == "y"):
		velDir = "y"
		
	if(velDir == "x"):
		wl = zf_0 - z0_0
		velFunc = lambda x: (x,np.zeros(x.shape))
	elif velDir == "y":
		wl = zf_1 - z0_1
		velFunc = lambda x: (np.zeros(x.shape),x)
	elif velDir == "d1":
		wl =  math.sqrt((zf_0 - z0_0)**2 + (zf_1 - z0_1)**2)
		k1 = wl/(zf_0 - z0_0)
		k2 = wl/(zf_1 - z0_1)
		print("k1=%d, k2=%d" % (k1,k2))
		velFunc = lambda x: (k1 * x ,k2 * x)
	if(argType == "x"):
		func = lambda z: z[0] 
	elif argType == "y":
		func = lambda z: z[1] 
	elif argType == "d1":
		func = lambda z: z[0] * k1 + z[1] * k2 
	elif argType == "r":
		func = lambda z: np.sqrt(z[0] **2 +  z[1] **2)

	
	if functionType == "sine":
		phi0 = math.pi / 6.0 
		#sine function
		def w1(z, func,  wl, z0):
			return np.sin(np.multiply((2.0 * math.pi/wl),func(z-z0)) + phi0 )
	
	if functionType == "gauss":
		R = 0.05
		def w1(z, func, wl, z0):
			return np.exp(-(func(z -z0) - 0.5 * wl)**2 / R) 

	def w(z):
		n = len(z[0])
		z0 = getArrayZShape(z0_0,z0_1,n)
		return w1(z, func, wl, z0)


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
	velFunc = lambda x: ((k0[0]/modk) * x , (k0[1]/modk) * x)
