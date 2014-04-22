import numpy as np
import sys,math
from constants import gamma



def getV00():
	from sound_wave_params import v00
	return v00

def getP00():
	from sound_wave_params import p00
	return p00

def getRho00():
	from sound_wave_params import rho00
	return rho00

def getCs0():
	from sound_wave_params import p00, rho00
	cs = math.sqrt(gamma *  p00 / rho00)
	return cs

def getInitialPresRhoVel(z):
	from sound_wave_params import A, p00, rho00, v00
	cs00 = math.sqrt(gamma * p00 / rho00)
	from sound_wave_params import w,velFunc
	r =  A * w(z)
	#initial velocity
	v1 = v00 + cs00 *  r
	vel = np.dstack(velFunc(v1, z))
	#print("getInitialPresRhoVel vel shape")
	#print(vel.shape)
	#print("getInitialPresRhoVel vel shape end")
	return {'pres': p00 + gamma * p00 * r  , 'rho': rho00 + rho00 * r , 'vel': vel } 



def getAnRhoPresVel(z, t):
	from sound_wave_params import A, p00, rho00, v00, wAnMonochromatic,velFunc
	cs00 = math.sqrt(gamma * p00 / rho00)
	#c is phaseVelocity
	def wAnMonochromatic(z, t, c):
		return np.cos(-k * c * t) * w(z)
	r =  A * wAnMonochromatic(z, t, cs00)
	v1 = v00 + cs00 *  r
	vel = np.dstack(velFunc(v1, z))
	return {'pres': p00 + gamma * p00 * r  , 'rho': rho00 + rho00 * r , 'vel': vel } 
	



#boundary conditions:	
from sound_wave_params import periodicType
		
if periodicType == "repeat":

	def lrBoundaryConditionsPresRho(array, skip=0):
		#print("lrBoundaryCond PresRho and Vel repeat")
		#print("array shape1")
		#print(array.shape)
		n = array.shape[0] - 1
		#insert rows	
		array = np.insert(array, 0,  array[n-skip,:], axis = 0)
		array = np.insert(array, n+2,  array[1+skip,:], axis = 0)
		#insert columns	
		array = np.insert(array, 0, array[:,n-skip], axis = 1)
		array = np.insert(array, n+2, array[:,1+skip], axis = 1)
		#print("array shape2")
		#print(array.shape)
		return array

	lrBoundaryConditionsVel = lrBoundaryConditionsPresRho

elif periodicType == "refl" or periodicType == "diff":
	
#	#degree + 1 points needed
#	def polyfitArray(arr, degree, direction, resIndex):
#		if(direction == "y"):
#			dim = arr.shape[0]
#		else:
#			dim = arr.shape[1]
#		res = np.zeros(dim)
#		for i in range(0, dim):
#			xvalues = range(1, degree+2)	
#			if(direction == "y"):
#				yvalues = arr[i,0:degree+1]
#			elif direction == "x":
#				yvalues = arr[0:degree+1, i]
#			c = np.polyfit(xvalues, yvalues, degree)
#			p = np.poly1d(c)
#			res[i]= p(resIndex)
#		return res
		

	def lrBoundaryConditionsPresRho(array, skip=0):
		n = array.shape[0] - 1
		fr = 2 * array[0,:] - array[1,:]
		lr = 2 * array[-1,:] - array[-2,:]
		array = np.insert(array, 0, fr, axis = 0)
		array = np.insert(array, n+2, lr, axis = 0)
		fc = 2 * array[:,0] - array[:,1]
		lc = 2 * array[:,-1] - array[:,-2]
		array = np.insert(array, 0, fc, axis = 1)
		array = np.insert(array, n+2, lc, axis = 1)
		#print("array shape2")
		#print(array.shape)
		#print("end")
		return array

#		#array may be 2d or 3d (in case of uc)
#		#print("refl lrBoundaryCondPresRho")
#		#print("array shape1")
#		#print(array.shape)
#		#I already know
#		#if(len(array)<2):
#		#	return
#		n = array.shape[0] - 1
#		#first insert rows
#		if(skip == 0):
#			fr = 2 * array[0,:] - array[1,:]
#			lr = 2 * array[-1,:] - array[-2,:]
#			#print("fr shape")
#			#print(fr.shape)
#		else: 
#			degree = 1+skip
#			fr = polyfitArray(array[0:1+degree,:], degree,  "x", 0)
#			lr = polyfitArray(array[-1-degree:,:], degree, "x", degree+2)
#		array = np.insert(array, 0, fr, axis = 0)
#		array = np.insert(array, n+2, lr, axis = 0)
#		if(skip == 0):
#			fc = 2 * array[:,0] - array[:,1]
#			lc = 2 * array[:,-1] - array[:,-2]
#			#print("fc shape")
#			#print(fc.shape)
#			#if(array.ndim == 2):
#				#print("fc[0]=%4.3f" % fc[0])
#				#print("lc[0]=%4.3f" % lc[0])
#			#else:
#				#print("fc_x[0]=%4.3f" % fc[0][0])
#				#print("fc_y[0]=%4.3f" % fc[0][1])
#				#print("lc_x[0]=%4.3f" % lc[0][0])
#				#print("lc_y[0]=%4.3f" % lc[0][1])
#
#		else: 
#			fc = polyfitArray(array[:,0:degree+1], degree,  "y", 0)
#			lc = polyfitArray(array[:,-1-degree:], degree, "y", degree+2)
#		array = np.insert(array, 0, fc, axis = 1)
#		array = np.insert(array, n+2, lc, axis = 1)
#		
#		#print("array shape2")
#		#print(array.shape)
#		#print("end")
#		return array

if periodicType == "refl":
	def lrBoundaryConditionsVel(array, skip=0):
		#print("refl lrBoundaryCondVel")
		#print("array shape1")
		#print(array.shape)
		n = array.shape[0] - 1
		if(skip==0):
			array = np.insert(array, 0,  -array[0,:], axis = 0)
			array = np.insert(array, n+2,  -array[-1,:], axis = 0)
			array = np.insert(array, 0,  -array[:,0], axis = 1)
			array = np.insert(array, n+2,  -array[:,-1], axis = 1)
		elif (skip==1):
			array[0,:] = 0
			array = np.insert(array, 0,  -array[2,:], axis = 0)
			array[:,-1] = 0			
			array = np.insert(array, n+2,  -array[-2,:], axis = 0)
			array[:,0] = 0
			array = np.insert(array, 0,  -array[:,2], axis = 1)
			array[:,-1] = 0			
			array = np.insert(array, n+2,  -array[:,-2], axis = 1)
		#print("array shape2")
		#print(array.shape)
		#print("end")
		return array

elif periodicType == "diff":
	lrBoundaryConditionsVel = lrBoundaryConditionsPresRho


