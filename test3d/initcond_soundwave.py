import numpy as np
import sys,math
from constants import gamma, z0_0, zf_0, z0_1, zf_1
from sound_wave_params import  argType, velDir



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

	#print("initcond_soundwave z[0]")
	#print(z[0])
	#print("initcond_soundwave z[1]")
	#print(z[1])
	
	if(velDir == "x"):
		wl = zf_0 - z0_0
		z0Per = z0_0
	elif velDir == "y":
		wl = zf_1 - z0_1
		z0Per = z0_1
	elif velDir == "d1":
		wl =  math.sqrt((zf_0 - z0_0)**2 + (zf_1 - z0_1)**2)
		k1 = wl/(zf_0 - z0_0)
		k2 = wl/(zf_1 - z0_1)
		z0Per = z0_0 * k1 + z0_1 * k2 
		print("z0Per = %E, wl = %E" % (z0Per, wl))

	if(argType == "x"):
		arg = z[0] #argument of periodic function only x   
	elif argType == "y":
		arg = z[1] #arg y  
	elif argType == "d1":
		arg = z[0] * k1 +  z[1] * k2
	elif argType == "r":
		arg = np.sqrt(z[0] **2 +  z[1] **2)

	from sound_wave_function import w
	r =  A * w(arg, wl, z0Per)


	#initial velocity
	v1 = v00 + cs00 *  r
	if(velDir == "x"):
		vel = np.dstack((v1,np.zeros(v1.shape)))  
	elif velDir == "y":
		vel = np.dstack((np.zeros(v1.shape),v1))  
	elif velDir == "d1":
		vel = np.dstack((k1 * v1 ,k2 * v1))   

	#print("******************Vel")
	#print(vel)	

	return {'pres': p00 + gamma * p00 * r  , 'rho': rho00 + rho00 * r , 'vel': vel } 

##TODO why np.diag is not working with 3D array
##array must be 2d or 3d (I think it would work with 1d Array as numpy.diag is ONLY defined for both 1d and 2d arrays)
#def secondDiag(array):
#	if(array.ndim ==2):
#		return  np.diag(np.fliplr(array))
#	elif (array.ndim ==3):
#		n = array.shape[0] - 1
#		res = np.zeros((n+1, array.shape[2]))
#		#I assume array.shape[0] == array.shape[1] -> square
#		for i in range(0, n+1):
#			for j in range(0, n+1):
#				res[i] = array[i][n-j] 
#		return res
#
##array must be 2d or 3d
#def firstDiag(array):
#	if(array.ndim ==2):
#		return  np.diag(array)
#	elif (array.ndim ==3):
#		n = array.shape[0]
#		res = np.zeros((n,array.shape[2]))
#		#I assume array.shape[0] == array.shape[1] -> square
#		for i in range(0, n):
#			for j in range(0, n):
#				res[i] = array[i][j] 
#		return res
#
##array must be 1d or 2d
#def interpMiddle(array):
#	#TODO use numpy.interp 
#	n = array.shape[0] 
#	if(array.ndim ==1):
#		res = np.zeros(2*n-1)  
#	elif (array.ndim ==2):
#		res = np.zeros((2*n-1,array.shape[1]))
#	for i in range(0, n-1):
#		res[2 * i] = array[i]
#		res[2 * i + 1] = 0.5 * (array[i] + array[i+1])
#	res[2 * n - 2] = array[n-1]
#	return res
	
	
from sound_wave_params import periodicType
		
if periodicType == "repeat":

	def lrBoundaryConditions(array, skip=0):
		n = array.shape[0] - 1
		#insert rows	
		array = np.insert(array, 0,  array[n-skip,:], axis = 0)
		array = np.insert(array, n+2,  array[1+skip,:], axis = 0)
		#insert columns	
		array = np.insert(array, 0, array[:,n-skip], axis = 1)
		array = np.insert(array, n+2, array[:,1+skip], axis = 1)
		return array

elif periodicType == "diff":
	
	#degree + 1 points needed
	def polyfitArray(arr, degree, dim, direction):
		res = np.zeros(dim)
		for i in range(0, dim):
			xvalues = np.range(1, degree+2)	
			if(direction == "y"):
				yvalues = arr[i,0:degree+1]
			elif direction == "x":
				yvalues = arr[0:degree+1, i]
			c = np.polyfit(xvalues, yvalues, degree)
			p = np.poly1d(c)
			res.append(p(0))
		return res
		

	def lrBoundaryConditions(array, skip=0):
		n = array.shape[0] - 1
		fr = polyfitArray(array[0:1+skip,:], skip, n, "x")
		array = np.insert(array, 0, fr, axis = 0)
		lr = polyfitArray(array[-1-skip:,:], skip, n, "x")
		array = np.insert(array, n, lr, axis = 0)
		fc = polyfitArray(array[:,0:1+skip], skip, n+1, "y")
		array = np.insert(array, 0, fc, axis = 1)
		lc = polyfitArray(array[:,-1-skip:], skip, n+1, "y")
		array = np.insert(array, n, lc, axis = 1)
		return array

elif periodicType == "same":	
	def lrBoundaryConditions(array, skip=0):
		n = array.shape[0] - 1
		array = np.insert(array, 0, array[0,:], axis = 0)
		array = np.insert(array, n, array[n,:], axis = 0)
		array = np.insert(array, 0, array[:,0], axis = 1)
		array = np.insert(array, n, array[:,n], axis = 1)
		return array
#elif periodic == "r":
#
#	def lrBoundaryConditions(array, skip=0):
#		n = array.shape[0] - 1
#		sdiag = secondDiag(array)	
#		#sdiagonal is diag right top  -> left bottom
#		sdiagrev = sdiag[::-1] #reverse order for the rows
#		#insert rows	
#		array = np.insert(array, 0,  sdiagrev, axis = 0)
#		array = np.insert(array, n+2,  sdiagrev, axis = 0)
#		#insert in the diagonal
#		print("sdaig = ")
#		print(sdiag) 
#		print("sd shape = ")
#		print(sdiag.shape) 
#		sdiag = np.insert(sdiag, 0,  sdiag[n-skip], axis=0)
#		sdiag = np.insert(sdiag, n+2,  sdiag[1+skip], axis=0)
#		print("sdaig2 = ")
#		print(sdiag) 
#		print("sd shape = ")
#		print(sdiag.shape) 
#		#insert columns	
#		array = np.insert(array, 0, sdiag, axis = 1)
#		array = np.insert(array, n+2, sdiag, axis = 1)
#		return array

#elif periodic == "d1":
#	def lrBoundaryConditions(array, skip=0):
#		print("ARRAY shape")
#		print(array.shape)
#		n = array.shape[0] - 1
#		print("n=%d"%n)
#		fdiag = firstDiag(array)	
#		#insert bc in the first diag
#		fdiag = np.insert(fdiag, 0,  fdiag[n-skip], axis=0)
#		fdiag = np.insert(fdiag, n+2,  fdiag[1+skip], axis=0)
#		print("diag shape")
#		print(fdiag.shape)
#		fdiagInterp = interpMiddle(fdiag)
#		print("diag shape interpolated")
#		print(fdiagInterp.shape)
#		print("diag interpolated")
#		print(fdiagInterp)
#		
#		frow = fdiagInterp[1:n+2, ...]
#		lrow = fdiagInterp[n+3:2*n+4, ...]
#		#lrow = fdiagInterp[n+4:-1, ...]
#		
#		print("lrBoundaryConditions array shape")
#		print(array.shape)
#		print("frow shape")
#		print(frow.shape)	
#	
#		#insert rows	
#		array = np.insert(array, 0,  frow, axis = 0)
#		array = np.insert(array, n+2,  lrow, axis = 0)
#
#		fcol = fdiagInterp[0:n+3, ...]
#		lcol = fdiagInterp[n+2:, ...]
#		print("array shape bfore ins cols")
#		print(array.shape)
#		print("fcol shape")
#		print(fcol.shape)
#		print("lcol shape")
#		print(lcol.shape)
#
#		
#		#insert columns	
#		array = np.insert(array, 0, fcol, axis = 1)
#		array = np.insert(array, n+2, lcol, axis = 1)
#		return array
#		return array
#

