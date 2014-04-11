import numpy as np
import sys,math
from constants import gamma
from sound_wave_params import functiontype, periodic




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
	
	if(periodic == "x"):
		arg = z[0] #argument of periodic function only x   
		if functiontype == "sine":
			from constants import z0_0 as z0Per
			from constants import zf_0 as zfPer
			wl = zfPer - z0Per
	elif periodic == "y":
		arg = z[1] #arg y  
		if functiontype == "sine":
			from constants import z0_1 as z0Per
			from constants import zf_1 as zfPer
			wl = zfPer - z0Per
#	elif periodic == "r":
#		#symmetric
#		arg = np.sqrt(z[0]**2 + z[1]**2)  
#	elif periodic == "d":
#		from sound_wave_diag_init import getPeriodicKx, getPeriodicKy
#		a = getPeriodicKx()
#		b = getPeriodicKy()
#		print("a = %E, b= %E" % (a, b))
#		arg = a * z[0] + b * z[1]  
	elif periodic == "d1":
		from constants import zf_0 as  zf #square
		k = 2.0 * math.pi / zf
		from common import getPeriodicXArray2
		arg = z[0] * k +  z[1] * k
		z0Per = 0 #coord had changed
		wl = zf * math.sqrt(2)

	if functiontype == "sine":
		from sound_wave_sine_params import phi0
		phi =  phi0 - 2.0 * math.pi * z0Per / wl
		r =  A * np.sin(np.multiply((2.0 * math.pi/wl),arg) + phi )

	elif functiontype == "defined":
		from sound_wave_defined_params import w
		r = w(arg)

	#initial velocity
	v1 = v00 + cs00 *  r
	if(periodic == "x"):
		vel = np.dstack((v1,np.zeros(v1.shape)))  
	elif periodic == "y":
		vel = np.dstack((np.zeros(v1.shape),v1))  
	elif periodic == "r" or periodic == "d1":
		#symmetric
		vel = np.dstack((v1 ,v1))   
	elif periodic == "d":
		vel = np.dstack((a * v1 ,b * v1))   
		
	elif periodic == "d2":
		vel = np.dstack((v1 ,-v1))   



	return {'pres': p00 + gamma * p00 * r  , 'rho': rho00 + rho00 * r , 'vel': vel } 

from sound_wave_params import periodic
#TODO why np.diag is not working with 3D array
#array must be 2d or 3d (I think it would work with 1d Array as numpy.diag is ONLY defined for both 1d and 2d arrays)
def secondDiag(array):
	if(array.ndim ==2):
		return  np.diag(np.fliplr(array))
	elif (array.ndim ==3):
		n = array.shape[0] - 1
		res = np.zeros((n+1, array.shape[2]))
		#I assume array.shape[0] == array.shape[1] -> square
		for i in range(0, n+1):
			for j in range(0, n+1):
				res[i] = array[i][n-j] 
		return res

#array must be 2d or 3d
def firstDiag(array):
	if(array.ndim ==2):
		return  np.diag(array)
	elif (array.ndim ==3):
		n = array.shape[0]
		res = np.zeros((n,array.shape[2]))
		#I assume array.shape[0] == array.shape[1] -> square
		for i in range(0, n):
			for j in range(0, n):
				res[i] = array[i][j] 
		return res

#array must be 1d or 2d
def interpMiddle(array):
	#TODO use numpy.interp 
	n = array.shape[0] 
	if(array.ndim ==1):
		res = np.zeros(2*n-1)  
	elif (array.ndim ==2):
		res = np.zeros((2*n-1,array.shape[1]))
	for i in range(0, n-1):
		res[2 * i] = array[i]
		res[2 * i + 1] = 0.5 * (array[i] + array[i+1])
	res[2 * n - 2] = array[n-1]
	return res
	
	
		
if periodic == "x" or periodic == "y" or periodic == "d1" or periodic == "r":

	def lrBoundaryConditions(array, skip=0):
		n = array.shape[0] - 1
		#insert rows	
		array = np.insert(array, 0,  array[n-skip,:], axis = 0)
		array = np.insert(array, n+2,  array[1+skip,:], axis = 0)
		#insert columns	
		array = np.insert(array, 0, array[:,n-skip], axis = 1)
		array = np.insert(array, n+2, array[:,1+skip], axis = 1)
		return array

#elif periodic == "r":
#	def lrBoundaryConditions(array, skip=0):
#		n = array.shape[0] - 1
#		array = np.insert(array, 0, array[0,:], axis = 0)
#		array = np.insert(array, n, array[n,:], axis = 0)
#		array = np.insert(array, 0, array[:,0], axis = 1)
#		array = np.insert(array, n, array[:,n], axis = 1)
#		return array

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

