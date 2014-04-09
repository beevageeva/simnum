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

	print("initcond_soundwave z[0]")
	print(z[0])
	print("initcond_soundwave z[1]")
	print(z[1])
	
	if(periodic == "x"):
		arg = z[0] #argument of periodic function only x   
	elif periodic == "y":
		arg = z[1] #arg y   2
	elif periodic == "r":
		#symmetric
		arg = np.sqrt(z[0]**2 + z[1]**2)  
	elif periodic == "d":
		from sound_wave_diag_init import getPeriodicKx, getPeriodicKy
		a = getPeriodicKx()
		b = getPeriodicKy()
		print("a = %E, b= %E" % (a, b))
		arg = a * z[0] + b * z[1]  

	if functiontype == "sine":
		from sound_wave_sine_params import wl, phi
		r =  A * np.sin(np.multiply((2.0 * math.pi/wl),arg) + phi )
	elif functiontype == "defined":
		from sound_wave_defined_params import w
		r = w(arg)
	v1 = v00 + cs00 *  r

	if(periodic == "x"):
		vel = np.dstack((v1,np.zeros(v1.shape)))  
	elif periodic == "y":
		vel = np.dstack((np.zeros(v1.shape),v1))  
	elif periodic == "r":
		#symmetric
		vel = np.dstack((v1 ,v1))   
	elif periodic == "d":
		vel = np.dstack((a * v1 ,b * v1))   



	return {'pres': p00 + gamma * p00 * r  , 'rho': rho00 + rho00 * r , 'vel': vel } 

from sound_wave_params import periodic
			
if periodic == "x" or periodic == "y" or periodic == "d":

	def lrBoundaryConditions(array, skip=0):
		n = array.shape[0] - 1
		#insert rows	
		array = np.insert(array, 0,  array[n-skip,:], axis = 0)
		array = np.insert(array, n+2,  array[1+skip,:], axis = 0)
		#insert columns	
		array = np.insert(array, 0, array[:,n-skip], axis = 1)
		array = np.insert(array, n+2, array[:,1+skip], axis = 1)
		return array

elif periodic == "r":

	def secondDiag(array):
		if(array.ndim ==2):
			return  np.diag(np.fliplr(array))
		elif (array.ndim ==3):
			n = array.shape[0] - 1
			res = np.zeros((n+1, 2))
			for i in range(0, n+1):
				for j in range(0, n+1):
					res[i] = array[i][n-j] 
			return res
				

	def lrBoundaryConditions(array, skip=0):
		n = array.shape[0] - 1
		sdiag = secondDiag(array)	
		#sdiagonal is diag right top  -> left bottom
		sdiagrev = sdiag[::-1] #reverse order for the rows
		#insert rows	
		array = np.insert(array, 0,  sdiagrev, axis = 0)
		array = np.insert(array, n+2,  sdiagrev, axis = 0)
		#insert in the diagonal
		print("sdaig = ")
		print(sdiag) 
		print("sd shape = ")
		print(sdiag.shape) 
		sdiag = np.insert(sdiag, 0,  sdiag[n-skip], axis=0)
		sdiag = np.insert(sdiag, n+2,  sdiag[1+skip], axis=0)
		print("sdaig2 = ")
		print(sdiag) 
		print("sd shape = ")
		print(sdiag.shape) 
		#insert columns	
		array = np.insert(array, 0, sdiag, axis = 1)
		array = np.insert(array, n+2, sdiag, axis = 1)
		return array

elif periodic == "d":
	def lrBoundaryConditions(array, skip=0):
		return array


