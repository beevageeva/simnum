import numpy as np
from constants import z0, zf


periodicType = "repeat" 
#periodicType = "refl" 
if periodicType == "repeat":

	def  lrBoundaryConditionsPresRho(array, skip=0):
		n = array.shape[0] - 1
		array = np.insert(array, 0,  array[n-skip])
		array = np.insert(array, n+2,  array[1+skip])
		return array

	lrBoundaryConditionsVel = lrBoundaryConditionsPresRho

	def getPeriodicX(xval, a=z0, b=zf):
		p = float(b - a)
		k = int((xval-a)/p)
		res = xval - k * p
		if(res < a):
			res+=p
		if(res > b):
			res-=p
		return res

elif periodicType == "refl":


	def lrBoundaryConditionsPresRho(array, skip=0):
		n = array.shape[0] - 1
		array = np.insert(array, 0,  2 * array[0] - array[1])
		array = np.insert(array, n+2, 2 * array[-1] - array[-2])
		return array

	def lrBoundaryConditionsVel(array, skip=0):
		#I already know
		#if(len(array)<1+skip):
		#	return
		n = array.shape[0] - 1
		if(skip==0):
			array = np.insert(array, 0,  -array[0])
			array = np.insert(array, n+2,  -array[-1])
		elif (skip==1):
			array[0] = 0
			array = np.insert(array, 0,  -array[2])
			array[-1] = 0
			array = np.insert(array, n+2,  -array[-2])
		return array
		

	def getPeriodicX2(xval, a=z0, b=zf):
		p = float(b - a)
		k = int((xval-a)/p)
		if k%2==1:
			res = 2.0*b - xval + k*p
		else:
			res = xval - k * p
		if(res < a):
			res+=p
		if(res > b):
			res-=p
		return res


	def getPeriodicX(xval, a=z0, b=zf):
		if xval > b:
			return 2.0*b - xval
		if xval < a:
			return 2.0*a -xval
		return xval


def getPeriodicXArray(xarray, a=z0, b=zf):
	res = []	
	for xval in xarray:
		res.append(getPeriodicX(xval, a, b))
	return np.array(res)
