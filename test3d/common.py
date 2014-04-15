from constants import nint, z0_0, z0_1, zf_0, zf_1
import numpy as np


def getDz0():
	return float(zf_0 - z0_0) / nint

def getDz1():
	return float(zf_1 - z0_1) / nint

def getZArray():
	dz0 = getDz0()
	dz1 = getDz1()
	a = np.linspace(z0_0-0.5 *dz0, zf_0+0.5*dz0, nint+2)
	b = np.linspace(z0_1-0.5 *dz1, zf_1+0.5*dz1, nint+2)
	return np.meshgrid(a, b)
	

def getPeriodicX(xval, a, b):
	p = float(b - a)
	k = int((xval-a)/p)
	res = xval - k * p
	if(res < a):
		res+=p
	if(res > b):
		res-=p
	return res


def getPeriodicXArray(xarray, a, b):
	res = []	
	for xval in xarray:
		res.append(getPeriodicX(xval, a, b))
	return np.array(res)

#xarray has 2 dim
def getPeriodicXArray2(xarray, a, b):
	print("xarray shape")
	print(xarray.shape)	
	res = np.zeros(xarray.shape)
	print("res shape")
	print(res.shape)	
	for i in range(0, xarray.shape[0]):
		for j in range(0, xarray.shape[1]):
			res[i][j] = getPeriodicX(xarray[i][j], a, b)
	return res

def getPeriodicX2(xval, a, b):
	p = b - a
	while (xval < a):
		xval+=p
	while (xval > b):
		xval-=p
	return xval

def getZIndex0(z):
	return int( float(nint)*(float(z) - z0_0)/(zf_0 - z0_0) )

def getZIndex1(z):
	return int( float(nint)*(float(z) - z0_1)/(zf_1 - z0_1) )

#assumes periodic function
def displacedPoint(z, c, t):
	newz = z + c * t
	periodicz = getPeriodicX(newz)
	return periodicz


def getSpeedPeriodic(newVal, oldVal, z0, zf, dt):
	if newVal < oldVal:
		newVal +=  zf - z0	
	return (newVal - oldVal) / dt

def getSpeedPeriodic0(newVal, oldVal, dt):
	return getSpeedPeriodic(newVal, oldVal, z0_0, zf_0, dt)
	
def getSpeedPeriodic1(newVal, oldVal, dt):
	return getSpeedPeriodic(newVal, oldVal, z0_1, zf_1, dt)


#creates an output directory called out_0, out_1, ... the first that does not exists
def createFolder(dirname_base="out"):
	import os
	dirExists = True
	i = 0
	dirname = "%s_%i" % (dirname_base, i)
	while os.path.exists(dirname):
		i +=1
		dirname = "%s_%i" % (dirname_base, i)
	os.mkdir(dirname)
	return dirname





