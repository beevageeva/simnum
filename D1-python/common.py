from constants import nint, z0, zf
import numpy as np


def getDz():
	return float(zf - z0) / nint

def getZArray():
	dz = getDz()
	return np.linspace(z0-0.5 *dz, zf+0.5*dz, nint+2)
	


def zeroArray():
	return np.zeros(nint+2)

def getPeriodicX(xval, a=z0, b=zf):
	p = float(b - a)
	k = int((xval-a)/p)
	res = xval - k * p
	if(res < a):
		res+=p
	if(res > b):
		res-=p
	return res


def getPeriodicXArray(xarray, a=z0, b=zf):
	res = []	
	for xval in xarray:
		res.append(getPeriodicX(xval, a, b))
	return np.array(res)

def getPeriodicX2(xval, a=z0, b=zf):
	p = b - a
	while (xval < a):
		xval+=p
	while (xval > b):
		xval-=p
	return xval

def getZIndex(z):
	return int( float(nint)*(z - z0)/(zf - z0) )

def displacedPoint(z, c, t):
	newz = z + c * t
	periodicz = getPeriodicX(newz)
	return periodicz


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

