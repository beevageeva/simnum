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
