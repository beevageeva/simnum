import numpy as np
import sys,math
from constants import problemtype


if(not problemtype in ["sound_wave"]):
	print("%s not implemented" % problemtype)
	sys.exit(0)


def getV00():
	if(problemtype == 'sound_wave'):
		from sound_wave_params import v00
		return v00


def getInitialFunctionValues(z):
	if(problemtype == 'sound_wave'):
		from sound_wave_params import w
		wz = w(z)
		#!!!first and last element are not in domain of w because of the cells 
		#boundary conditions
		wz[0] = wz[len(wz) - 2]
		wz[len(wz)-1] = wz[1]
		return wz

def getInitialFunctionMaxZ(z):
	if(problemtype == 'sound_wave'):
		wz = getInitialFunctionValues(z)
		#all points are the same in this case for pressure, density and velocity: all in phase
		markPoint = z[np.argmax(wz)]
		return {'pres': markPoint, 'vel': markPoint, 'rho': markPoint}
		


def getInitialPresRhoVel(z):
	if(problemtype == 'sound_wave'):
		from constants import gamma
		from sound_wave_params import A, p00, rho00, v00, w
		cs00 = math.sqrt(gamma * p00 / rho00)
		wz = getInitialFunctionValues(z)
		return {'pres': p00 + gamma * p00 * A * wz, 'rho': rho00 + rho00 * A * wz, 'vel': v00 + cs00 * A * wz }
			
def getCs0():
	from constants import gamma
	if(problemtype == 'sound_wave'):
		from sound_wave_params import p00, rho00
		cs = math.sqrt(gamma *  p00 / rho00)
		return cs

def getRhoCurve(z, t):
	if(problemtype == 'sound_wave'):
		from sound_wave_params import A, w, v00
		from constants import z0, zf
		from math import pi
		cs = getCs0()
		newz = z - (cs + v00)* t
		perz = []
		from common import getPeriodicX
		for zval in newz:
			perz.append(getPeriodicX(zval))				
		res = A *  w(perz)	
		return res
			
def getVelCurve(z, t):
	return getRhoCurve(z, t)
			

def getPresCurve(z, t):
	from constants import gamma
	return gamma * getRhoCurve(z,t)	
	

def getRhoCurveNumeric(rho):
	if(problemtype == 'sound_wave'):
		from sound_wave_params import rho00
		return np.divide(np.subtract(rho,rho00),rho00)


def getPresCurveNumeric(p):
	if(problemtype == 'sound_wave'):
		from sound_wave_params import p00
		return np.divide(np.subtract(p,p00),p00)

def getVelCurveNumeric(v):
	if(problemtype == 'sound_wave'):
		cs00 = getCs0() 
		return np.divide(v,cs00)
