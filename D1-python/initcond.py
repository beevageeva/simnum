import numpy as np
import sys,math
from constants import problemtype, gamma


wave = None
if(not problemtype in ["sound_wave"]):
	print("problemtype %s not implemented" % problemtype)
	sys.exit(0)
else:
	from constants import functiontype
	from sound_wave_params import p00, rho00, csSign, v00
	from sound_wave import SoundWave, SoundWaveSine
	cs00 = math.sqrt(gamma * p00 / rho00)
	if functiontype == "sine":
		from sound_wave_sine_params import wl, phi
		wave = SoundWaveSine(cs00 + v00, csSign, wl, phi)
	elif functiontype == "defined":
		from sound_wave_defined_params import w
		wave = SoundWave(cs00+v00, csSign, w)	
	else:	
		print("functiontype %s not implemented" % functiontype)
		sys.exit(0)


def getV00():
	if(problemtype == 'sound_wave'):
		from sound_wave_params import v00
		return v00


def getInitialFunctionValues(z):
#	if(problemtype == 'sound_wave'):
#		from sound_wave_params import w
#		wz = w(z)
#		#!!!first and last element are not in domain of w because of the cells 
#		#boundary conditions
#		wz[0] = wz[len(wz) - 2]
#		wz[len(wz)-1] = wz[1]
#		return wz
	return wave.getWaveShape(z)	

def getInitialFunctionMaxZ(z):
	if(problemtype == 'sound_wave'):
		wz = getInitialFunctionValues(z)
		#all points are the same in this case for pressure, density and velocity: all in phase
		markPoint = z[np.argmax(wz)]
		return {'pres': markPoint, 'vel': markPoint, 'rho': markPoint}
		


def getInitialPresRhoVel(z):
	if(problemtype == 'sound_wave'):
		from constants import gamma
		from sound_wave_params import A, p00, rho00, v00
		cs00 = math.sqrt(gamma * p00 / rho00)
		wz = getInitialFunctionValues(z)
		return {'pres': p00 + gamma * p00 * A * wz, 'rho': rho00 + rho00 * A * wz, 'vel': v00 + csSign * cs00 * A * wz } 
		#return {'pres': p00 + gamma * p00 * A * wz, 'rho': rho00 + rho00 * A * wz, 'vel': v00 + cs00 * A * wz } #INITIAL
		#return {'pres': p00 - gamma * p00 * A * wz, 'rho': rho00 - rho00 * A * wz, 'vel': v00 + cs00 * A * wz }
		#return {'pres': p00 + gamma * p00 * A * wz, 'rho': rho00 + rho00 * A * wz, 'vel': v00 - cs00 * A * wz }
		#return {'pres': p00 - gamma * p00 * A * wz, 'rho': rho00 - rho00 * A * wz, 'vel': v00 - cs00 * A * wz }
		#return {'pres': p00 + gamma * p00 * A * wz, 'rho': rho00 - rho00 * A * wz, 'vel': v00 - cs00 * A * wz }  #INVALID, (rho1>rho2 <=> p1 > p2) <=> d_p / d_rho >0 (=cs**2)
			
def getCs0():
	from constants import gamma
	if(problemtype == 'sound_wave'):
		from sound_wave_params import p00, rho00
		cs = math.sqrt(gamma *  p00 / rho00)
		return cs

def getRhoCurve(z, t):
#	if(problemtype == 'sound_wave'):
#		from sound_wave_params import A, w, v00
#		from constants import z0, zf
#		from math import pi
#		cs = getCs0()
#		newz = z - (cs + v00)* t
#		perz = []
#		from common import getPeriodicX
#		for zval in newz:
#			perz.append(getPeriodicX(zval))				
#		res = A *  w(perz)	
	from sound_wave_params import A
	res =  A * wave.getWaveShape(z, t)	
	return res
			
def getVelCurve(z, t):
	return csSign * getRhoCurve(z, t)
			

def getPresCurve(z, t):
	from constants import gamma
	return gamma * getRhoCurve(z,t)	

#we can send curves calculated previously as parameters in order not to calculate them twice: I need to represent both on the graph
def getRhoAn(z, t , rhoCurve=None):
	if(problemtype == 'sound_wave'):
		if rhoCurve is None:
			rhoCurve = getRhoCurve(z, t)
		from sound_wave_params import rho00
		return np.add(np.multiply(rhoCurve ,rho00), rho00)

def getPresAn(z, t, presCurve=None):
	if(problemtype == 'sound_wave'):
		if presCurve is None:
			presCurve = getPresCurve(z, t)
		from sound_wave_params import p00
		return presCurve * p00 + p00

def getVelAn(z, t, velCurve=None):
	if(problemtype == 'sound_wave'):
		if velCurve is None:
			velCurve = getVelCurve(z, t)
		from sound_wave_params import v00
		result = np.add(np.multiply(velCurve, getCs0()), v00)	
		return result

	

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
		from sound_wave_params import v00
		cs00 = getCs0() 
		return np.divide(np.subtract(v, v00),cs00)
