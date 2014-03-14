import numpy as np
import sys,math
from constants import gamma

from sound_wave import SoundWave, SoundWaveSine

def createWave(csSign , A):
	if functiontype == "sine":
		from sound_wave_sine_params import wl, phi
		return SoundWaveSine(A, cs00 *csSign + v00,  wl, phi)
	elif functiontype == "defined":
		from sound_wave_defined_params import w
		return SoundWave(A, cs00 *csSign + v00, w)	
	else:	
		print("functiontype %s not implemented" % functiontype)
		sys.exit(0)

wavePresRho = None
waveVel = None
from constants import functiontype
from sound_wave_params import p00, rho00, v00, A, init_functions_generation
cs00 = math.sqrt(gamma * p00 / rho00)

if len(init_functions_generation)==1:
	wavePresRho= createWave(init_functions_generation[0]['csSign'] , init_functions_generation[0]['A'])
	waveVel= createWave(init_functions_generation[0]['csSign'] , init_functions_generation[0]['csSign'] *  init_functions_generation[0]['A'])
else:
	waveArrayPresRho = []
	waveArrayVel = []
	for fdef in init_functions_generation:
		waveArrayPresRho.append(createWave(fdef['csSign'], fdef['A']))
		waveArrayVel.append(createWave(fdef['csSign'], fdef['csSign'] * fdef['A']))
	from sound_wave import SuperpositionSoundWave
	wavePresRho = SuperpositionSoundWave(waveArrayPresRho)	
	waveVel = SuperpositionSoundWave(waveArrayVel)	
	


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

def getInitialFunctionMaxZ(z):
	wVel0= waveVel.getWaveShape(z)
	wPresRho0= wavePresRho.getWaveShape(z)
	markPointPresRho = z[np.argmax(wPresRho0)]
	markPointVel = z[np.argmax(wVel0)]
	return {'pres': markPointPresRho, 'vel': markPointVel, 'rho': markPointPresRho}
		


def getInitialPresRhoVel(z):
	from sound_wave_params import A, p00, rho00, v00
	cs00 = math.sqrt(gamma * p00 / rho00)
	#wPresRho0 = wave.getInitialShape(z)
	#wVel0 = wave.getInitialShape(z, "vel")
	#print(wPresRho0)
	#print(wVel0)
	#return {'pres': p00 + gamma * p00 * wPresRho0  , 'rho': rho00 + rho00 *  wPresRho0 , 'vel': v00 + cs00 *  wVel0 } 
	wVel0= waveVel.getWaveShape(z)
	wPresRho0= wavePresRho.getWaveShape(z)
	return {'pres': p00 + gamma * p00 * wPresRho0  , 'rho': rho00 + rho00 *  wPresRho0 , 'vel': v00 + cs00 *  wVel0 } 

			

def getRhoCurve(z, t):
	res =  wavePresRho.getWaveShape(z, t)	
	return res
			
def getVelCurve(z, t):
	return  waveVel.getWaveShape(z, t)
			

def getPresCurve(z, t):
	return gamma * getRhoCurve(z,t)	

#we can send curves calculated previously as parameters in order not to calculate them twice: I need to represent both on the graph
def getRhoAn(z, t , rhoCurve=None):
	if rhoCurve is None:
		rhoCurve = getRhoCurve(z, t)
	from sound_wave_params import rho00
	return np.add(np.multiply(rhoCurve ,rho00), rho00)

def getPresAn(z, t, presCurve=None):
	if presCurve is None:
		presCurve = getPresCurve(z, t)
	from sound_wave_params import p00
	return presCurve * p00 + p00

def getVelAn(z, t, velCurve=None):
	if velCurve is None:
		velCurve = getVelCurve(z, t)
	from sound_wave_params import v00
	result = np.add(np.multiply(velCurve, getCs0()), v00)	
	return result

	

def getRhoCurveNumeric(rho):
	from sound_wave_params import rho00
	return np.divide(np.subtract(rho,rho00),rho00)


def getPresCurveNumeric(p):
	from sound_wave_params import p00
	return np.divide(np.subtract(p,p00),p00)

def getVelCurveNumeric(v):
	from sound_wave_params import v00
	cs00 = getCs0() 
	return np.divide(np.subtract(v, v00),cs00)
