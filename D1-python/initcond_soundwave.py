import numpy as np
import sys,math
from constants import gamma,zf, z0

from sound_wave_params import p00, rho00, v00, A, init_functions_generation, functiontype, periodicType
from sound_wave import SoundWave, SoundWaveSine

def createWave(csSign , A):
	if functiontype == "sine":
		from sound_wave_sine_params import phi0
		wl = zf - z0
		phi = phi0 - 2.0 * math.pi * z0 / wl
		return SoundWaveSine(A, cs00 *csSign + v00,  wl, phi)
	elif functiontype == "defined":
		from sound_wave_defined_params import w
		return SoundWave(A, cs00 *csSign + v00, w)	
	else:	
		print("functiontype %s not implemented" % functiontype)
		sys.exit(0)

wavePresRho = None
waveVel = None
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

print("PresRhoWave")
wavePresRho.printString()		
print("VelWave")
waveVel.printString()		


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



if periodicType == "repeat":

	def  lrBoundaryConditionsPresRho(array, skip=0):
		n = len(array) - 1
		array.insert(0, array[n - skip])
		array.append(array[1 + skip])

	lrBoundaryConditionsVel = lrBoundaryConditionsPresRho

elif periodicType == "diff":

		

	def lrBoundaryConditionsPresRho(array, skip=0):
		array.insert(0, 2 * array[0] - array[1])
		array.append(2 * array[-1] - array[-2])
#		#I already know
#		#if(len(array)<2):
#		#	return
#		if(skip == 0):
#			array.insert(0, 2 * array[0] - array[1])
#			array.append(2 * array[-1] - array[-2])
#			return 
#		#TODO no more polyfit from np
#		#print("array 1")
#		#print(array)
#		degree = 1+skip
#		xvalues = range(1,degree+2)
#		#print("xvalues = ")
#		#print(xvalues)
#		#print("degree = %d" % degree)
#		yvalues = array[0:degree+1]
#		#print("yvalues1 = ")
#		#print(yvalues)
#		c = np.polyfit(xvalues,yvalues, degree)
#		#print("coef1 = ")
#		#print(c)
#		p = np.poly1d(c)
#		#print("first val = %4.3f" % p(0))
#		array.insert(0, p(0))	
#		yvalues = array[-degree-1:]
#		#print("yvalues2 = ")
#		#print(yvalues)
#		c = np.polyfit(xvalues,yvalues, degree)
#		#print("coef2 = ")
#		#print(c)
#		p = np.poly1d(c)
#		#print("last val = %4.3f" % p(degree+2))
#		array.append(p(degree+2))	
#		#print("array 2")
#		#print(array)

	def lrBoundaryConditionsVel(array, skip=0):
		#I already know
		#if(len(array)<1+skip):
		#	return
		if(skip==0):
			array.insert(0, -array[0])
			array.append(-array[-1]) 
		elif (skip==1):
			array[0] = 0
			array.insert(0, -array[2])
			array[-1] = 0
			array.append(-array[-2]) 




def getVelFFTAn(k):
	if functiontype == "defined":
		from sound_wave_params import A, p00, rho00
		cs00 = math.sqrt(gamma * p00 / rho00)
		from sound_wave_defined_params import wFFTAn
		return A * cs00 * wFFTAn(k)
	return None
		

