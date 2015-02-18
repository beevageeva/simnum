import numpy as np
import sys,math
from constants import gamma,zf, z0

from sound_wave_params import p00, rho00, v00, A, functiontype, periodicType, mediumType

def getWFunction():
	if functiontype == "sine":
		from sound_wave import getSoundWaveSineFunction
		from sound_wave_sine_params import phi0
		wl = zf - z0
		phi = phi0 - 2.0 * math.pi * z0 / wl
		return getSoundWaveSineFunction(wl, phi)
	elif functiontype == "gauss":
		from sound_wave_gauss_params import R
		from sound_wave import getSoundWaveGaussFunction
		return  getGaussFunction(R) 	
	elif functiontype == "bessel":
		from sound_wave import getSoundWaveBesselFunction
		return getSoundWaveBesselFunction()
	elif functiontype == "wavepacket":
		from sound_wave import getSoundWavePacketFunction
		from sound_wave_packet_params import k0, zc, W
		return getSoundWavePacketFunction(k0, zc, W)
	else:	
		print("functiontype %s not implemented" % functiontype)
		sys.exit(0)




if mediumType == "homog":
	from sound_wave import SoundWave
	
	def createWave(csSign , A):
		return SoundWave(A, cs00 *csSign + v00, getWFunction() )	
	
	wavePresRho = None
	waveVel = None
	cs00 = math.sqrt(gamma * p00 / rho00)
	
	from sound_wave_params import init_functions_generation 
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
	
	def getCs00():
		from sound_wave_params import p00, rho00
		cs = math.sqrt(gamma *  p00 / rho00)
		return cs
	def getVelFFTAn(k):
		if functiontype == "wavepacket":
			from sound_wave_params import A, p00, rho00
			cs00 = math.sqrt(gamma * p00 / rho00)
			from sound_wave import getSoundWavePacketFFTAnalyticalAbs as wFFTAn
			from sound_wave_packet_params import k0, zc, W
			return A * cs00 * wFFTAn(k0, zc, W)(k)
		return None
	
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
	
	
	
	
	def getPresCurveNumeric(p):
		from sound_wave_params import p00
		return np.divide(np.subtract(p,p00),p00)
	
	def getVelCurveNumeric(v):
		from sound_wave_params import v00
		cs00 = getCs0() 
		return np.divide(np.subtract(v, v00),cs00)
	
	#analitycal values
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
		result = np.add(np.multiply(velCurve, getCs00()), v00)	
		return result
	
	#analitycal values end

else:

	def getCs00(z):
		from sound_wave_params import densFunc
		rhoIni =  densFunc(z)
		return np.sqrt(np.divide(gamma * p00,rhoIni))


	def getInitialPresRhoVel(z):
		from sound_wave_params import densFunc
		rhoIni =  densFunc(z)
		cs00 = np.sqrt(np.divide(gamma * p00,rhoIni))
		w = getWFunction()
		f = w(z)
		return {'pres': p00 + gamma * p00 * A* f  , 'rho': rhoIni + rho00 *A* f , 'vel': v00 + cs00 * A* f }

	def getInitialFunctionMaxMinZIndex(z):
		w = getWFunction()(z)
		return [np.argmin(w), np.argmax(w)]


	if functiontype == "wavepacket":
		def wAnal(z, t , cs):	
			omega0 = np.mean(k0 * cs)
			t2 = np.subtract(z,zc+ omega0*t) ** 2
			return np.multiply(np.exp(-np.divide(t2, W**2)), np.cos(k * (z-z0 - t * cs)))
							
		#analitycal values
		def getRhoCurve(z, t):
			res =  wavePresRho.getWaveShape(z, t)	
			return res
					
		def getVelCurve(z, t):
			return  wAnal(z,t,getCs00(z))
					
		
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
			result = np.add(np.multiply(velCurve, getCs00(z)), v00)	
			return result
		
		#analitycal values end


def getRhoCurveNumeric(rho,z):
	from sound_wave_params import mediumType, rho00
	if(mediumType == "homog"):	
		rhoIni = rho00
	else:
		from sound_wave_params import densFunc
		rhoIni=densFunc(z)
	return np.divide(np.subtract(rho,rhoIni),rho00)




if periodicType == "repeat":

	def  lrBoundaryConditionsPresRho(array, skip=0):
		n = array.shape[0] - 1
		array = np.insert(array, 0,  array[n-skip])
		array = np.insert(array, n+2,  array[1+skip])
		return array

	lrBoundaryConditionsVel = lrBoundaryConditionsPresRho

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
		




		

