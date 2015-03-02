import numpy as np
import sys,math
from constants import gamma,zf, z0


from soundwave_perturbation_params import perturbationType
from soundwave_medium_params import mediumType

def getWFunction(functiontype):
	if functiontype == "sine":
		from sound_wave_sine_params import phi, k0, getSoundWaveFunction
		return getSoundWaveFunction(k0, phi)
	elif functiontype == "gauss":
		from sound_wave_gauss_params import R, getSoundWaveFunction
		return  getSoundWaveFunction(R) 	
	elif functiontype == "wavepacket":
		from sound_wave_packet_params import k0, zc, W, getSoundWaveFunction
		return getSoundWaveFunction(k0, zc, W)
	elif functiontype == "defined":
		from sound_wave_defined_params import getSoundWaveFunction
		return getSoundWaveBesselFunction()
	else:	
		print("functiontype %s not implemented" % functiontype)
		sys.exit(0)

def fromCurvesToVals(pEq, rhoEq, vEq, presPert, rhoPert, velPert):
	csEq = np.sqrt(gamma * pEq / rhoEq)
	return {'pres': pEq + gamma * pEq * presPert  , 'rho': rhoEq + rhoEq * rhoPert , 'vel': vEq + csEq * velPert }


def fromValsToCurvePres(p):
	from soundwave_medium_params import p00
	np.set_printoptions(threshold='nan')
	print("numeric pres curve")	
	curve = np.divide(np.subtract(p,p00),p00)
	print("MAX CURVE = %e , IND MAX CURVE =%d, min CURVE = %e , ind MIN curve = %d" % (np.max(curve), np.argmax(curve), np.min(curve), np.argmin(curve) ) )
	#print(np.divide(np.subtract(p,p00),p00))	
	return np.divide(np.subtract(p,p00),p00)

def getInitialPresRhoVel(z):
	from soundwave_medium_params import p00, v00
	if mediumType == "homog":
		from soundwave_medium_params import rho00
		rhoIni = rho00
	elif mediumType == "inhomog":
		from soundwave_medium_params import densFunc
		rhoIni = densFunc(z)
	if perturbationType == "superposition":
		from soundwave_perturbation_params import  init_functions_generation
		velPert = np.zeros(z.shape)
		presRhoPert = np.zeros(z.shape)
		for fdef in init_functions_generation:
			wFuncVals = getWFunction(fdef['functiontype'])(z)	
			velPert += fdef['csSign'] * fdef['A'] * wFuncVals
			presRhoPert += fdef['A'] * wFuncVals
		return fromCurvesToVals(p00, rhoIni, v00, presRhoPert, presRhoPert, velPert)
	elif perturbationType == "one":
		from soundwave_perturbation_params import  A, functiontype
		velPresRhoPert = A * getWFunction(functiontype)(z)
		return fromCurvesToVals(p00, rhoIni, v00, velPresRhoPert, velPresRhoPert, velPresRhoPert)
	

def fromValsToCurveRho(rho,z):
	if(mediumType == "inhomog"):
		from soundwave_medium_params import  densFunc
		rhoIni=densFunc(z)
	else:
		from soundwave_medium_params import  rho00 as rhoIni
	return np.divide(np.subtract(rho,rhoIni),rhoIni)

	def fromValsToCurveVel(v, z):
		from soundwave_medium_params import v00,cs00
		if mediumType == "homog":
			return np.divide(np.subtract(v, v00),cs00)
		else:
			return np.divide(np.subtract(v, v00),cs00(z))

if mediumType == "homog":
	

	def getVelFFTAn(k):
		from soundwave_medium_params import  p00, rho00
		from soundwave_perturbation_params import perturbationType
		if(perturbationType == "one"):
			from soundwave_perturbation_params import A,functiontype
			if functiontype == "wavepacket":
				cs00 = math.sqrt(gamma * p00 / rho00)
				from sound_wave_packet_params import k0, zc, W, getSoundWaveFFTAnalytical
				return A * cs00 * getSoundWaveFFTAnalytical(k0, zc, W)(k)
		print("velFFTAn not implemented")	
		return None
	
			
	
	
	




