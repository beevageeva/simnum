import numpy as np
from soundwave_medium_params import mediumType
from initcond_soundwave import getWFunction, fromCurvesToVals
from soundwave_boundary_conditions import getPeriodicXArray

def getVals(z, curves):
	from soundwave_medium_params import p00, v00
	if mediumType == "homog":
		from soundwave_medium_params import rho00
		rhoIni = rho00
	elif mediumType == "inhomog":
		from soundwave_medium_params import densFunc
		rhoIni = densFunc(z)
	presCurve = curves['pres']
	rhoCurve = curves['rho']
	velCurve = curves['vel']
	return fromCurvesToVals(p00, rhoIni, v00, presCurve, rhoCurve, velCurve)


if mediumType == "homog":
	def getWFunctionVals(functiontype, csSign, z, t):
		from soundwave_medium_params import cs00		
		newZ = getPeriodicXArray(z - csSign * cs00 * t)
		return getWFunction(functiontype)(newZ)


elif mediumType == "inhomog":
	from common import getZArray	
	newZ = getZArray()

	print("creating newZ for analytic sol and  " )

	
	def getZ0Index(zval):
		delta = 0.0001
		for index in range(newZ.shape[0]):
			if (abs(newZ[index] - zval)<delta):
				return index


	def getWFunctionVals(functiontype, csSign, z, time):
		if(functiontype != "wavepacket"):
			print("analytical method not implemented")
			import sys
			sys.exit(0)
		from constants import z0, zf
		indexZ0 = np.zeros(newZ.shape)		
		for index in range(newZ.shape[0]):
			indexZ0[index] = getZ0Index(z[index])
			#print(indexZ0[index])
		nans, x= np.isnan(indexZ0), lambda z: z.nonzero()[0]
		indexZ0[nans]= np.interp(x(nans), x(~nans), indexZ0[~nans])
		indexZ0 = indexZ0.astype(int)
		np.set_printoptions(threshold='nan')	
		print(indexZ0)
		firstZ = z[indexZ0]
		from soundwave_medium_params import cs00
		from soundwave_perturbation_params import A	
		csZ0 = cs00(firstZ)
		csZt = cs00(z)
		from sound_wave_packet_params import getSoundWaveGaussFunction, zc, W, k0
		ampIni = getSoundWaveGaussFunction(zc, W)(firstZ)
		#print("AMP INI")
		#print(ampIni)
		curve = np.zeros(firstZ.shape)
		from math import pi
		for index in range(firstZ.shape[0]):
			curve[index] =  ampIni[index] * csZ0[index] ** (0.5) * csZt[index] ** (-0.5) *  A * np.cos( 2 * pi * k0 * (csZ0[index]/csZt[index]) / (zf - z0) * ( z[index] - csZt[index] * time) - 2 * pi * k0 * z0 / (zf - z0))
	
		print("curve")

		#print(curve)
		print("MAX CURVE = %e , IND MAX CURVE =%d, min CURVE = %e , ind MIN curve = %d" % (np.max(curve), np.argmax(curve), np.min(curve), np.argmin(curve) ) )
		
		

		return curve
	

	def updateNewZ(modelObj, dt, csNumerical):
		for index in range(newZ.shape[0]):
			newZ[index] = modelObj.getNewPoint(newZ[index],dt, False)





def getCurves(z, t):
	from soundwave_perturbation_params import perturbationType
	if perturbationType == "superposition":
		from soundwave_perturbation_params import  init_functions_generation
		velPert = np.zeros(z.shape)
		presRhoPert = np.zeros(z.shape)
		for fdef in init_functions_generation:
			wFuncVals =  getWFunctionVals(fdef['functiontype'],fdef['csSign'] , z, t)
			velPert += fdef['csSign'] * fdef['A'] * wFuncVals
			presRhoPert += fdef['A'] * wFuncVals
		return {'pres':presRhoPert, 'rho':presRhoPert, 'vel': velPert}
	elif perturbationType == "one":
		from soundwave_perturbation_params import  A, functiontype
		velPresRhoPert = A * getWFunctionVals(functiontype ,1 , z, t)    
		return {'pres':velPresRhoPert, 'rho':velPresRhoPert, 'vel': velPresRhoPert}


	


