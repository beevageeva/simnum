import numpy as np
from soundwave_medium_params import mediumType
from initcond_soundwave import getWFunction, fromCurvesToVals
from soundwave_boundary_conditions import getPeriodicXArray

def getValsCurvesNotShifted(z, curves):
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

	getVals = getValsCurvesNotShifted

elif mediumType == "inhomog":
	from common import getZArray	
	firstZ = getZArray()
	newZ = getZArray()

	shiftNow = False
	print("creating newZ for analytic sol and  shiftNow is %s" % (str(shiftNow)))

	def getWFunctionVals(functiontype, csSign, z, t):
		if(functiontype != "wavepacket"):
			print("analytical method not implemented")
			import sys
			sys.exit(0)
		
		csZ0 = cs00(firstZ)
		csZt = cs00(newZ)
		
		ampIni = getSoundWaveGaussFunction(zc, W)(firstZ)
		curve = np.full(firstZ.shape, np.nan)
		for index in range(firstZ.shape[0]):
			if shiftNow:
				cIndex = getZIndex(newZ[index])
				if(np.isnan(curve[cIndex])):
					curve[cIndex] = 0
				curve[cIndex] += ampIni[index] * csZ0[index] ** (0.5) * csZt[index] ** (-0.5) *  A * np.cos( 2 * pi * k0 * (csZ0[index]/csZt[index]) / (zf - z0) * ( newZ[index] - csZt[index] * time) - 2 * pi * k0 * z0 / (zf - z0))
			else:
				curve[index] = ampIni[index] * csZ0[index] ** (0.5) * csZt[index] ** (-0.5) *  A * np.cos( 2 * pi * k0 * (csZ0[index]/csZt[index]) / (zf - z0) * ( newZ[index] - csZt[index] * time) - 2 * pi * k0 * z0 / (zf - z0))
	
		if shiftNow:
			#interpolate nan values!
			nans, x= np.isnan(curve), lambda z: z.nonzero()[0]
			curve[nans]= np.interp(x(nans), x(~nans), curve[~nans])
		return {'pres': curve,'rho': curve,'vel': curve}

	if shiftNow:
		getVals = getValsCurvesNotShifted
	else:	
		def getVals(z, curves):
			rhoIni = densFunc(newZ) #= gamma * p00  cs(z(t))** (-2) shift afterwards
			velAn =   v00 + csZt * curves['vel'] 
			presAn =  p00 + p00 * gamma * curves['pres']
			rhoAn =   rhoIni + rhoIni * curves['rho']
			return {'pres' : presAn, 'rho' : rhoAn, 'vel': velAn}
	

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


	


