import numpy as np
import sys,math
from constants import gamma,zf, z0

from sound_wave_params import p00, rho00, v00, A, functiontype, periodicType, mediumType

def getWFunction():
	if functiontype == "sine":
		from sound_wave_sine_params import phi, k0, getSoundWaveFunction
		return getSoundWaveFunction(k0, phi)
	elif functiontype == "gauss":
		from sound_wave_gauss_params import R, getSoundWaveFunction
		return  getFunction(R) 	
	elif functiontype == "wavepacket":
		from sound_wave_packet_params import k0, zc, W, getSoundWaveFunction
		return getSoundWaveFunction(k0, zc, W)
	elif functiontype == "defined":
		from sound_wave_defined_params import getSoundWaveFunction
		return getSoundWaveBesselFunction()
	else:	
		print("functiontype %s not implemented" % functiontype)
		sys.exit(0)


if periodicType == "repeat":

	def  lrBoundaryConditionsPresRho(array, skip=0):
		n = array.shape[0] - 1
		array = np.insert(array, 0,  array[n-skip])
		array = np.insert(array, n+2,  array[1+skip])
		return array

	lrBoundaryConditionsVel = lrBoundaryConditionsPresRho

	def getPeriodicX(xval, a=z0, b=zf):
		p = float(b - a)
		k = int((xval-a)/p)
		res = xval - k * p
		if(res < a):
			res+=p
		if(res > b):
			res-=p
		return res

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
		

	def getPeriodicX2(xval, a=z0, b=zf):
		p = float(b - a)
		k = int((xval-a)/p)
		if k%2==1:
			res = 2.0*b - xval + k*p
		else:
			res = xval - k * p
		if(res < a):
			res+=p
		if(res > b):
			res-=p
		return res


	def getPeriodicX(xval, a=z0, b=zf):
		if xval > b:
			return 2.0*b - xval
		if xval < a:
			return 2.0*a -xval
		return xval


#TODO eliminate common.py
def getPeriodicXArray(xarray, a=z0, b=zf):
	res = []	
	for xval in xarray:
		res.append(getPeriodicX(xval, a, b))
	return np.array(res)
		



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
			from sound_wave_packet_params import k0, zc, W, getSoundWaveFFTAnalytical
			return A * cs00 * getSoundWaveFFTAnalytical(k0, zc, W)(k)
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
		return {'pres': p00 + gamma * p00 *  wPresRho0  , 'rho': rho00 + rho00 *  wPresRho0 , 'vel': v00 + cs00 *  wVel0 } 
	
	
	
	
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
	from sound_wave_params import densFunc, sqrtDensPowMinusOneDer, sqrtDensInt

	def getCs00(z):
		rhoIni =  densFunc(z)
		return np.sqrt(np.divide(gamma * p00,rhoIni))


	def getInitialPresRhoVel(z):
		#repeated def of cs00 in order not to calculate twice rhoIni
		rhoIni =  densFunc(z)
		cs00 =  np.sqrt(np.divide(gamma * p00,rhoIni))
		w = getWFunction()
		f = w(z)
		return {'pres': p00 + gamma * p00 * A* f  , 'rho': rhoIni + rho00 *A* f , 'vel': v00 + cs00 * A* f }

	def getInitialFunctionMaxMinZIndex(z):
		w = getWFunction()(z)
		return [np.argmin(w), np.argmax(w)]

	csderAnal = lambda z: np.sqrt(gamma * p00) * sqrtDensPowMinusOneDer(z)
	#F in eq dx / dt = cs(x) after integrating: F(x) - t = F(x0) F(x) = Int(1/cs(x))
	#TODO boundary conditions	
	def getXIntF(z):
		return (1.0 / np.sqrt(gamma * p00)) * sqrtDensInt(z)
#	#and solve implicit F(x) = v, for monoton F
#	def getXImpl(z, F,v):
##		n = len(z) - 1
##		sign = 1
##		k1 = 0
##		k2 = n
##		if(z[0]>z[n]):
##			sign = -1
##		kmid = int(0.5 * (k1 + k2))
##		if(sign * F(z[kmid])< sign * v):
##			up = 1
##			k1 = kmid
##		else:
##			up = -1
##			k2 = kmid  
##		while (up * F(z[kmid]) < up * v ):
##			if up:	
##				k1 = kmid
##			else:
##				k2 = kmid
##			kmid = int(0.5 * (k1 + k2))
#			indexZ = np.searchsorted(F(z), v)
#			#interpolation
#			return 0.5 * (z[indexZ] + z[indexZ + 1])

	#from eq F(x0) = F(x) -t  F = getXIntF
	#not used in wAnal as I won't calculate FValues every time
	def getX0Index(z, t, zval):
		from cache import getValue, putValue
		FValues = getValue("FValues")
		#print(FValues)
		if(FValues is None):
			#from common import getPeriodicXArray, getPeriodicX
			print("FValues not cached..")
			FValues = getPeriodicXArray(getXIntF(z))
			print("putting")
			print(FValues)
			putValue("FValues", FValues)

		#TODO boundary conditions eleiminate common!
		#from common import getPeriodicXArray, getPeriodicX
		znp =  getXIntF(zval) - t #ADDED z0 for per cond not working
		#znp = z0 + getXIntF(zval) - t
		val = getPeriodicX(znp)
		print("not per = %e, per = %e" % (znp, val))

		#TODO indices
		#return np.searchsorted(FValues, val)
		indices = []	
		delta = 0.001
		for index in range(FValues.shape[0]):
			if(abs(FValues[index] - val)<delta):
				indices.append(index)
		print("indices")
		print(indices)
		return indices[len(indices)/2]
		


	if functiontype == "wavepacket":
		from sound_wave_packet_params import k0, zc, W, getSoundWaveGaussFunction
		def k0Func(z):
			from math import pi
			return k0 *  2.0 * pi /  (zf - z0)
		def a0Func(z):
			return getSoundWaveGaussFunction(zc, W)(z)
		def getphi0():
			from math import pi
			from constants import z0, zf
			return - 2 * pi * k0 * z0 / (zf - z0)
	elif functiontype == "sine":
		from sound_wave_sine_params import phi, k0, getSoundWaveFunction
		def k0Func(z):
			from math import pi
			return 2 * pi * k0/ (zf - z0)
		def getphi0():
			return phi
		def a0Func(z):
			return np.ones(z.shape)

	def kAnal2(z, t):
		from cache import getValue, putValue
		cs = getCs00(z)
		res = np.zeros(len(z))

	
		#with searchval
		#from common import getPeriodicXArray, getPeriodicX
		FValues = getValue("FValues")

		#print(FValues)
		if(FValues is None):
			print("FValues not cached..")
			FValues = getPeriodicXArray(getXIntF(z))
			print("putting")
			print(FValues)
			putValue("FValues", FValues)
		#endsearchval

		#TODO use numpy operations on whole array
		for index in range(len(z)):
			zval = z[index]
	
			#with searchVal
			FSearchVal = getPeriodicX(getXIntF(zval) -t)
			#x0Index = np.searchsorted(FValues, FSearchVal)
			#TODO indices
			indices = []	
			delta = 0.01
			for index in range(FValues.shape[0]):
				if(abs(FValues[index] - FSearchVal)<delta):
					indices.append(index)
			#print("indices")
			#print(indices)
			x0Index = indices[len(indices)/2]
			#end searchval

			#HOMOG 
#			from common import getZIndex
#			x0Index = getZIndex(zval - t * cs[index])

			#print("x0Index = %d" % x0Index)
			#k0x0 = k0Func(z[x0Index])
			k0x0 = k0
			#csPrime = csderAnal(zval)
			csPrime = csderAnal(z[x0Index])
			#endsearchval
			#with derivative BUT the following is not equiv
#			k0x0 =  k0
#			csPrime = csderAnal(zval) * cs[index]
			#endderivative

			res[index] = k0x0 * np.exp(-csPrime * t)
		return res



	def kAnal(z0, t, z):
		from cache import getValue, putValue
		FValues = getValue("FValues")
		#print(FValues)
		if(FValues is None):
			#from common import getPeriodicXArray, getPeriodicX
			print("FValues not cached..")
			FValues = getPeriodicXArray(getXIntF(z))
			print("putting")
			print(FValues)
			putValue("FValues", FValues)
		cs = getCs00(z)
		csder = csderAnal(z)
		def intT(t):
			indices = []	
			delta = 0.01
			FSearchVal = t + F(z0)
			for index in range(FValues.shape[0]):
				if(abs(FValues[index] - FSearchVal)<delta):
					indices.append(index)
			#print("indices")
			#print(indices)
			x0Index = indices[len(indices)/2]
			return cs[x0Index] + csder[x0Index]
		from scipy import integrate
		int1 = integrate.quad(intT, 0, t)
		return np.exp(-int1[0])	
		

		



	def wAnal(z, t):
		from cache import getValue, putValue
		if getValue("timeWAnal") == t:
			return getValue("wAnal") 
		cs = getCs00(z)
		#print("cs is")
		#print(cs)
		res = np.zeros(len(z))
		a0 = a0Func(z)
		np.set_printoptions(threshold='nan')
		#print("a0 = ")
		#print(a0)
		#print("MAXINDEX is %d " % np.argmax(a0))
		phi0 = getphi0()
		from common import getPeriodicXArray, getPeriodicX
		FValues = getValue("FValues")
		#print(FValues)
		if(FValues is None):
			print("FValues not cached..")
			FValues = getPeriodicXArray(getXIntF(z))
			print("putting")
			print(FValues)
			putValue("FValues", FValues)
		#TODO use numpy operations on whole array
		for index in range(len(z)):
			zval = z[index]
			FSearchVal = getPeriodicX(getXIntF(zval) -t)
			x0Index = np.searchsorted(FValues, FSearchVal)

			#HOMOG 
#			from common import getZIndex
#			x0Index = getZIndex(zval - t * cs[index])

			if x0Index >= z.shape:
				x0Index -=1
			#print("x0Index = %d" % x0Index)
			k0x0 = k0Func(z[x0Index])
			omega0  = cs[x0Index] * k0x0
			csPrime = csderAnal(zval)
			k = k0x0 * np.exp(-csPrime * t)
			a = cs[index] * a0[x0Index] / cs[x0Index] * np.exp(-0.5 * csPrime * t)
			#print("index = %d, FSearchVal = %e, indexX0=%d, a = %e, k = %e, omega = %e, a01 = %e, a02 = %e" % (index, FSearchVal, x0Index, a, k, omega0, a0[x0Index], cs[index] * a0[x0Index] / cs[x0Index]))
			res[index] = a * np.cos(k * zval + omega0 * t +phi0 )
		from sound_wave_params import A
		finalVals =  A * res
		putValue("wAnal", finalVals)
		putValue("timeWAnal", t)
		return finalVals
		
	def wAnal2(z, t):
		if functiontype == "wavepacket":
			from math import pi
			from sound_wave_packet_params import k0, zc, W, getSoundWaveGaussFunction
			cs = getCs00(z)
			omega = np.mean(cs) * 2.0 * pi * k0 / (zf - z0)
			from sound_wave_params import A
			return A*np.cos(2.0 * pi * k0 * (z-z0)/ (zf - z0) - omega * t ) * getSoundWaveGaussFunction(zc , W)(z - cs * t)
		

						
	#analitycal values TO KEEP them like in homog case 
	getRhoCurve = wAnal2
	getVelCurve = wAnal2
	
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




