
import numpy as np
from constants import gamma
from math import sqrt, atan, atanh, cos, log

p00 = 1.0
v00 = 0.0

#try with homog!
#cs00 = math.sqrt(gamma * p00 / rho00)
#v00 = - cs00 /  5.5
#v00 =  0.5 * cs00
#v00 = - cs00
#v00 =  cs00

#mediumType = "homog"  
mediumType = "inhomog"  #variable density rho00 to test with wave packet

if mediumType == "homog":
	rho00 = 1.0
	#rho00 = 0.3  #second exp of inhom
	cs00 = sqrt(gamma * p00 / rho00)

elif(mediumType=="inhomog"):
	inhomogSubtype = 1	
	#inhomogSubtype = 2	
	from constants import z0, zf
	if inhomogSubtype == 1:	
		rho00 = 1.0
		rho01 = 0.01
		#ze = 0.5*(z0 + zf)   #start of change in rho
		ze = z0 + 0.2*(zf - z0) #first exp of inhom new , more at the beginning
		#ze = z0 + 0.15*(zf - z0) #first exp of inhom new , even more 
		we = 0.4
	else:	
		rho01 = 1.2#second exp of inhom
		rho00 = 0.3
		#inverse
		#rho01 = 1.0#second exp of inhom
		#rho00 = 0.01
		#ze = z0 + 0.7*(zf - z0) #second exp of inhom
		#ze = z0 + 0.2*(zf - z0) #first exp of inhom new , more at the beginning
		ze = 0.5*(zf + z0) #middle
		we = 0.5
	#we = 0.4
	#we = 0.2
	def densFunc(z):
		return rho00 + 0.5 * (rho01-rho00) * (1 + np.tanh((z-ze)/we))
	#desympy
	sqrtDensPowMinusOneDer = lambda z:(-(-0.5*rho00 + 0.5*rho01)*(-np.tanh((z - ze)/we)**2 + 1)/(2*we*(rho00 + (-0.5*rho00 + 0.5*rho01)*(np.tanh((z - ze)/we) + 1))**(3/2)))

	#de mathematica
	#but the value evaluated in 3.1 for example is complex!!!
	sqrtDensIntMathematica = lambda z: sqrt(rho01)*we*np.arctanh((sqrt(2)*sqrt(rho01 + rho00 + (rho01 - rho00)* np.tanh((z - ze)/we)))/sqrt(rho01)) -  sqrt(rho00)*we*np.arctanh((sqrt(2)*sqrt(rho01 + rho00 + (rho01 - rho00)*np.tanh((z - ze)/we)))/sqrt(rho00))




	#calculate numerically!
	def sqrtDensIntNumeric(z):
		from constants import z0
		intFunc = lambda z: np.sqrt(densFunc(z))
		from scipy import integrate
		def getIntOneValue(zval):
			int1 = integrate.quad(intFunc, z0, zval)	
			return int1[0]
		if(hasattr(z,'__len__')):
			res = np.zeros(len(z))
			for i in range(len(z)):
				res[i] = getIntOneValue(z[i])	
			return res
		return getIntOneValue(z)

	sqrtDensInt = sqrtDensIntNumeric			

	csderAnal = lambda z: np.sqrt(gamma * p00) * sqrtDensPowMinusOneDer(z)
	def getXIntF(z):
		return (1.0 / np.sqrt(gamma * p00)) * sqrtDensInt(z)
	

	def getDensVarLength(z):
		d = densFunc(z)
		absGrad = abs(np.gradient(d, z[1]-z[0]))
		indNotZero = (absGrad!=0)
		return np.min(abs(d)[indNotZero] / absGrad[indNotZero])


	def cs00(z):
		return np.sqrt(gamma * p00 / densFunc(z))


	#HOMOG
	#try homgenous distribution!
#	densFunc = lambda z: rho00 * np.ones(z.shape)
#	sqrtDensPowMinusOneDer 	= lambda z: np.zeros(z.shape)
#	sqrtDensInt = lambda z: sqrt(rho00) * z



