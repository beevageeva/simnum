import numpy as np
from constants import gamma
from math import sqrt, atan, atanh, cos, log

rho00 = 1.0
#rho00 = 0.3  #second exp of inhom

#mediumType = "homog"
mediumType = "inhomog"  #variable density rho00 to test with wave packet
if(mediumType=="inhomog"):
	inhomogSubtype = 1	
	#inhomogSubtype = 2	
	from constants import z0, zf
	if inhomogSubtype == 1:	
		rho01 = 0.01
		#ze = 0.5*(z0 + zf)   #start of change in rho
		ze = z0 + 0.2*(zf - z0) #first exp of inhom new , more at the beginning
	else:	
		rho01 = 1.2#second exp of inhom
		ze = z0 + 0.7*(zf - z0) #second exp of inhom
	#we = 0.4
	we = 0.2
	densFunc = lambda z: rho00 + 0.5 * (rho01-rho00) * (1 + np.tanh((z-ze)/we))
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

	#try homgenous distribution!
	densFunc = lambda z: rho00 * np.ones(z.shape)
	sqrtDensPowMinusOneDer 	= lambda z: np.zeros(z.shape)
	sqrtDensInt = lambda z: sqrt(rho00) * z

#functiontype = 'sine'
#functiontype = 'gauss'
functiontype = 'wavepacket'
#functiontype = 'defined'

periodicType = "repeat" 
#periodicType = "refl" 

p00 = 1.0


v00 = 0.0

#cs00 = math.sqrt(gamma * p00 / rho00)
#v00 = - cs00 /  5.5
#v00 =  0.5 * cs00
#v00 = - cs00
#v00 =  cs00

A = 3.0 * 10.0 ** (-4)
#A = 5.0 * 10.0 ** (-2)
#A = 10.0 ** (-2) #this will work fine as well with fg scheme

if(mediumType == "homog"):
	#init_functions_generation = [{'csSign':1, 'A': A}, {'csSign':-1, 'A': 0.5*A}] #SUPERPOSITION wave travelling right with amp A and left with amp 0.5 * A
	#init_functions_generation = [{'csSign':-1, 'A': A}, {'csSign':1, 'A': 0.5*A}] #SUPERPOSITION wave travelling left with amp A and right with amp 0.5 * A
	#init_functions_generation = [{'csSign':-1, 'A': A}, {'csSign':1, 'A': A}] #SUPERPOSITION wave travelling left with amp A and right with amp  A
	#init_functions_generation = [{'csSign':-1, 'A': A}] #travelling left
	init_functions_generation = [{'csSign':1, 'A': A}] #travelling right



#The following are used for the notifier(taken from notifier_params as they were wavesound specific)
plotPresCurve = False
#plotRhoCurve = False
plotRhoCurve = True
plotVelCurve = False
#plotPresAn = False
#plotRhoAn = False
#plotVelAn = False
plotPresAn = True
plotRhoAn = True
plotVelAn = True
#in inhomogeneous medium maximum does not travel at cs speed(phase velocity) , but  at group velocity, see initcond_sounwave the functio
#to get max for ini pres and rho is only defined for homog medium
markPoints = True
plotVelFFT = True
#plotVelFFT = False
#plotVelFFTAnal=False
plotVelFFTAnal=True
plotPresFFT = True
#plotPresFFT = False

#if(mediumType == "inhomog"):
if(mediumType == "inhomog" ):
	plotPresCurve = False
	plotVelCurve = False
	plotVelFFTAnal = False


#adding packet analytical solution for inhomog type
#if periodicType == "refl" or mediumType == "inhomog":
if periodicType == "refl" or mediumType == "inhomog" and functiontype!="wavepacket":
	#analytical function does not make sense
	plotPresAn = False
	plotRhoAn = False
	plotVelAn = False




















