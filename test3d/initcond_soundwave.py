import numpy as np
import sys,math
from constants import gamma
from sound_wave_params import mediumType, wType



if(mediumType == "homog"):
	def getRho00(z):
		from sound_wave_params import rho00
		return rho00
	
else:
	def getRho00(z):
		from sound_wave_params import rho01, ze, densFunc, densargfunc, rho00
		from common import getArrayZShape
		return rho00 + 0.5 * (rho01-rho00) * densFunc(densargfunc(z - getArrayZShape(ze[0], ze[1], len(z[0]))))
		

		
if wType == "pot":


	def getGradientFunction(f, z):
		try:
			from sound_wave_params import symDerivW
			return symDerivW(z)
		except ImportError:
			from common import getDz0, getDz1
			dz0 = getDz0()
			dz1 = getDz1()
			return np.gradient(f,dz0, dz1 )
			


	def getAnValuesFromVelPot(z, t):
		#only for homogenous medium, otherwise exists(see sound_wave_params)
		from sound_wave_params import A, p00, rho00, v00, w, velFunc
		from sound_wave_params import k
		from common import derivZ0, derivZ1
		cs00 = math.sqrt(gamma * p00 / rho00) #only for homogenous
		f = w(z)
		omega = k * cs00
		r = np.sqrt(z[0]**2 + z[1]**2)
		gradW = getGradientFunction(f, z)
		v1 = np.real(np.multiply(gradW, np.exp(-1j * omega * t)))
		velPert = A* np.dstack(velFunc(v1, z))
		presPert = rho00 * omega * A * np.real(1j * f * np.exp(-1j * omega * t))
		rhoPert = presPert / cs00 ** 2
		return {'pres': p00 + presPert  , 'rho': rho00 + rhoPert , 'vel': v00 + velPert } 





def getInitialPresRhoVel(z):
	from sound_wave_params import A, p00, rho00, v00, w, wType,timesZArgW
	f =  w(z)
	if  wType == "all":
		from sound_wave_params import velFunc, mediumType
		rhoIni = getRho00(z)
		cs00 = np.sqrt(np.divide(gamma * p00,rhoIni))
		#initial velocity
		v1 = v00 + cs00 * A* f
		vel = np.dstack(velFunc(v1, z))
		return {'pres': p00 + gamma * p00 * A* f  , 'rho': rhoIni + rho00 *A* f , 'vel': vel }
	elif wType == "pot":
		return getAnValuesFromVelPot(z, 0)



def getAnRhoPresVel(z, t):
	from sound_wave_params import mediumType, timesZArgW
	if mediumType !="homog":
		print("analitical solution not implemented FOR inhomogen  medium  set plotAnalitical = False in notifier_params.py")
		sys.exit(0)
	from sound_wave_params import A, p00, rho00, v00, wType, w, timesZArgW
	cs00 = math.sqrt(gamma * p00 / rho00) #only for homogenous
	#TODO chapuza
	if(timesZArgW==1):
		from sound_wave_params import argType
	if(wType == "pot"):
		return getAnValuesFromVelPot(z, t)
	#no analitical solution for wave packet
	
	elif(wType == "all" and timesZArgW == 1 and (argType == "x" or argType == "y") and periodicType == "repeat"):
		from sound_wave_params import velFunc
		from common import getPeriodicXArray2
		from constants import z0_0, zf_0, z0_1, zf_1
		if(argType == "x"):
			newz0 = z[0] - cs00 * t
			z = [getPeriodicXArray2(newz0, z0_0, zf_0), z[1]]
		elif(argType == "y"):
			newz1 = z[1] - cs00 * t
			z = [z[0], getPeriodicXArray2(newz1, z0_1, zf_1)]
		f = w(z)
		v1 = v00 + cs00 * A* f
		vel = np.dstack(velFunc(v1, z))
		return {'pres': p00 + gamma * p00 * A* f  , 'rho': rho00 + rho00 *A* f , 'vel': vel }
	else:
		print("analitical solution not implemented set plotAnalitical = False in notifier_params.py")
		sys.exit(0)

		
	



#boundary conditions:	
from sound_wave_params import periodicType
		
if periodicType == "repeat":

	def lrBoundaryConditionsPresRho(array, skip=0):
		n = array.shape[0] - 1
		#insert rows	
		array = np.insert(array, 0,  array[n-skip,:], axis = 0)
		array = np.insert(array, n+2,  array[1+skip,:], axis = 0)
		#insert columns	
		array = np.insert(array, 0, array[:,n-skip], axis = 1)
		array = np.insert(array, n+2, array[:,1+skip], axis = 1)
		return array

	lrBoundaryConditionsVel = lrBoundaryConditionsPresRho

elif periodicType == "refl" or periodicType == "diff":
		

	def lrBoundaryConditionsPresRho(array, skip=0):
		n = array.shape[0] - 1
		fr = 2 * array[0,:] - array[1,:]
		lr = 2 * array[-1,:] - array[-2,:]
		array = np.insert(array, 0, fr, axis = 0)
		array = np.insert(array, n+2, lr, axis = 0)
		fc = 2 * array[:,0] - array[:,1]
		lc = 2 * array[:,-1] - array[:,-2]
		array = np.insert(array, 0, fc, axis = 1)
		array = np.insert(array, n+2, lc, axis = 1)
		#print("array shape2")
		#print(array.shape)
		#print("end")
		return array


if periodicType == "refl":
	def lrBoundaryConditionsVel(array, skip=0):
		n = array.shape[0] - 1
		if(skip==0):
			array = np.insert(array, 0,  -array[0,:], axis = 0)
			array = np.insert(array, n+2,  -array[-1,:], axis = 0)
			array = np.insert(array, 0,  -array[:,0], axis = 1)
			array = np.insert(array, n+2,  -array[:,-1], axis = 1)
		elif (skip==1):
			array[0,:] = 0
			array = np.insert(array, 0,  -array[2,:], axis = 0)
			array[:,-1] = 0			
			array = np.insert(array, n+2,  -array[-2,:], axis = 0)
			array[:,0] = 0
			array = np.insert(array, 0,  -array[:,2], axis = 1)
			array[:,-1] = 0			
			array = np.insert(array, n+2,  -array[:,-2], axis = 1)
		return array

elif periodicType == "diff":
	lrBoundaryConditionsVel = lrBoundaryConditionsPresRho


