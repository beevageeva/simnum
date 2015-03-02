import numpy as np
import sys
from constants import gamma
from base_model import BaseModel
from sound_wave_params import plotPresCurve, plotVelCurve, plotRhoCurve, markPoints, plotPresAn, plotRhoAn, plotVelAn, plotVelFFT, plotVelFFTAnal, mediumType, plotPresFFT
from initcond_soundwave import getCs00



#showErr = True
showErr = False
calcKc = True
#calcKc = False

#addMarkPoint = None
#uncomment this to add a new mark point
#the following for the wave packet
from sound_wave_packet_params import zc,W
addMarkPoint = zc 
#addMarkPoint = zc - W*2 #other point at the beginning of the packet


#plotCsMaxMin = True
plotCsMaxMin = False

#calcAmp = True
calcAmp = False


calcNewZ = True

class Model(BaseModel):
	

	def __init__(self):
		BaseModel.__init__(self)
		if addMarkPoint:
			self.addMarkPoint = addMarkPoint
			print("addMarkPoint = %E, plotting on pres axis" % self.addMarkPoint)


	#csNumerical = False won't work with superposition same point and time (wave packet is not in this case)
	def getNewPoint(self, zval, dt, csNumerical = True):
		from common import displacedPoint, getZIndex
		from math import sqrt
		from sound_wave_params import v00, p00, periodicType
		if(csNumerical):
			zIndex = getZIndex(zval)	
			#from sound_wave_params import csSign 
			#I should not import from here: this should be used for generating initial conditions only
			# I have to calculate it from actual values
			#Imagine that it should work for a superposition of wave travelling right and left
			if (self.pres[zIndex]< p00 and self.vel[zIndex] > v00 ) or (self.pres[zIndex]> p00 and self.vel[zIndex] < v00 ):
				csSign = -1
			else:
				csSign = 1
			cs = csSign * sqrt(gamma * self.pres[zIndex] / self.rho[zIndex])
		else:
			from initcond_soundwave import getCs00
			if(mediumType == "homog"):
				cs = getCs00()
			else:
				cs = getCs00(zval)	
		v = v00 +  cs	
		newz = displacedPoint(zval, v, dt, periodicType)
		return newz


	def updateValuesModel(self, dt, time):
		#markPoints used to plot maximum and minumum for rho, pres and vel , in hom medium they travel at cs speed (phase vel), but in inhom at group speed
		if(markPoints):

			self.maxPresZ = self.getNewPoint(self.maxPresZ,dt)
			graphPresMaxZ = self.z[np.argmax(self.pres)]
			#print("pres max: %e==%e, dif=%e" % (self.maxPresZ, graphPresMaxZ, abs(self.maxPresZ - graphPresMaxZ)))


#			if mediumType == "homog":
#				self.maxPresZ = self.getNewPoint(self.maxPresZ,dt)
#				print("Homog medium : calculating max pres analytically(movind at cs speed) at z = %E" % self.maxPresZ)
#			else:
#				newZ = self.z[np.argmax(self.pres)]	
#				maxSpeed = (newZ - self.maxPresZ)/dt
#				#print("Inhomog medium : getting max pres at z = %E, travelling at speed = %E" % (newZ, maxSpeed))
#				self.maxPresZ = newZ
		if hasattr(self, "addMarkPoint"):
#			mp = self.addMarkPoint
			self.addMarkPoint = self.getNewPoint(self.addMarkPoint,dt)
#			if mediumType == "inhomog":
#				print("markPoint num = %e, markPoint from cs = %e" % (self.addMarkPoint, mp + dt * getCs00(mp)  ))
		
		#CALCNEWZ		
#		if calcNewZ:
#			from common import getZIndex
#			from initcond_soundwave import csderAnal
#			#TODO weave
#			#cs = getCs00(self.z)
#			#csder = csderAnal(self.z)
#			for index in range(self.z.shape[0]):
#				#zval = self.z[index]
#				
#				#newZ = self.getNewPoint(self.newZ[index],dt)
#				#newZIndex = getZIndex(newZ)
#				#print("CRAP %e" %  (self.K[newZIndex] * csderAnal(self.newZ[index]) * getCs00(self.newZ[index]) * dt))
#				#self.K[index] = self.K[index] - self.K[newZIndex] * csderAnal(self.newZ[index]) * getCs00(self.newZ[index]) * dt
#				#self.K[index] = self.K[index] - self.K[index] * csderAnal(self.newZ[index]) * getCs00(self.newZ[index]) * dt
#				#self.K[index] = self.K[index] - self.K[index] * csderAnal(self.newZ[index])  * dt
#				#self.K[index] +=(-np.exp(-csderAnal(self.newZ[index]) * getCs00(self.newZ[index])) + np.exp(-csderAnal(newZ) * getCs00(newZ)))*dt
#				#self.K[index] *=np.exp((csderAnal(self.newZ[index]) * getCs00(self.newZ[index]) - csderAnal(newZ) * getCs00(newZ))*dt)
#				#print("%e" % (getCs00(newZ) / getCs00(self.newZ[index]) ))
#				#self.K[index] = self.K[index] * getCs00(self.newZ[index]) / getCs00(newZ)
#				self.newZ[index] = self.getNewPoint(self.newZ[index],dt)
#				#self.K[index] = self.K[index] - self.K[index] * csderAnal(zval) * getCs00(zval) * dt
#				#self.K[index] = self.K[index] - self.K[index] * csderAnal(zval) * dt
#

		
			

		


		if(calcKc):
			#print("upd")
			#use initial markpoint
			if mediumType == "inhomog" and not addMarkPoint is None:
				#from initcond_soundwave import kAnal
				#k = kAnal(self.z, time)
				#print("kc = %e, mean anKc = %e, max anKc = %e" % (kc, np.mean(k), np.max(k)))
				#from initcond_soundwave import csderAnal, getX0Index
				#from common import getZIndex
				#indexMarkPoint = getZIndex(addMarkPoint)
				#cp = np.exp(-csderAnal(addMarkPoint) * getCs00(addMarkPoint) * time)   #NO
				#cp = np.exp(-csderAnal(addMarkPoint) * getCs00(self.addMarkPoint) * time) #NO
				#csd = csderAnal(self.z)
				#csdIndMin = np.argmin(csd)
				#show find
				#k = 60.0 * np.exp(-csderAnal(self.z) * time)
				#k = kAnal(self.z, time)
				#print("kc = %e, mean anKc = %e, max anKc = %e" % (kc, np.mean(k), np.max(k)))
				#x0AddMarkPoint =  getX0Index(self.z, time, self.addMarkPoint)
				#print("%e == %e" % (addMarkPoint, self.z[x0AddMarkPoint]))
				#cp = np.exp(-csderAnal(self.addMarkPoint) * time) 
				#cp = np.exp(-csderAnal(self.addMarkPoint) * time)
				#cs = getCs00(self.addMarkPoint) 
				#print("kc=%e,cp=%e,kc/cp=%e,cs=%e,cs*cp=%e,kc/(cp*cs)=%e,kc*cs=%e" % (kc, cp, kc / cp, cs, cs*cp, kc/(cp * cs), cs*kc))
				from common import getZIndex			
				#np.set_printoptions(threshold='nan')
				#print("self.newZIndex B")
				#print(self.newZIndex)
				#print("self.newZ B")
				#print(self.newZ)

				tempNewZ = np.zeros(self.z.shape)
				tempNewZIndex = np.zeros(self.z.shape)
				#going index forward same as cs : will overwrite the array TODO no temp vars
				for index in range(self.z.shape[0]):
					tempNewZ[index] = self.getNewPoint(self.newZ[index],dt, False)
					#resolution!
					tempNewZIndex[index] = getZIndex(self.getNewPoint(self.z[self.newZIndex[index]],dt, False))
#					if (index > 1000 and index <1020):
#						print("index=%d, firstz=%e, oldnewz=%e, newnewz=%e" % (index ,self.z[index], self.newZ[index],  tempNewZ[index] ))
				self.newZ = tempNewZ
				self.newZIndex = tempNewZIndex
				print("first mp = %d, newz = %e, num=%e, ini index = %d , new z index of mark point newZIndex: %d, newz: %d , num: %d" % (addMarkPoint, self.newZ[getZIndex(addMarkPoint)], self.addMarkPoint, getZIndex(addMarkPoint), self.newZIndex[getZIndex(addMarkPoint)], getZIndex(self.newZ[getZIndex(addMarkPoint)]), getZIndex(self.addMarkPoint) ))

				from initcond_soundwave import getCs00
				from common import getDz
				print("dz: %e, dt*minCs=%e" % (getDz(), dt * np.min(getCs00(self.z))))
				
				#print("zindex of addMarkPoint %e, zindex of new markpoint %e"   % (getZIndex(addMarkPoint), getZIndex(self.addMarkPoint)))

				#print("self.newZ[addMarkpoint] %e == %e self.addMarkPoint " % (self.newZ[getZIndex(addMarkPoint)],self.addMarkPoint ))	
				


				#print("%e == %e" % (csderAnal(addMarkPoint) * getCs00(addMarkPoint), csderAnal(self.addMarkPoint))) NO
				#print("%e == %e" % (csderAnal(addMarkPoint) * getCs00(self.addMarkPoint), csderAnal(self.addMarkPoint))) NO
				#print("%e" % (csderAnal(addMarkPoint) * getCs00(self.addMarkPoint)/ csderAnal(self.addMarkPoint)))
#				from initcond_soundwave import kAnal
#				k = kAnal(self.z, time)
#				print("kanal mean = %e , kAn max = %e" % (np.mean(k), np.max(k)))			

	
				
		if(calcAmp):
			print("pres amp = %E" % (np.max(self.pres) - np.min(self.pres)))
			print("vel amp = %E" % (np.max(self.vel) - np.min(self.vel)))
			print("rho amp = %E" % (np.max(self.rho) - np.min(self.rho)))


	def updateValuesNotifier(self, dt, time):
		if(mediumType == "inhomog" and (plotPresAn or plotRhoAn or plotVelAn)):
			from initcond_soundwave import parametricCurves, csderAnal
			csZ0 = getCs00(self.z)
			csZt = getCs00(self.newZ)
			#csZt2 = getCs00(self.z[self.newZIndex])
			from sound_wave_packet_params import getSoundWaveGaussFunction, zc, W, k0
			from sound_wave_params import p00, A, densFunc, v00
			#gamma * p00/densFunc(z(t)) = cs(z(t)) ** 2	
			from constants import gamma, z0, zf
			from math import pi	
			from common import getZIndex

			def newZIndexInverse(zvalIndex):
				for	index in range(self.z.shape[0]):
					if self.newZIndex[index] == zvalIndex :
						return index
			
			def reshiftInv(array):
				res = np.zeros(array.shape)
				for index in range(0, array.shape[0]):
					res[newZIndexInverse(index)] = array[index]
				return res
			def reshift(array):
				res = np.zeros(array.shape)
				for index in range(0, array.shape[0]):
					res[self.newZIndex[index]] = array[index]
				return res
			def reshift2(array):
				res = np.zeros(array.shape)
				for index in range(0, array.shape[0]):
					res[getZIndex(self.newZ[index])] = array[index]
				return res
	
			self.presFFT = self.getPresFFTVals(False)
			F = self.presFFT[1]
			Y = self.presFFT[0]
			#first value is the mean
			#print("F=")
			#print(F)
			kc = abs(F[np.argmax(Y[1:])+1])
			ampIni = getSoundWaveGaussFunction(zc, W)(self.z)
			shiftNow = True
			#presAmp = ampIni(self.z) * csZ0 ** (0.5) * csZt ** (-0.5) * (p00 * gamma * A)
			#rhoAmp = presAmp * csZt**(-2)
			#velAmp  = presAmp * (1.0 / (p00 * gamma)) * csZt
			#k = k0 * csZ0 / csZt
			#k2 = k0 * csZ0 / csZt2
			#kSh = reshift(k)
			#kSh2 = reshift2(k)
			#print("update not K mark point %e , ksh self add  markpoint = %e " % (k[getZIndex(addMarkPoint)],kSh[getZIndex(self.addMarkPoint)] ))
			#print("update not2 K mark point %e , ksh self add  markpoint = %e " % (k2[getZIndex(addMarkPoint)],kSh2[getZIndex(self.addMarkPoint)] ))
			#curve = np.zeros(self.z.shape) INTERP NEXT!
			curve = np.full(self.z.shape, np.nan)
			for index in range(self.z.shape[0]):
				#curve[newZIndexInverse(index)] = ampIni[index] * csZ0[index] ** (0.5) * csZt[index] ** (-0.5) *  A * np.cos( 2 * pi / (zf - z0) *(k[index] * self.z[self.newZIndex[index]] - k0 * csZ0[index] * time - k0 * z0))
				#curve[self.newZIndex[index]] = ampIni[index] * csZ0[index] ** (0.5) * csZt[index] ** (-0.5) *  A * np.cos( 2 * pi / (zf - z0) *(k[index] * self.z[self.newZIndex[index]] - k0 * csZ0[index] * time - k0 * z0))
				#curve[index] = ampIni[index] * csZ0[index] ** (0.5) * csZt2[index] ** (-0.5) *  A * np.cos( 2 * pi / (zf - z0) *(k2[index] * self.newZ[index] - k0 * csZ0[index] * time - k0 * z0))
				if shiftNow:
					cIndex = getZIndex(self.newZ[index])
					if(np.isnan(curve[cIndex])):
						curve[cIndex] = 0
					curve[cIndex] += ampIni[index] * csZ0[index] ** (0.5) * csZt[index] ** (-0.5) *  A * np.cos( 2 * pi * k0 * (csZ0[index]/csZt[index]) / (zf - z0) * ( self.newZ[index] - csZt[index] * time) - 2 * pi * k0 * z0 / (zf - z0))
				else:
					curve[index] = ampIni[index] * csZ0[index] ** (0.5) * csZt[index] ** (-0.5) *  A * np.cos( 2 * pi * k0 * (csZ0[index]/csZt[index]) / (zf - z0) * ( self.newZ[index] - csZt[index] * time) - 2 * pi * k0 * z0 / (zf - z0))
			
			#argSort = np.argsort(self.newZ)	
			#curve = curve[argSort]

			if shiftNow:
				#interpolate nan values!
				nans, x= np.isnan(curve), lambda z: z.nonzero()[0]
				curve[nans]= np.interp(x(nans), x(~nans), curve[~nans])

				rhoIni = densFunc(self.z) #= gamma * p00  cs(z(t))** (-2)  no shift
				velAn =   v00 + csZ0 * curve 
				newZ = None
			else:
				rhoIni = densFunc(self.newZ) #= gamma * p00  cs(z(t))** (-2) shift afterwards
				velAn =   v00 + csZt * curve 
				newZ = [self.z, self.newZ]
			presAn =  p00 + p00 * gamma * curve
			rhoAn =   rhoIni + rhoIni * curve

			#anCurves = parametricCurves(self.z, self.newZ, time)

			#print("MP num, a a2 %e == %e == %e" % (self.addMarkPoint, self.newZ[getZIndex(addMarkPoint)], self.z[self.newZIndex[getZIndex(addMarkPoint)]] ))	
			cszp0 = getCs00(addMarkPoint)
			cszp = 	getCs00(self.addMarkPoint)
			kt = 60.0 *  cszp0/cszp
			ktnew =  k0 * (csZ0[getZIndex(addMarkPoint)]/csZt[getZIndex(addMarkPoint)]) 
			#ampMaxPres =p00 + p00* A*ampIni(addMarkPoint) * cszp ** (-0.5) * cszp0 ** 0.5
			#print("maxPres=%e, ampPres=%e, az0/cszp=%e" % (np.max(self.pres) , ampMaxPres, ampIni(addMarkPoint)/cszp**(-0.5)))
			#kt = self.K[getZIndex(self.addMarkPoint)]
			#kt1 = self.K[getZIndex(self.addMarkPoint)]
			#print("kc=%e,kt=%e,kc/kt = %e,kt1=%e,kc/kt1=%e,kc*cs(zp(t))=%e" % (kc, kt, kc / kt, kt1, kc/kt1,kc * cs))
			print("kc=%e,kt=%e,ktnew=%e,ktExp=%e,kc/kt = %e, kc*cs(zp(t))=%e,kt*cszp=%e,k0*csz0=%e" % (kc, kt,ktnew, k0 * np.exp(-csderAnal(self.addMarkPoint) * time),kc / kt, kc * cszp, kt*cszp, k0*cszp0))
			print("MP num, a  %e == %e " % (self.addMarkPoint, self.newZ[getZIndex(addMarkPoint)]))	

			print("pres mark point %e == %e" % (self.pres[getZIndex(self.addMarkPoint)], presAn[getZIndex(self.addMarkPoint)] ))
			#curve = anCurves['curve']	
	

		#TODO simpl
		if(plotPresCurve):
			from initcond_soundwave import getPresCurveNumeric
		if(plotVelCurve):
			from initcond_soundwave import getVelCurveNumeric
		if(plotRhoCurve):
			from initcond_soundwave import getRhoCurveNumeric
		if(plotPresAn):
			if(mediumType == "homog"):
				from initcond_soundwave import getPresCurve,getPresAn
				presc = getPresCurve(self.z, time)
				anPres =  getPresAn(self.z, time, presc)
			else:
				presc = curve
				anPres = presAn
			presNewVals = [self.pres, anPres]
			if(showErr):
				err = np.max(np.absolute(np.subtract(self.pres, anPres)))
				print("time=%E,max abs (self.pres - anPres) = %E" % (time,err))
			if(plotPresCurve):
				presCurveNewVals = [getPresCurveNumeric(self.pres), presc]
		else:
			presNewVals = self.pres
			if(plotPresCurve):
				presCurveNewVals = getPresCurveNumeric(self.pres)
		if(plotRhoAn):
			if(mediumType == "homog"):	
				from initcond_soundwave import getRhoCurve,getRhoAn
				rhoc = getRhoCurve(self.z, time)
				anRho =  getRhoAn(self.z, time, rhoc)
			else:
				rhoc = curve
				anRho = rhoAn
			rhoNewVals = [self.rho, anRho]
			if(showErr):
				err = np.max(np.absolute(np.subtract(self.rho, anRho)))
				print("max abs (self.rho - rhoAn) = %E" % err)
			if(plotRhoCurve):
				rhoCurveNewVals = [getRhoCurveNumeric(self.rho, self.z), rhoc]
		else:
			rhoNewVals = self.rho
			if(plotRhoCurve):
				rhoCurveNewVals = getRhoCurveNumeric(self.rho, self.z)
		if(plotVelAn):
			if(mediumType == "homog"):
				from initcond_soundwave import getVelCurve,getVelAn
				velc = getVelCurve(self.z, time)
				anVel =  getVelAn(self.z, time, velc)
			else:
				velc = curve
				anVel = velAn
			velNewVals = [self.vel, anVel]
			if(showErr):
				err = np.max(np.absolute(np.subtract(self.vel, anVel)))
				print("max abs (self.vel - velAn) = %E" % err)
			if(plotVelCurve):
				velCurveNewVals = [getVelCurveNumeric(self.vel), velc]
		else:
			velNewVals = self.vel
			if(plotVelCurve):
				velCurveNewVals = getVelCurveNumeric(self.vel)

		if(mediumType == "homog"):
			newRhoZ = None
			newPresZ = None
			newVelZ = None
		else:	
			newRhoZ = newZ
			newPresZ = newZ
			newVelZ = newZ
		

		self.notifier.updateValues("rho", rhoNewVals, newRhoZ)
		self.notifier.updateValues("pres", presNewVals, newPresZ)
		self.notifier.updateValues("vel", velNewVals, newVelZ)
		if(plotPresCurve):
			self.notifier.updateValues("presCurve", presCurveNewVals)
		if(plotVelCurve):
			self.notifier.updateValues("velCurve", velCurveNewVals)
		if(plotRhoCurve):
			self.notifier.updateValues("rhoCurve", rhoCurveNewVals)
		if(plotVelFFT):
			self.notifier.updateValues("velFFT", self.getVelFFTVals(True)[0:-1])
		if(plotPresFFT):
			#if ('presFFT' in vars()):
			if (hasattr(self, 'presFFT')):
				presFFT = self.presFFT
			else:
				print("Calculate presFFT ")
				presFFT = self.getPresFFTVals(True)
			self.notifier.updateValues("presFFT", presFFT[0])

		if(markPoints):
			#only for pres	
			self.notifier.markPoint("pres", "maxPresZ", self.maxPresZ)
		if hasattr(self, "addMarkPoint"):
			self.notifier.markPoint("pres", "addMarkPoint", self.addMarkPoint)
#		from analyze_functions import getFirstIndexDifferentLeft, getFirstIndexDifferentRight
#		delta = 0.1
#		self.notifier.markPoint("pres", "LB",getFirstIndexDifferentLeft(self.pres, delta, zc))
#		self.notifier.markPoint("pres", "RB",getFirstIndexDifferentRight(self.pres, delta, None))
		

	def getInitialValues(self):
		if plotPresAn:
			presIniVal = [self.pres, self.pres]
		else:
			presIniVal = self.pres
		if plotRhoAn:
			rhoIniVal = [self.rho, self.rho]
		else:
			rhoIniVal = self.rho
		if plotVelAn:
			velIniVal = [self.vel, self.vel]
		else:
			velIniVal = self.vel
		#return [[self.pres, self.pres], [self.rho, self.rho], [self.vel, self.vel]]
		#If I don't want analitical solution plotted for velocity:
		#return  [[self.pres, self.pres], [self.rho, self.rho], self.vel]
		return  [presIniVal, rhoIniVal, velIniVal]

		
	def getVelFFTVals(self, middlePoints):
		if(middlePoints):
			vals = (self.vel[1:] + self.vel[:-1]) / 2
		else:
			vals = self.vel	
		from constants import z0, zf
		intlen = zf - z0
		from scipy.fftpack import fft,fftfreq
		Y=fft(vals)/len(vals)
		F=fftfreq(len(vals), self.z[1] - self.z[0]) 
		#print("max freq in plotVelFFT %e" % intlen * F[np.argmax(Y[1:]) + 1 ])
		if(plotVelFFTAnal):
			from initcond_soundwave import getVelFFTAn
			vals = [intlen * abs(Y), abs(getVelFFTAn(F)), F * intlen]
			#vals = [intlen * abs(Y), getVelFFTAn(F), F ]
		else:
			vals = [intlen * abs(Y), F * intlen]	
			#vals = [intlen * abs(Y), F ]	
		return vals

	def getPresFFTVals(self, middlePoints):
		if(middlePoints):
			vals = (self.pres[1:] + self.pres[:-1]) / 2
		else:
			vals = self.pres	
		from constants import z0, zf
		intlen = zf - z0
		from scipy.fftpack import fft,fftfreq
		Y=fft(vals)/len(vals)
		F=fftfreq(len(vals), self.z[1] - self.z[0])
		#print("getPresFFT")
		#print(abs(Y)) 
		#print("max freq in plotPresFFT %e " % intlen * F[np.argmax(Y[1:]) + 1 ] )
		return  [intlen * abs(Y), F * intlen]	
		#return [intlen * abs(Y), F ]	
	
	

	def additionalInit(self):
		#TODO all initial values of curves are calcutaed 2 times: 2 function calls for each
		#plot Curves of pression , vel, density
		if(plotPresCurve):
			from initcond_soundwave import getPresCurveNumeric
			self.notifier.addGraph(self.z, "presCurve",[getPresCurveNumeric(self.pres), getPresCurveNumeric(self.pres)] if plotPresAn else getPresCurveNumeric(self.pres))
		if(plotVelCurve):
			from initcond_soundwave import getVelCurveNumeric
			self.notifier.addGraph(self.z, "velCurve", [getVelCurveNumeric(self.vel), getVelCurveNumeric(self.vel)] if plotVelAn else getVelCurveNumeric(self.vel) )
		if(plotRhoCurve):
			from initcond_soundwave import getRhoCurveNumeric
			self.notifier.addGraph(self.z, "rhoCurve", [getRhoCurveNumeric(self.rho,self.z),getRhoCurveNumeric(self.rho, self.z)] if plotRhoAn else getRhoCurveNumeric(self.rho, self.z))
		if(plotVelFFT):
			vals = self.getVelFFTVals(True)
			self.notifier.addGraph(vals[-1], "velFFT", vals[0:-1] , "k" , 'None')
		if(plotPresFFT):
			vals = self.getPresFFTVals(False)
			print("PRES FFT")
			print(vals)
			self.notifier.addGraph(vals[-1], "presFFT", vals[0] , "k", 'None')
	
		if(markPoints):
			if mediumType == "homog":
				from initcond_soundwave import  getInitialFunctionMaxZ
				r = getInitialFunctionMaxZ(self.z)
				#mark only pres point
				self.maxPresZ = r["pres"]
			else:
				from initcond_soundwave import  getInitialFunctionMaxMinZIndex, getCs00
				r = getInitialFunctionMaxMinZIndex(self.z)
				self.maxPresZ = self.z[r[1]]
				print("Group velocity??????  = %e" % getCs00(r[1]))

			self.notifier.markPoint("pres", "maxPresZ", self.maxPresZ)
		if(mediumType=="inhomog"):
			self.notifier.plotAxisTwin(self.z, "vel",getCs00(self.z) , "cs00")
			#VARLENGTH
			from sound_wave_params import functiontype
			if functiontype == "wavepacket":
					from sound_wave_packet_params import k0
					from sound_wave_params import getDensVarLength	
					cl =	getDensVarLength(self.z)
					print("wl = %e, cl = %e" % (1.0 /k0, cl))
		if(plotCsMaxMin):
			if(mediumType=="inhomog"):
				from initcond_soundwave import  getInitialFunctionMaxMinZIndex, getCs00
				from sound_wave_params  import  p00, densFunc
				r = getInitialFunctionMaxMinZIndex(self.z)
				minZIndex = r[0]
				maxZIndex = r[1]
				#print("minZIndex")
				#print(minZIndex)
				#print("maxZIndex")
				#print(maxZIndex)
				cs00 = getCs00(self.z)
				#for alpha in np.arange(-2.5, 2.5, 0.5):
				for alpha in [-0.5]:
					vals = np.power(cs00, alpha)
					self.notifier.plotAxis(self.z,"pres",  p00 + ((self.pres[minZIndex] - p00)/ vals[minZIndex]) * vals, "(min)%1.1f" % alpha)
					self.notifier.plotAxis(self.z,"pres",  p00 + ((self.pres[maxZIndex] - p00)/ vals[maxZIndex]) * vals, "(max)%1.1f" % alpha)
				for alpha in [ 0.5]:
					vals = np.power(cs00, alpha)
					self.notifier.plotAxis(self.z,"vel", np.multiply(self.vel[minZIndex] / vals[minZIndex],vals), "(min)%1.1f" % alpha)
					self.notifier.plotAxis(self.z,"vel", np.multiply(self.vel[maxZIndex] / vals[maxZIndex],vals), "(max)%1.1f" % alpha)
					if(plotRhoCurve):
						for alpha in [-0.5]:
							vals = np.power(cs00, alpha)
							self.notifier.plotAxis(self.z,"rhoCurve", ((self.rho[minZIndex] - densFunc(self.z[minZIndex])) / (vals[minZIndex] * densFunc(self.z[minZIndex])))* vals , "(min)%1.1f" % alpha)
							self.notifier.plotAxis(self.z,"rhoCurve", ((self.rho[maxZIndex] - densFunc(self.z[maxZIndex])) / (vals[maxZIndex] * densFunc(self.z[maxZIndex])))* vals , "(max)%1.1f" % alpha)

			else:
				print("csmin max only inhom")
		if calcNewZ and mediumType == "inhomog":
			self.newZIndex = np.arange(self.z.shape[0])
			self.newZ = np.array(self.z)
		



