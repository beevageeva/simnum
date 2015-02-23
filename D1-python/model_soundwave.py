import numpy as np
import sys
from constants import gamma
from base_model import BaseModel
from sound_wave_params import plotPresCurve, plotVelCurve, plotRhoCurve, markPoints, plotPresAn, plotRhoAn, plotVelAn, plotVelFFT, plotVelFFTAnal, mediumType, plotPresFFT
from initcond_soundwave import getCs00



showErr = True
#showErr = False
calcKc = True
#calcKc = False

#addMarkPoint = None
#uncomment this to add a new mark point
#the following for the wave packet
from sound_wave_packet_params import zc,W
#addMarkPoint = zc 
addMarkPoint = zc - W*1.4 #other point at the beginning of the packet


#plotCsMaxMin = True
plotCsMaxMin = False

#calcAmp = True
calcAmp = False


class Model(BaseModel):
	

	def __init__(self):
		BaseModel.__init__(self)
		if addMarkPoint:
			self.addMarkPoint = addMarkPoint
			print("addMarkPoint = %E, plotting on pres axis" % self.addMarkPoint)


	def getNewPoint(self, zval, dt):
		from common import displacedPoint, getZIndex
		from math import sqrt
		from sound_wave_params import v00, p00, periodicType
		zIndex = getZIndex(zval)	
		#from sound_wave_params import csSign 
		#I should not import from here: this should be used for generating initial conditions only
		# I have to calculate it from actual values
		#Imagine that it should work for a superposition of wave travelling right and left
		if (self.pres[zIndex]< p00 and self.vel[zIndex] > v00 ) or (self.pres[zIndex]> p00 and self.vel[zIndex] < v00 ):
			csSign = -1
		else:
			csSign = 1
		cs = sqrt(gamma * self.pres[zIndex] / self.rho[zIndex])
		v = v00 + csSign * cs	
		newz = displacedPoint(zval, v, dt, periodicType)
		return newz


	def updateValuesModel(self, dt, time):
		#markPoints used to plot maximum and minumum for rho, pres and vel , in hom medium they travel at cs speed (phase vel), but in inhom at group speed
		if(markPoints):
			if mediumType == "homog":
				self.maxPresZ = self.getNewPoint(self.maxPresZ,dt)
				print("Homog medium : calculating max pres analytically(movind at cs speed) at z = %E" % self.maxPresZ)
			else:
				newZ = self.z[np.argmax(self.pres)]	
				maxSpeed = (newZ - self.maxPresZ)/dt
				#print("Inhomog medium : getting max pres at z = %E, travelling at speed = %E" % (newZ, maxSpeed))
				self.maxPresZ = newZ
		if hasattr(self, "addMarkPoint"):
			self.addMarkPoint = self.getNewPoint(self.addMarkPoint,dt)
		if(calcKc):
			#print("upd")
			self.presFFT = self.getPresFFTVals(False)
			F = self.presFFT[1]
			Y = self.presFFT[0]
			#first value is the mean
			#print("F=")
			#print(F)
			kc = abs(F[np.argmax(Y[1:])+1])
			#kc *=intlen #alreday multiplied all array
			#print("%E\t%E\t%E" % (cs,kc, cs*kc))   #!
			#TODO ONLY if markPoints = True ??
			if(markPoints):
				from initcond_soundwave import getCs00
				if mediumType == "homog":
					cs = getCs00()
					gradCs = 0
				else:
					cs = getCs00(self.maxPresZ)
					#TODO argmax calculated twice (see before in this function)
					gradientCs =  np.gradient(getCs00(self.z))
					gradCs = gradientCs[np.argmax(self.pres)]
					from sound_wave_packet_params import k0
					#print("%E\t%E\t%E" % (cs,kc, kc / cos(k0 * self.maxPresZ - self.omega0 * time)))   #!
					#print("%E\t%E\t%E" % (cs,kc, kc * self.maxPresZ))   #!
					#print("%E\t%E\t%E" % (cs,kc, (dt * gradCs)/(cs * (self.oldKc - kc))   ))   #!
					#print("%E\t%E\t%E" % (cs,kc, gradCs ))   #!
					#if(abs(self.oldKc - kc)>1e-10):
						#print("%E\t%E\t%E\t%E\t%E\t%E" % (cs,kc, gradCs, (self.oldKc - kc)/dt,  -kc * gradCs , (self.oldKc - kc)/dt * (cs / gradCs)  ))   #!
						#print("%E\t%E\t%E\t%E\t%E\t%E" % (cs,kc, gradCs, (self.oldKc - kc)/dt,  -self.oldKc * gradCs , (self.oldKc - kc)/dt * (cs / gradCs)  ))   #!
					from initcond_soundwave import csderAnal
					print("numkc = %E\ttkc = %E" % (kc, k0 * np.mean(np.exp(-csderAnal(self.z)  * time) )) )   #!
					if ('maxSpeed' in vars()) :
						print("Max speed = %e" % maxSpeed)
						#pass	

				#print("%E\t%E\t%E" % (cs,kc, cs + gradCs * kc))   #!
			#from common import getDz
			#from analyze_functions import getIndexRightAlmost0
			#indR = getIndexRightAlmost0(abs(Y), getDz()*0.005, 1)
			#print("width fourier pres = %E" % (2.0*kc))  #!
				
		if(calcAmp):
			print("pres amp = %E" % (np.max(self.pres) - np.min(self.pres)))
			print("vel amp = %E" % (np.max(self.vel) - np.min(self.vel)))
			print("rho amp = %E" % (np.max(self.rho) - np.min(self.rho)))


	def updateValuesNotifier(self, dt, time):
		#TODO simpl
		if(plotPresCurve):
			from initcond_soundwave import getPresCurveNumeric
		if(plotVelCurve):
			from initcond_soundwave import getVelCurveNumeric
		if(plotRhoCurve):
			from initcond_soundwave import getRhoCurveNumeric
		if(plotPresAn):
			from initcond_soundwave import getPresCurve,getPresAn
			presc = getPresCurve(self.z, time)
			anPres =  getPresAn(self.z, time, presc)
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
			from initcond_soundwave import getRhoCurve,getRhoAn
			rhoc = getRhoCurve(self.z, time)
			anRho =  getRhoAn(self.z, time, rhoc)
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
			from initcond_soundwave import getVelCurve,getVelAn
			velc = getVelCurve(self.z, time)
			anVel =  getVelAn(self.z, time, velc)
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

		self.notifier.updateValues("rho", rhoNewVals)
		self.notifier.updateValues("pres", presNewVals)
		self.notifier.updateValues("vel", velNewVals)
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
				print("Alreday calculated")
				presFFT = self.presFFT
			else:
				presFFT = self.getPresFFTVals(True)
			self.notifier.updateValues("presFFT", presFFT[0])

		if(markPoints):
			#only for pres	
			self.notifier.markPoint("pres", "maxPresZ", self.maxPresZ)
		if hasattr(self, "addMarkPoint"):
			self.notifier.markPoint("pres", "addMarkPoint", self.addMarkPoint)

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
		if(plotCsMaxMin):
			if(mediumType=="inhomog"):
				from initcond_soundwave import  getInitialFunctionMaxMinZIndex, getCs00
				r = getInitialFunctionMaxMinZIndex(self.z)
				minZIndex = r[0]
				maxZIndex = r[1]
				#print("minZIndex")
				#print(minZIndex)
				#print("maxZIndex")
				#print(maxZIndex)
				cs00 = getCs00(self.z)
				for alpha in np.arange(-2.5, 2.5, 0.5):
					vals = np.power(cs00, alpha)
					self.notifier.plotAxis("pres", np.multiply(self.pres[minZIndex] / vals[minZIndex],vals), "(min)%1.1f" % alpha)
					self.notifier.plotAxis("pres", np.multiply(self.pres[maxZIndex] / vals[maxZIndex],vals), "(max)%1.1f" % alpha)
					self.notifier.plotAxis("vel", np.multiply(self.vel[minZIndex] / vals[minZIndex],vals), "(min)%1.1f" % alpha)
					self.notifier.plotAxis("vel", np.multiply(self.vel[maxZIndex] / vals[maxZIndex],vals), "(max)%1.1f" % alpha)
					if(plotRhoCurve):
						from initcond_soundwave import getRhoCurveNumeric
						self.notifier.plotAxis("rhoCurve", getRhoCurveNumeric(np.multiply(self.rho[minZIndex] / vals[minZIndex],vals), self.z), "(min)%1.1f" % alpha)
						self.notifier.plotAxis("rhoCurve", getRhoCurveNumeric(np.multiply(self.rho[maxZIndex] / vals[maxZIndex],vals), self.z), "(max)%1.1f" % alpha)

			else:
				print("plotCsMin does not make sense with homogeneous medium")
		



