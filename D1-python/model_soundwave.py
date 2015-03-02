import numpy as np
import sys
from constants import gamma
from base_model import BaseModel
from notifier_params import plotPresCurve, plotVelCurve, plotRhoCurve, markPoints, plotPresAn, plotRhoAn, plotVelAn, plotVelFFT, plotVelFFTAnal,  plotPresFFT
from soundwave_medium_params import mediumType, cs00



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


plotCsMaxMin = True
#plotCsMaxMin = False



calcNewZ = True

class Model(BaseModel):
	

	def __init__(self):
		BaseModel.__init__(self)
		if addMarkPoint:
			self.addMarkPoint = addMarkPoint
			print("addMarkPoint = %E, plotting on pres axis" % self.addMarkPoint)


	#csNumerical = False won't work with superposition same point and time (wave packet is not in this case)
	def getNewPoint(self, zval, dt, csNumerical = True):
		from common import getZIndex
		from math import sqrt
		from soundwave_medium_params import v00, p00
		from soundwave_boundary_conditions import getPeriodicX
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
			from soundwave_medium_type import cs00
			if(mediumType == "homog"):
				cs = cs00
			else:
				cs = cs00(zval)	
		v = v00 +  cs	
		newz = getPeriodicX(zval + v * dt)
		return newz


	def updateValuesModel(self, dt, time):
		#markPoints used to plot maximum and minumum for rho, pres and vel , in hom medium they travel at cs speed (phase vel), but in inhom at group speed
		if(markPoints):
			self.maxPresZ = self.getNewPoint(self.maxPresZ,dt)
			graphPresMaxZ = self.z[np.argmax(self.pres)]
		if hasattr(self, "addMarkPoint"):
			self.addMarkPoint = self.getNewPoint(self.addMarkPoint,dt)
		from soundwave_perturbation_params import perturbationType
		if(perturbationType == "one"):
			from soundwave_perturbation_params import functiontype
			if(functiontype == "wavepacket" and mediumType == "inhomog"):
				if(plotPresAn or plotVelAn or plotRhoAn):
					from analytic_solution import updateNewZ
					updateNewZ(self, dt, False)
				if(calcKc and not addMarkPoint is None):
					from sound_wave_packet_params import k0
					from soundwave_medium_params import cs00, csderAnal
					from common import getDz
					print("dz: %e, dt*minCs=%e" % (getDz(), dt * np.min(cs00(self.z))))
					self.presFFT = self.getPresFFTVals(False)
					F = self.presFFT[1]
					Y = self.presFFT[0]
					kc = abs(F[np.argmax(Y[1:])+1])
					cszp0 = cs00(addMarkPoint)
					cszp = 	cs00(self.addMarkPoint)
					kt = k0 *  cszp0/cszp
					ktnew =  k0 * np.exp(-csderAnal(self.addMarkPoint) * time)
					print("kc=%e,kt=%e,ktnew=%e,kc/kt = %e, kc*cs(zp(t))=%e,kt*cszp=%e,k0*csz0=%e" % (kc, kt,ktnew,kc / kt, kc * cszp, kt*cszp, k0*cszp0))

				


	def updateValuesNotifier(self, dt, time):
		#TODO simpl
		if(plotPresCurve):
			from initcond_soundwave import fromValsToCurvePres
		if(plotVelCurve):
			from initcond_soundwave import fromValsToCurveVel
		if(plotRhoCurve):
			from initcond_soundwave import fromValsToCurveRho
		if(plotPresAn or plotVelAn or plotRhoAn):
			from analytic_solution import getCurves, getVals
			#TODO
			curves = getCurves(self.z, time)
			anVals = getVals(self.z, curves)
		
		if(plotPresAn):
			presc = curves['pres']
			anPres =  anVals['pres']
			presNewVals = [self.pres, anPres]
			if(showErr):
				err = np.max(np.absolute(np.subtract(self.pres, anPres)))
				print("time=%E,max abs (self.pres - anPres) = %E" % (time,err))
			if(plotPresCurve):
				presCurveNewVals = [fromValsToCurvePres(self.pres), presc]
		else:
			presNewVals = self.pres
			if(plotPresCurve):
				presCurveNewVals = getPresCurveNumeric(self.pres)
		if(plotRhoAn):
			rhoc = curves["rho"]
			anRho =  anVals["rho"]
			rhoNewVals = [self.rho, anRho]
			if(showErr):
				err = np.max(np.absolute(np.subtract(self.rho, anRho)))
				print("max abs (self.rho - rhoAn) = %E" % err)
			if(plotRhoCurve):
				rhoCurveNewVals = [fromValsToCurveRho(self.rho, self.z), rhoc]
		else:
			rhoNewVals = self.rho
			if(plotRhoCurve):
				rhoCurveNewVals = fromValsToCurveRho(self.rho, self.z)
		if(plotVelAn):
			velc = curves["vel"]
			anVel =  anVals["vel"]
			velNewVals = [self.vel, anVel]
			if(showErr):
				err = np.max(np.absolute(np.subtract(self.vel, anVel)))
				print("max abs (self.vel - velAn) = %E" % err)
			if(plotVelCurve):
				velCurveNewVals = [fromValsToCurveVel(self.vel), velc]
		else:
			velNewVals = self.vel
			if(plotVelCurve):
				velCurveNewVals = fromValsToCurveVel(self.vel)

		newRhoZ = None
		newPresZ = None
		newVelZ = None
		if(mediumType == "inhomog"):
			from analytic_solution import shiftNow
			if not shiftNow:
				from analytic_solution import newZ
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
				presFFT = self.getPresFFTVals(False)

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
			from initcond_soundwave import fromValsToCurvePres
			prescvals = fromValsToCurveVel(self.pres)
			self.notifier.addGraph(self.z, "presCurve",[prescvals, prescvals] if plotPresAn else prescvals)
		if(plotVelCurve):
			from initcond_soundwave import fromValsToCurveVel
			velcvals = fromValsToCurveVel(self.vel, self.z)
			self.notifier.addGraph(self.z, "velCurve", [velcvals, velcvals] if plotVelAn else velcvals )
		if(plotRhoCurve):
			from initcond_soundwave import fromValsToCurveRho
			rcvals = fromValsToCurveRho(self.rho, self.z)
			self.notifier.addGraph(self.z, "rhoCurve", [rcvals,rcvals] if plotRhoAn else rcvals)
		if(plotVelFFT):
			vals = self.getVelFFTVals(True)
			self.notifier.addGraph(vals[-1], "velFFT", vals[0:-1] , "k" , 'None')
		if(plotPresFFT):
			#use the same boolean parameter as in update	
			vals = self.getPresFFTVals(False)
			print("PRES FFT")
			print(vals)
			self.notifier.addGraph(vals[-1], "presFFT", vals[0] , "k", 'None')
	
		if(markPoints):
			self.maxPresZ = self.z[np.argmax(self.pres)]	
			self.notifier.markPoint("pres", "maxPresZ", self.maxPresZ)
		if(mediumType=="inhomog"):
			from soundwave_medium_params import cs00
			self.notifier.plotAxisTwin(self.z, "vel",cs00(self.z) , "cs00")
			#VARLENGTH
			from soundwave_perturbation_params import perturbationType
			if perturbationType == "one":
				from soundwave_perturbation_params import functiontype
				if functiontype == "wavepacket":
					from sound_wave_packet_params import k0
					from soundwave_medium_params import getDensVarLength	
					cl =	getDensVarLength(self.z)
					print("wl = %e, cl = %e" % (1.0 /k0, cl))
		if(plotCsMaxMin):
			if(mediumType=="inhomog"):

				maxZIndex = np.argmax(self.pres)
				minZIndex = np.argmin(self.pres)

				from soundwave_medium_params import cs00 as cs00Func
				from soundwave_medium_params import p00, densFunc

				cs00 = cs00Func(self.z)
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
		



