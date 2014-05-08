import numpy as np
import sys
from constants import gamma
from base_model import BaseModel
from sound_wave_params import plotPresCurve, plotVelCurve, plotRhoCurve, markPoints, plotPresAn, plotRhoAn, plotVelAn, plotVelFFT




#showErr = True
showErr = False
#calcKc = True
calcKc = False


class Model(BaseModel):
	

	def __init__(self):
		BaseModel.__init__(self)
		#uncomment this to add a new mark point
		#the following for the wave packet
		from sound_wave_defined_params import zc,W
		self.addMarkPoint = zc 
		#self.addMarkPoint = zc - 0.98*W #other point at the beginning of the packet
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
		if(markPoints):
			self.maxRhoZ = self.getNewPoint(self.maxRhoZ,dt)
			self.maxPresZ = self.getNewPoint(self.maxPresZ,dt)
			self.maxVelZ = self.getNewPoint(self.maxVelZ,dt)
		if hasattr(self, "addMarkPoint"):
			if(calcKc):
				from scipy.fftpack import fft,fftfreq#forFourierTransform
				numPoints = len(self.z)
				Y=fft(self.pres)/(numPoints)
				kc = np.max(abs(Y))
				from sound_wave_params import mediumType
				from initcond_soundwave import getCs00
				if mediumType == "homog":
					cs = getCs00()
				else:
					cs = getCs00(self.addMarkPoint)
				print("Product of kc*cs const?? %E" % cs * kc)
			self.addMarkPoint = self.getNewPoint(self.addMarkPoint,dt)


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
			#TODO it is calculated every time 
			from initcond_soundwave import getVelFFTAn
			getVelFFTAn = None
			self.notifier.updateFFTAxis("velFFT", self.vel, getVelFFTAn)
		if(markPoints):
			self.notifier.markPoint("rho", "maxRhoZ", self.maxRhoZ)
			self.notifier.markPoint("pres", "maxPresZ", self.maxPresZ)
			self.notifier.markPoint("vel", "maxVelZ", self.maxVelZ)
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

		


	def additionalInit(self):
		#TODO all initial values of curves are calcutaed 2 times: 2 function calls for each
		#plot Curves of pression , vel, density
		if(plotPresCurve):
			from initcond_soundwave import getPresCurveNumeric
			self.notifier.addGraph("presCurve",[getPresCurveNumeric(self.pres), getPresCurveNumeric(self.pres)] if plotPresAn else getPresCurveNumeric(self.pres))
		if(plotVelCurve):
			from initcond_soundwave import getVelCurveNumeric
			self.notifier.addGraph("velCurve", [getVelCurveNumeric(self.vel), getVelCurveNumeric(self.vel)] if plotVelAn else getVelCurveNumeric(self.vel) )
		if(plotRhoCurve):
			from initcond_soundwave import getRhoCurveNumeric
			self.notifier.addGraph("rhoCurve", [getRhoCurveNumeric(self.rho,self.z),getRhoCurveNumeric(self.rho, self.z)] if plotRhoAn else getRhoCurveNumeric(self.rho, self.z))
		if(plotVelFFT):
			from initcond_soundwave import getVelFFTAn
			getVelFFTAn = None
			self.notifier.addFFTAxis("velFFT", self.vel, getVelFFTAn)
		if(markPoints):
			from initcond_soundwave import  getInitialFunctionMaxZ
			r = getInitialFunctionMaxZ(self.z)
			self.maxRhoZ = r["rho"]
			self.maxPresZ = r["pres"]
			self.maxVelZ = r["vel"]
			self.notifier.markPoint("rho", "maxRhoZ", self.maxRhoZ)
			self.notifier.markPoint("pres", "maxPresZ", self.maxPresZ)
			self.notifier.markPoint("vel", "maxVelZ", self.maxVelZ)
		from sound_wave_params import mediumType
		if(mediumType=="inhomog"):
			from initcond_soundwave import getCs00
			self.notifier.plotAxisTwin("vel",getCs00(self.z) , "cs00")




