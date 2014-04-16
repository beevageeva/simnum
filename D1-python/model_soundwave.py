import numpy as np
import sys
from constants import gamma
from base_model import BaseModel
from sound_wave_params import plotPresCurve, plotVelCurve, plotRhoCurve, markPoints, plotPresAn, plotRhoAn, plotVelAn, plotVelFFT




class Model(BaseModel):
	
	def __init__(self):
		BaseModel.__init__(self)


	def getNewPoint(self, zval, dt):
		from common import displacedPoint, getZIndex
		from math import sqrt
		from initcond_soundwave import getV00, getRho00, getP00
		v00 = getV00()
		p00 = getP00()
		rho00 = getRho00()
		zIndex = getZIndex(zval)	
		#from sound_wave_params import csSign 
		#I should not import from here: this should be used for generating initial conditions only
		# I have to calculate it from actual values
		#Imagine that it should work for a superposition of wave travelling right and left
		if (self.rho[zIndex]< rho00 and self.vel[zIndex] > v00 ) or (self.rho[zIndex]> rho00 and self.vel[zIndex] < v00 ):
			csSign = -1
		else:
			csSign = 1
		cs = sqrt(gamma * self.pres[zIndex] / self.rho[zIndex])
		v = v00 + csSign * cs	
		newz = displacedPoint(zval, v, dt)
		return newz


	def updateValuesModel(self, dt, time):
		if(markPoints):
			self.maxRhoZ = self.getNewPoint(self.maxRhoZ,dt)
			self.maxPresZ = self.getNewPoint(self.maxPresZ,dt)
			self.maxVelZ = self.getNewPoint(self.maxVelZ,dt)

	def updateValuesNotifier(self, dt, time):
		from initcond_soundwave import getRhoCurve, getPresCurve, getVelCurve, getRhoCurveNumeric, getPresCurveNumeric, getVelCurveNumeric, getRhoAn, getPresAn, getVelAn
		if(plotPresAn):
			presc = getPresCurve(self.z, time)
			presNewVals = [self.pres, getPresAn(self.z, time, presc)]
			if(plotPresCurve):
				presCurveNewVals = [getPresCurveNumeric(self.pres), presc]
		else:
			presNewVals = self.pres
			if(plotPresCurve):
				presCurveNewVals = getPresCurveNumeric(self.pres)
		if(plotRhoAn):
			rhoc = getRhoCurve(self.z, time)
			rhoNewVals = [self.rho, getRhoAn(self.z, time, rhoc)]
			if(plotRhoCurve):
				rhoCurveNewVals = [getRhoCurveNumeric(self.rho), rhoc]
		else:
			rhoNewVals = self.rho
			if(plotRhoCurve):
				rhoCurveNewVals = getRhoCurveNumeric(self.rho)
		if(plotVelAn):
			velc = getVelCurve(self.z, time)
			velNewVals = [self.vel, getPresAn(self.z, time, velc)]
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
			self.notifier.updateFFTAxis("velFFT", self.vel)
		if(markPoints):
			self.notifier.markPoint("rho", "maxRhoZ", self.maxRhoZ)
			self.notifier.markPoint("pres", "maxPresZ", self.maxPresZ)
			self.notifier.markPoint("vel", "maxVelZ", self.maxVelZ)

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
		#plot Curves of pression , vel, density
		from initcond_soundwave import getRhoCurveNumeric, getPresCurveNumeric, getVelCurveNumeric
		if(plotPresCurve):
			self.notifier.addGraph("presCurve", [getPresCurveNumeric(self.pres), getPresCurveNumeric(self.pres)])
		if(plotVelCurve):
			self.notifier.addGraph("velCurve", [getVelCurveNumeric(self.vel), getVelCurveNumeric(self.vel)])
		if(plotRhoCurve):
			self.notifier.addGraph("rhoCurve", [getRhoCurveNumeric(self.rho),getRhoCurveNumeric(self.rho)])
		if(plotVelFFT):
			self.notifier.addFFTAxis("velFFT", self.vel)
		if(markPoints):
			from initcond_soundwave import  getInitialFunctionMaxZ
			r = getInitialFunctionMaxZ(self.z)
			self.maxRhoZ = r["rho"]
			self.maxPresZ = r["pres"]
			self.maxVelZ = r["vel"]
			self.notifier.markPoint("rho", "maxRhoZ", self.maxRhoZ)
			self.notifier.markPoint("pres", "maxPresZ", self.maxPresZ)
			self.notifier.markPoint("vel", "maxVelZ", self.maxVelZ)



