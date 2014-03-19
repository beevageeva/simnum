import numpy as np
import sys
from constants import gamma
from base_model import BaseModel




class Model(BaseModel):
	
	def __init__(self):
		BaseModel.__init__(self)


	def getNewPoint(self, zval, dt):
		from common import displacedPoint, getZIndex
		from math import sqrt
		zIndex = getZIndex(zval)	
		cs = sqrt(gamma * self.pres[zIndex] / self.rho[zIndex])
		from initcond_soundwave import getV00, getRho00, getP00
		v00 = getV00()
		p00 = getP00()
		rho00 = getRho00()
		#from sound_wave_params import csSign 
		#I should not import from here: this should be used for generating initial conditions only
		# I have to calculate it from actual values
		#Imagine that it should work for a superposition of wave travelling right and left
		if (self.rho[zIndex]< rho00 and self.vel[zIndex] > v00 ) or (self.rho[zIndex]> rho00 and self.vel[zIndex] < v00 ):
			csSign = -1
		else:
			csSign = 1	
		v = getV00() + csSign * cs	
		newz = displacedPoint(zval, v, dt)
		return newz


	def updateValues(self, dt, time):
		from initcond_soundwave import getRhoCurve, getPresCurve, getVelCurve, getRhoCurveNumeric, getPresCurveNumeric, getVelCurveNumeric, getRhoAn, getPresAn, getVelAn
		rhoc = getRhoCurve(self.z, time)
		presc = getPresCurve(self.z, time)
		velc = getVelCurve(self.z, time)
		self.notifier.updateValues("rho", [self.rho, getRhoAn(self.z, time, rhoc)])
		self.notifier.updateValues("pres", [self.pres, getPresAn(self.z, time, presc)])
		self.notifier.updateValues("vel", [self.vel, getVelAn(self.z, time, velc)])
		#If I don't want analitical solution plotted for velocity:
		#self.notifier.updateValues("vel", self.vel)
		if(self.plotPresCurve):
			self.notifier.updateValues("presCurve", [getPresCurveNumeric(self.pres), presc])
		if(self.plotVelCurve):
			self.notifier.updateValues("velCurve", [getVelCurveNumeric(self.vel), velc])
		if(self.plotRhoCurve):
			self.notifier.updateValues("rhoCurve", [getRhoCurveNumeric(self.rho), rhoc])
		if(self.markPoints):
			self.maxRhoZ = self.getNewPoint(self.maxRhoZ,dt)
			self.notifier.markPoint("rho", "maxRhoZ", self.maxRhoZ)
			self.maxPresZ = self.getNewPoint(self.maxPresZ,dt)
			self.notifier.markPoint("pres", "maxPresZ", self.maxPresZ)
			self.maxVelZ = self.getNewPoint(self.maxVelZ,dt)
			self.notifier.markPoint("vel", "maxVelZ", self.maxVelZ)


	def getInitialValues(self):
		return [[self.pres, self.pres], [self.rho, self.rho], [self.vel, self.vel]]
		#If I don't want analitical solution plotted for velocity:
		#return  [[self.rho, self.rho], [self.pres, self.pres], self.vel]

	def additionalInit(self):
		#plot Curves of pression , vel, density
		from notifier_params import plotPresCurve, plotVelCurve, plotRhoCurve, markPoints
		from initcond_soundwave import getRhoCurveNumeric, getPresCurveNumeric, getVelCurveNumeric
		self.plotPresCurve = plotPresCurve
		self.plotVelCurve = plotVelCurve
		self.plotRhoCurve = plotRhoCurve
		self.markPoints = markPoints
		if(self.plotPresCurve):
			self.notifier.addGraph("presCurve", [getPresCurveNumeric(self.pres), getPresCurveNumeric(self.pres)])
		if(self.plotVelCurve):
			self.notifier.addGraph("velCurve", [getVelCurveNumeric(self.vel), getVelCurveNumeric(self.vel)])
		if(self.plotRhoCurve):
			self.notifier.addGraph("rhoCurve", [getRhoCurveNumeric(self.rho),getRhoCurveNumeric(self.rho)])
		if(self.markPoints):
			from initcond_soundwave import  getInitialFunctionMaxZ
			r = getInitialFunctionMaxZ(self.z)
			self.maxRhoZ = r["rho"]
			self.maxPresZ = r["pres"]
			self.maxVelZ = r["vel"]
			self.notifier.markPoint("rho", "maxRhoZ", self.maxRhoZ)
			self.notifier.markPoint("pres", "maxPresZ", self.maxPresZ)
			self.notifier.markPoint("vel", "maxVelZ", self.maxVelZ)



