import numpy as np
import sys
from constants import gamma



def getNotifier(notifierType, z, titles, iniValues):
	if(notifierType=="visual"):
		#this can be any other class implementing corresponding methods
		from visual_plot import VisualPlot
		return VisualPlot(z, titles, iniValues)
	elif(notifierType=="file"):
		from save_to_file import SaveToFile 
		return SaveToFile(z, titles, iniValues)

class Model:
	
	def __init__(self):
		from common import getZArray
		from initcond import getInitialPresRhoVel, getRhoCurveNumeric, getPresCurveNumeric, getVelCurveNumeric
		from alg import getInitialUcUe
		self.z = getZArray()
		r = getInitialPresRhoVel(self.z)	
		self.pres = r['pres']
		self.rho = r['rho']
		self.vel = r['vel']
		r = getInitialUcUe(self.rho, self.vel, self.pres)
		self.uc = r['uc']
		self.ue = r['ue']
		#flux arrays
		self.fm = None
		self.fc = None
		self.fe = None
		#plot Curves of pression , vel, density
		from notifier_params import plotPresCurve, plotVelCurve, plotRhoCurve, markPoints, notifierType
		self.plotPresCurve = plotPresCurve
		self.plotVelCurve = plotVelCurve
		self.plotRhoCurve = plotRhoCurve
		self.markPoints = markPoints
		self.notifier = getNotifier(notifierType, self.z, ["rho", "pres", "vel"], [[self.rho, self.rho], [self.pres, self.pres], [self.vel, self.vel]])
		#If I don't want analitical solution plotted for velocity:
		#self.notifier = getNotifier(notifierType, self.z, ["rho", "pres", "vel"], [[self.rho, self.rho], [self.pres, self.pres], self.vel)
		
		if(self.plotPresCurve):
			self.notifier.addGraph("presCurve", [getPresCurveNumeric(self.pres), getPresCurveNumeric(self.pres)])
		if(self.plotVelCurve):
			self.notifier.addGraph("velCurve", [getVelCurveNumeric(self.vel), getVelCurveNumeric(self.vel)])
		if(self.plotRhoCurve):
			self.notifier.addGraph("rhoCurve", [getRhoCurveNumeric(self.rho),getRhoCurveNumeric(self.rho)])
		if(self.markPoints):
			from initcond import  getInitialFunctionMaxZ
			r = getInitialFunctionMaxZ(self.z)
			self.maxRhoZ = r["rho"]
			self.maxPresZ = r["pres"]
			self.maxVelZ = r["vel"]
			self.notifier.markPoint("rho", "maxRhoZ", self.maxRhoZ)
			self.notifier.markPoint("pres", "maxPresZ", self.maxPresZ)
			self.notifier.markPoint("vel", "maxVelZ", self.maxVelZ)
		self.notifier.afterInit()


	def printVars(self, time):
		from constants import verbose
		if verbose:
			print("Before calculating at time=%4.3f\nz" % time)
			print(self.z)
			print("rho")
			print(self.rho)
			print("pres")
			print(self.pres)
			print("vel")
			print(self.vel)
			print("uc")
			print(self.uc)
			print("ue")
			print(self.ue)
			print("fm")
			print(self.fm)
			print("fc")
			print(self.fc)
			print("fe")
			print(self.fe)
			print("END")


	def getNewPoint(self, zval, dt):
		from common import displacedPoint, getZIndex
		from math import sqrt
		from sound_wave_params import csSign
		zIndex = getZIndex(zval)	
		cs = sqrt(gamma * self.pres[zIndex] / self.rho[zIndex])
		from initcond import getV00
		v = getV00() + csSign * cs	
		newz = displacedPoint(self.maxRhoZ, v, dt)
		return newz


	def updateValues(self, dt, time):
		from initcond import getRhoCurve, getPresCurve, getVelCurve, getRhoCurveNumeric, getPresCurveNumeric, getVelCurveNumeric, getRhoAn, getPresAn, getVelAn
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
		self.notifier.afterUpdateValues(time)



	def mainLoop(self, timeEnd):
		from alg import recalculateFluxes, getTimestep, recalculateU, recalculateVelPres,getInitialUcUe
		time = 0.0
		nstep = 0
		while(time<timeEnd):
			r = recalculateFluxes(self.rho, self.uc, self.ue, self.vel, self.pres)
			self.fm = r['fm']
			self.fc = r['fc']
			self.fe = r['fe']
			dt = getTimestep(self.vel, self.pres, self.rho)
			time+=dt
			nstep+=1
			#recalculate u at next step - time	
			self.rho = recalculateU(self.rho, self.fm, dt)	
			self.uc = recalculateU(self.uc, self.fc, dt)	
			self.ue = recalculateU(self.ue, self.fe, dt)
			#recalculate velocity and pression
			r = recalculateVelPres(self.rho, self.uc, self.ue)
			self.vel = r["vel"]	
			self.pres = r["pres"]

			self.printVars(time)
			from notifier_params import nstepsPlot
			if(nstep % nstepsPlot == 0):
				self.updateValues(nstepsPlot *  dt, time)
		self.notifier.finish()
