import numpy as np
import sys
from constants import gamma, problemType

if problemType == "soundwave":
	import initcond_soundwave as initcond
elif problemType == "riemann":
	import initcond_riemann as initcond
else:
	print("problemtype %s not implemented" % problemType)
	sys.exit(0)


def getNotifier(notifierType, z, titles, iniValues):
	if(notifierType=="visual"):
		#this can be any other class implementing corresponding methods
		from visual_plot import VisualPlot
		return VisualPlot(z, titles, iniValues)
	elif(notifierType=="file"):
		from save_to_file import SaveToFile 
		return SaveToFile(z, titles, iniValues)

class BaseModel:
	
	def __init__(self):
		from common import getZArray
		from alg import getInitialUcUe
		self.z = getZArray()
		r = initcond.getInitialPresRhoVel(self.z)	
		self.pres = r["pres"]
		self.rho = r["rho"]
		self.vel = r["vel"]
		r = getInitialUcUe(self.rho, self.vel, self.pres)
		self.uc = r['uc']
		self.ue = r['ue']
		#flux arrays
		self.fm = None
		self.fc = None
		self.fe = None
		from notifier_params import notifierType
		self.notifier = getNotifier(notifierType, self.z, ["pres", "rho", "vel"], self.getInitialValues())
		#If I don't want analitical solution plotted for velocity:
		#self.notifier = getNotifier(notifierType, self.z, ["rho", "pres", "vel"], [[self.rho, self.rho], [self.pres, self.pres], self.vel)
		self.additionalInit()		
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
			#check if dt is 0 -> because  pres/rho might have become negative see  getTimestep in alg.py
			if dt==0:
				print("STOP")
				import time
				time.sleep(5)
				break
			time+=dt
			nstep+=1
			#recalculate u at next step - time	
			#print("Time %4.3f" % time)
			result = recalculateU(self.rho, self.uc, self.ue, self.fm, self.fc, self.fe, dt)
			self.rho = result["rho"]
			self.uc = result["uc"]
			self.ue = result["ue"]
			#recalculate velocity and pression
			r = recalculateVelPres(self.rho, self.uc, self.ue)
			self.vel = r["vel"]	
			self.pres = r["pres"]
			#print("NSTEP %d" % nstep)
			self.printVars(time)
			#splitted updateValues in updateValuesModel and updateValuesNotifier because I want to catch stop condition computed in updateValuesModel in each step
			self.updateValuesModel(dt, time)	
			from notifier_params import nstepsPlot
			if(nstep % nstepsPlot == 0):
				#print("upd")
				self.updateValuesNotifier(dt, time)
				self.notifier.afterUpdateValues(time)
		self.notifier.finish()
