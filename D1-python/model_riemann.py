import numpy as np
import sys
from constants import gamma
from base_model import BaseModel




class Model(BaseModel):
	
	def __init__(self):
		BaseModel.__init__(self)




	def updateValues(self, dt, time):
		self.notifier.updateValues("rho", self.rho)
		self.notifier.updateValues("pres", self.pres)
		self.notifier.updateValues("vel", self.vel)
		self.rwPoint = self.getNewPointLeft(self.rwPoint,dt)	
		self.shPoint = self.getNewPointRight(self.shPoint,dt)	
		#self.shPoint = self.getNewPointRight2(self.shPoint,dt)	
		self.notifier.markPoint("pres", "rwPres", self.rwPoint)
		self.notifier.markPoint("pres", "shPres", self.shPoint)
		self.notifier.markPoint("rho", "rwRho", self.rwPoint)
		self.notifier.markPoint("rho", "shRho", self.shPoint)
		self.notifier.markPoint("vel", "rwVel", self.rwPoint)
		self.notifier.markPoint("vel", "shVel", self.shPoint)


	def getInitialValues(self):
		return [self.pres, self.rho, self.vel]

	def additionalInit(self):
		from riemann_params import zC
		#for the moment the point name is the key in the hash - I should take into account the ax title and be able to repeat the pointNames in multiple axes - I have to have different names for different axes!
		self.notifier.markPoint("rho", "zCRho", zC)
		self.notifier.markPoint("pres", "zCPres", zC)
		self.notifier.markPoint("vel", "zCVel", zC)
		delta = 0.5#EMPIRICAL!!!
		self.rwPoint = zC - delta
		self.shPoint = zC + delta
		self.notifier.markPoint("pres", "rwPres", self.rwPoint)
		self.notifier.markPoint("pres", "shPres", self.shPoint)
		self.notifier.markPoint("rho", "rwRho", self.rwPoint)
		self.notifier.markPoint("rho", "shRho", self.shPoint)
		self.notifier.markPoint("vel", "rwVel", self.rwPoint)
		self.notifier.markPoint("vel", "shVel", self.shPoint)


	def getNewPointLeft(self, zval ,dt):
		from common import displacedPoint
		from initcond_riemann import velLeft
		return self.getNewPoint(zval, dt, velLeft, -1)	

	def getNewPointRight(self, zval ,dt):
		from common import displacedPoint
		from initcond_riemann import velRight
		return self.getNewPoint(zval, dt, velRight, 1)	

	def getNewPointRight2(self, zval, dt):
		from common import displacedPoint, getZIndex
		from initcond_riemann import getCsRight, velRight
		from math import sqrt
		v = velRight +  getCsRight()	
		newz = displacedPoint(zval, v, dt)
		return newz
