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


	def getInitialValues(self):
		return [self.pres, self.rho, self.vel]

	def additionalInit(self):
		from riemann_params import zC
		#for the moment the point name is the key in the hash - I should take into account the ax title and be able to repeat the pointNames in multiple axes - I have to have different names for different axes!
		self.notifier.markPoint("rho", "zCRho", zC)
		self.notifier.markPoint("pres", "zCPres", zC)
		self.notifier.markPoint("vel", "zCVel", zC)




