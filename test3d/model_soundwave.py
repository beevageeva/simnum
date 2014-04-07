import numpy as np
import sys
from constants import gamma
from base_model import BaseModel




class Model(BaseModel):
	
	def __init__(self):
		BaseModel.__init__(self)



	def updateValuesModel(self, dt, time):
		pass


	def updateValuesNotifier(self, dt, time):
		self.notifier.updateValues("rho", self.rho)
		self.notifier.updateValues("pres", self.pres)
		self.notifier.updateValues("vel", self.vel)


	def additionalInit(self):
		pass


