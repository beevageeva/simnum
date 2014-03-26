import numpy as np
import sys
from visual_plot import VisualPlot


nstepsPlot = 1

nint = 64
x0 = -5.0
xf = 5.0
K = 10.0

#initial function
def u0(x):
	return 1 + np.exp(- (np.power(x, 2)/4))

#boundary conditions
def bc(array):
	array.insert(0, 1)
	array.append(1)	


class Model:
	
	def __init__(self):
		x = np.linspace(x0, xf, nint+1)
		self.u = u0(x)
		self.notifier = VisualPlot(x, ["u"], [self.u])
		self.notifier.afterInit()

	def recalcU(self, lu):
		r = []
		for i in range(1, len(self.u) - 1 ):
			val = self.u[i] + K * lu * (self.u[i+1] + self.u[i-1] -  2 * self.u[i])
			r.append(val)
		#boundary conditions
		bc(r)
		return r
		


	def mainLoop(self, timeEnd):
		time = 0.0
		nstep = 0
		dx = (xf - x0) / float(nint)
		fvn = 0.98
		dt = fvn * (dx ** 2)/(2.0 * K)
		print("dt=%E" % dt)
		lu = fvn / (2.0 * K)

		while(time<timeEnd):
			time+=dt
			self.u = self.recalcU(lu)
			if(nstep % nstepsPlot == 0):
				self.notifier.updateValues("u", self.u)
				self.notifier.afterUpdateValues(time)
		self.notifier.finish()
