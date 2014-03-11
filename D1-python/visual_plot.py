import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import sys


def getRandomColor():
	from random import random
	red = random()
	blue = random()
	green = random()	
	return "#%02x%02x%02x" % (red*255, green*255, blue*255)

#because methods are different in python3 and python2
def testKeyInDict(key, dictionary):
	import sys	
	if (sys.version_info[0]==2):
		return dictionary.has_key(key)
	else:
		return key in dictionary	

	


class VisualPlot:

	def __init__(self, z, titles, iniValues):
		nplots = len(titles)
		self.z = z
		fig, ax =  plt.subplots(nplots,1,True)
		self.figures = [fig]
		self.lines = {}
		self.axes = {}
		for i in range(0, nplots):
			title = titles[i]
			self.addAxis(ax[i], title, iniValues[i])
		self.plotTitle = ax[0].set_title("Time 0")
		plt.draw()
		plt.show(block=False)

	def addAxis(self, ax, title, vals):
		ax.set_xlabel("z")
		ax.set_ylabel(title)
		ax.grid(True)
		self.axes[title] = ax
		shape = np.shape(vals)
		#we can plot multiple graphs on the same axis : example numerical and analytical
		if(len(shape)==1):
			l, = ax.plot(self.z, vals, lw=2, color='b')
			self.lines[title] = l
		elif(len(shape)==2):
			self.lines[title] = []
			for i in range(0, shape[0]):
				l, = ax.plot(self.z, vals[i], lw=2, color=getRandomColor(), label="%d" % i)
				self.lines[title].append(l)
			ax.legend()
		

	def markPoint(self, axTitle, pointName, value):
		if not hasattr(self, 'markPoints'):
			self.markPoints = {}
		if(testKeyInDict(pointName, self.markPoints)):
			self.markPoints[pointName].remove()
			del  self.markPoints[pointName]
		l = self.lines[axTitle]
		if hasattr(l, '__len__'):
			l = l[0]
		yvals = l.get_ydata()
		self.markPoints[pointName] = self.axes[axTitle].vlines(value, np.min(yvals), np.max(yvals), 'b')	
		

	def addGraph(self, title, vals):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		self.addAxis(ax, title, vals)
		self.figures.append(fig)
		plt.draw()
		plt.show(block=False)

	def updateValues(self, title, newValues):
		shape = np.shape(newValues)
		#we can plot multiple graphs on the same axis : example numerical and analytical: see addAxis before!!
		if(len(shape)==1):
			self.lines[title].set_ydata(newValues)
		elif(len(shape)==2):
			nlines = shape[0]
			if(hasattr(self.lines[title], "__len__") and len(self.lines[title])==nlines):
				for i in range(0, nlines):
					self.lines[title][i].set_ydata(newValues[i])
		self.axes[title].relim()
		self.axes[title].autoscale_view(True,True,True)
		
	def afterUpdateValues(self, newTime):
		self.plotTitle.set_text("Time %4.3f" % newTime)
		for fig in self.figures:
			fig.canvas.draw()
		#import time
		#time.sleep(1)








