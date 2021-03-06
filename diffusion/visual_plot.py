import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import sys, os


fullscreenMainfigure = True

saveImages = False

def getRandomColor():
	from random import random
	red = random()
	blue = random()
	green = random()	
	return "#%02x%02x%02x" % (red*255, green*255, blue*255)

def getColorFromArray(array):
	return "#%02x%02x%02x" % (array[0][0]*255, array[0][1]*255, array[0][2]*255)



#because methods are different in python3 and python2
def testKeyInDict(key, dictionary):
	import sys	
	if (sys.version_info[0]==2):
		return dictionary.has_key(key)
	else:
		return key in dictionary	


	


class VisualPlot:

	def __init__(self, z, titles, iniValues):
		if saveImages:
			from common import createFolder	
			self.dirname = createFolder("outImages")
		nplots = len(titles)
		self.z = z
		fig, ax =  plt.subplots(nplots,1,True)
		if(nplots == 1):
			ax = [ax]
		self.figures = [fig]
		self.lines = {}
		self.axes = {}
		for i in range(0, nplots):
			title = titles[i]
			self.addAxis(ax[i], title, iniValues[i])
		self.plotTitle = ax[0].set_title("Time 0")
		if fullscreenMainfigure:
			#I think this only works with TkAgg backend
			wm = plt.get_current_fig_manager()
			wm.full_screen_toggle()
		plt.draw()
		plt.show(block=False)

	def afterInit(self):
		#import time
		#time.sleep(5)
		#save initial figures to files
		if saveImages:
			numFig = 0
			for fig in self.figures:
				os.mkdir(os.path.join(self.dirname, "Fig%d" % numFig))
				fig.savefig(os.path.join(self.dirname, "Fig%d" % numFig , "img000000.png"))
				numFig +=1


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
			ax.relim()
			ax.autoscale_view(True,True,True)
			if fullscreenMainfigure:
				ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
			else:
				ax.legend()
		

	def markPoint(self, axTitle, pointName, value):
		if not hasattr(self, 'markPoints'):
			self.markPoints = {}
		if(testKeyInDict(pointName, self.markPoints)):
			self.markPoints[pointName].remove()
			#keep color
			color = getColorFromArray(self.markPoints[pointName].get_color())
			del  self.markPoints[pointName]
		else:
			#generate random color
			color = getRandomColor()
		l = self.lines[axTitle]
		if hasattr(l, '__len__'):
			l = l[0]
		yvals = l.get_ydata()
		minValue =  np.min(yvals)
		maxValue =  np.max(yvals)
		#I have to make the following test because
		#sometimes (in the case of riemann problem and initial velocity  0 )
		#because zC is only marked once at the beginning  when velocity is 0 for all z
		#TODO I choose 1 but it might be too small
		if(maxValue == minValue):
			maxValue = minValue + 1
		self.markPoints[pointName] = self.axes[axTitle].vlines(value, minValue, maxValue, color=color, label=pointName)
		if fullscreenMainfigure:
			self.axes[axTitle].legend(loc='center left', bbox_to_anchor=(1, 0.5))
		else:
			self.axes[axTitle].legend()
			
		

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
		self.plotTitle.set_text("Time %4.4f" % newTime)
		numFig = 0
		for fig in self.figures:
			fig.canvas.draw()
			if saveImages:
				#make name compatible with ffmpeg
				#ffmpeg -r 1 -i img%05d.png -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4
				#the above does not work
				#ffmpeg -r 1 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p out.mp4
				#convert HANGS!!
				#convert -antialias -delay 1x2 *.png mymovie.mp4
				imgname = "%4.4f" % newTime
				if(len(imgname) == 6):
					imgname = "0"+imgname
				imgname = imgname.replace(".", "")	
				fig.savefig(os.path.join(self.dirname, "Fig%d" % numFig , "img%s.png"%imgname))
			numFig +=1
		#import time
		#time.sleep(5)

	def finish(self):
		pass







