import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import sys, os
from scipy.fftpack import fft,fftfreq#forFourierTransform



#saveImages = False
saveImages = True

def getXLimits(title):
	if(title == "velFFT" or title == "presFFT"):
		return {"minX":-80, "maxX":80}
	from constants import z0, zf
	return {"minX" : z0, "maxX" : zf}

def getYLimits(title):
	from soundwave_medium_params import mediumType
	if mediumType == "homog":
		medType = "homog"
	else:	
		from soundwave_medium_params import inhomogSubtype
		if inhomogSubtype == 1:
			medType = "inhomog1"
		else:
			medType = "inhomog2"
	if(title == "pres"):
		if medType == "homog":
			return { "maxY": 1.0005, "minY": 0.9995} #homog
		elif medType == "inhomog1":
			return { "maxY": 1.0006, "minY": 0.9995} #inhomog1
		elif medType == "inhomog2":
			return {  "maxY": 1.0008, "minY": 0.9992} #inhomog2
	elif(title == "vel"):
		if medType == "homog":
			return	{ "maxY": 0.00035, "minY": -0.00035} 
		elif medType == "inhomog1":
			return	{ "maxY": 0.0015, "minY": -0.0015} 
		elif medType == "inhomog2":
			return	 { "maxY": 0.0015, "minY": -0.0015}
	elif(title == "rho"):
		if medType == "homog":
			return { "maxY": 1.0004, "minY": 0.9996} 
		elif medType == "inhomog1":
			return { "maxY": 1.0004, "minY": 0}	 
		elif medType == "inhomog2":
			return { "maxY": 1.3, "minY": 0} 	 
	elif(title == "rhoCurve"):
		if medType == "homog":
			return { "maxY": 0.00025, "minY": -0.00025}		
		elif medType == "inhomog1":
			return	{ "maxY": 0.00035, "minY": -0.00035}	
		elif medType == "inhomog2":
			return	{ "maxY": 0.0020, "minY": -0.0020}
	#I should always relim y for pres fft beacuse of the central value in k = 0 (= 1 ,=  p00, = mean value of p)
	#there is no need for vel fft beacuse mean vel = 0
	elif(title == "presFFT"):
		if medType == "inhomog1":
			return	{ "maxY": 0.00015, "minY": 0}
		elif medType == "homog":
			return	{ "maxY": 0.000025, "minY": 0}

	return None


def relimAxis(ax, title, setLimits = False):
	ylim = getYLimits(title)
	xlim = getXLimits(title)
	if(xlim or ylim):
		ax.relim()
	if not xlim and not ylim:
		ax.autoscale_view(True,True,True)
	else:
		if setLimits:	
			if ylim:
				from matplotlib.ticker import FormatStrFormatter
				ax.set_ylim(ylim["minY"],ylim["maxY"])
				ax.yaxis.set_major_formatter(FormatStrFormatter('%f'))
			if xlim:
				ax.set_xlim(xlim["minX"],xlim["maxX"])
		if not xlim:
			ax.autoscale_view(True,True,False)
		else:
			ax.autoscale_view(True,False,True)


#ylim = {"pres":{ "maxY": 1.0005, "minY": 0.9995} , "vel" : { "maxY": 0.00035, "minY": -0.00035}, "rho":{ "maxY": 1.0004, "minY": 0.9996}, 'rhoCurve': { "maxY": 0.00025, "minY": -0.00025} , 'velFFT' : {"maxX" : 80, "minX" : -80} } 
#inhom1
#ylim = {"pres":{ "maxY": 1.0006, "minY": 0.9995} , "vel" : { "maxY": 0.0015, "minY": -0.0015}, "rho":{ "maxY": 1.0004, "minY": 0}, 'rhoCurve': { "maxY": 0.00035, "minY": -0.00035} , 'velFFT' : {"maxX" : 80, "minX" : -80} } 
#inhom2
#ylim = {"pres":{ "maxY": 1.0008, "minY": 0.9992} , "vel" : { "maxY": 0.0015, "minY": -0.0015}, "rho":{ "maxY": 1.3, "minY": 0}, 'rhoCurve': { "maxY": 0.0020, "minY": -0.0020}, 'velFFT' : {"maxX" : 80, "minX" : -80} } 
#ylim = None


def getRandomColor():
	from random import random
	red = random()
	blue = random()
	green = random()	
	return "#%02x%02x%02x" % (red*255, green*255, blue*255)

def getColorFromArray(array):
	return "#%02x%02x%02x" % (array[0][0]*255, array[0][1]*255, array[0][2]*255)



	


class VisualPlot:

	def __init__(self, z, titles, iniValues):
		if saveImages:
			from common import createFolder	
			self.dirname = createFolder("outImages")
		nplots = len(titles)
		fig, ax =  plt.subplots(nplots,1,True)
		#fig.set_size_inches(300,200)
		#fig.set_figwidth(300)
		#fig.set_figheight(200)
		if(nplots == 1):
			ax = [ax]
		self.figures = [fig]
		self.lines = {}
		self.axes = {}
		for i in range(0, nplots):
			title = titles[i]
			self.addAxis(z, ax[i], title, iniValues[i])
		self.plotTitle = ax[0].set_title("Time 0")
		wm = plt.get_current_fig_manager()
		wm.window.wm_geometry("1000x900+50+50")
		fig.subplots_adjust(right=0.8)
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


	def addAxis(self, z, ax, title, vals, xlabel="z", linestyle="-"):
		ax.set_xlabel(xlabel)
		ax.set_ylabel(title)
		ax.grid(True)
		self.axes[title] = ax
		shape = np.shape(vals)
		plotLegend = False
		#we can plot multiple graphs on the same axis : example numerical and analytical
		if(len(shape)==1):
			if(linestyle == '-'):
				l, = ax.plot(z, vals, lw=2, color='b')
			else:
				l, = ax.plot(z, vals,  markersize=3, linestyle="None", marker="o",  color='b')
			#l, = ax.plot(z, vals[i], lw=2, color='b', markersize=5, linestyle="-", marker="o")
			self.lines[title] = l
		elif(len(shape)==2):
			self.lines[title] = []
			for i in range(0, shape[0]):
				if(linestyle == '-'):
					if(len(z.shape) == 2):
						l, = ax.plot(z[i], vals[i], lw=2, color=getRandomColor(),  label="%d" % i)
					else:
						l, = ax.plot(z, vals[i], lw=2, color=getRandomColor(),  label="%d" % i)
				else:
					if(len(z.shape) == 2):
						l, = ax.plot(z[i], vals[i], color=getRandomColor(),  markersize=3, linestyle="None", marker="o",  label="%d" % i)
					else:
						l, = ax.plot(z, vals[i], color=getRandomColor(),  markersize=3, linestyle="None", marker="o",  label="%d" % i)
				plotLegend = True
				#l, = ax.plot(z, vals[i], lw=2, color=getRandomColor(), markersize=5, linestyle="-", marker="o", label="%d" % i)
				self.lines[title].append(l)
#		ylim = getYLimits(title)
#		xlim = getXLimits(title)
#		if not xlim and not ylim:
#			ax.relim()
#			ax.autoscale_view(True,True,True)
#		else:
#			if (ylim):
#				from matplotlib.ticker import FormatStrFormatter
#				ax.set_ylim(ylim["minY"],ylim["maxY"])
#				ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))
#			if xlim:
#				ax.set_xlim(xlim["minX"],xlim["maxX"])
		if(plotLegend):
			ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small')
		relimAxis(ax, title, True)	
		#ax.relim()
		#ax.autoscale_view(True,True,True)

	def markPoint(self, axTitle, pointName, value):
		#print("markPoint on axTitle = %s, pointName = %s, value = %E" % (axTitle, pointName, value))
		from common import testKeyInDict
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
		delta = 0.2
		#I have to make the following test because
		#sometimes (in the case of riemann problem and initial velocity  0 )
		#because zC is only marked once at the beginning  when velocity is 0 for all z
		#TODO I choose 1 but it might be too small
		if(maxValue == minValue):
			maxValue = minValue + 1
		self.markPoints[pointName] = self.axes[axTitle].vlines(value, minValue - delta, maxValue + delta, color=color, label=pointName)
		self.axes[axTitle].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small')
			

	def addGraph(self, z, title, vals, xlabel="z", linestyle="-"):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		self.addAxis(z, ax, title, vals, xlabel, linestyle)
		self.figures.append(fig)
		fig.subplots_adjust(right=0.8)
		plt.draw()
		plt.show(block=False)

	def plotAxisTwin(self, z, title, vals, newtitle):
		ax2 = self.axes[title].twinx()
		ax2.set_ylabel(newtitle)
		ax2.plot(z, vals, color=getRandomColor())
		plt.draw()
		plt.show(block=False)

	def plotAxis(self, z, title, vals, label=None):
		ax = self.axes[title]
		ax.plot(z, vals, color=getRandomColor(),label=label)
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small')
		plt.draw()
		plt.show(block=False)

	def updateValues(self, title, newValues, newZ = None, sortValues = False):
		#TODO newZ as tuple
		if not newZ is None:
			sortValues = True
		#print("updateValues %s" % title)
		shape = np.shape(newValues)
		#we can plot multiple graphs on the same axis : example numerical and analytical: see addAxis before!!
		if(len(shape)==1):
			if not newZ is None:
				if sortValues:
					argSort = np.argsort(newZ)
					newZP = newZ[argSort]	
					newValuesP = newValues[argSort]
				else:
					newZP = newZ
					newValuesP = newValues	
				self.lines[title].set_xdata(newZP)
			else:
				newValuesP = newValues


			self.lines[title].set_ydata(newValuesP)

			#print(" ".join(map(str, newValues)))
		elif(len(shape)==2):
			#print(" ".join(map(str, newValues[0])))
			nlines = shape[0]
			if(hasattr(self.lines[title], "__len__") and len(self.lines[title])==nlines):
				for i in range(0, nlines):
					if not newZ is None:
						if sortValues:
							argSort = np.argsort(newZ[i])
							newZi = newZ[i][argSort]	
						else:
							newZi = newZ[i]
						self.lines[title][i].set_xdata(newZi)
					if sortValues:
						newValuesi = newValues[i][argSort]
					else:
						newValuesi = newValues[i]
					self.lines[title][i].set_ydata(newValuesi)
		relimAxis(self.axes[title], title)	
		
		
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
		import time
		time.sleep(10)
		#pass







