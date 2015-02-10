import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import sys, os
from scipy.fftpack import fft,fftfreq#forFourierTransform



saveImages = False
#saveImages = True

#ylim = {"pres":{ "maxY": 1.0005, "minY": 0.9995} , "vel" : { "maxY": 0.00035, "minY": -0.00035}, "rho":{ "maxY": 1.0004, "minY": 0.9996}, 'rhoCurve': { "maxY": 0.00025, "minY": -0.00025}} 
#inhom1
#ylim = {"pres":{ "maxY": 1.0006, "minY": 0.9995} , "vel" : { "maxY": 0.0015, "minY": -0.0015}, "rho":{ "maxY": 1.0004, "minY": 0}, 'rhoCurve': { "maxY": 0.00035, "minY": -0.00035}} 
#inhom2
#ylim = {"pres":{ "maxY": 1.0008, "minY": 0.9992} , "vel" : { "maxY": 0.0015, "minY": -0.0015}, "rho":{ "maxY": 1.3, "minY": 0}, 'rhoCurve': { "maxY": 0.0020, "minY": -0.0020}} 
from constants import z0, zf
xlim = {"minX" : z0, "maxX" : zf}
ylim = None


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
		self.z = z
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
			self.addAxis(ax[i], title, iniValues[i])
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


	def addAxis(self, ax, title, vals):
		ax.set_xlabel("z")
		ax.set_ylabel(title)
		ax.grid(True)
		self.axes[title] = ax
		shape = np.shape(vals)
		plotLegend = False
		#we can plot multiple graphs on the same axis : example numerical and analytical
		if(len(shape)==1):
			l, = ax.plot(self.z, vals, lw=2, color='b')
			#l, = ax.plot(self.z, vals[i], lw=2, color='b', markersize=5, linestyle="-", marker="o")
			self.lines[title] = l
		elif(len(shape)==2):
			self.lines[title] = []
			for i in range(0, shape[0]):
				l, = ax.plot(self.z, vals[i], lw=2, color=getRandomColor(), label="%d" % i)
				plotLegend = True
				#l, = ax.plot(self.z, vals[i], lw=2, color=getRandomColor(), markersize=5, linestyle="-", marker="o", label="%d" % i)
				self.lines[title].append(l)
		if(ylim):
			from matplotlib.ticker import FormatStrFormatter
			ax.set_ylim(ylim[title]["minY"],ylim[title]["maxY"])
			ax.set_xlim(xlim["minX"],xlim["maxX"])
			ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))
		else:
			ax.relim()
			ax.autoscale_view(True,True,True)
		if(plotLegend):
			ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small')
		ax.relim()
		ax.autoscale_view(True,True,True)

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
			
	def fftplot(self, ax, vals, aFunc=None):
		numPoints = len(self.z)
		ax.grid(True)
		Y=fft(vals)/(numPoints)
		F=fftfreq(numPoints, self.z[1] - self.z[0])
		intlen = self.z[len(self.z)-1] - self.z[0]
		#print("in visual plot fft vel kc = %E" % (intlen * abs(F[np.argmax(abs(Y[1:]))+1])) )
		#ax.set_xlim(-330,330)
		ax.set_xlim(-130,130) #inhomog first
		#ax.set_xlim(-40,40)# second exp of inhom
		#ax.set_xticks(np.arange(-130, 131, 20))
		ax.set_xticks(np.arange(-130, 131, 20))
		#ax.set_xticks(np.arange(-70, 70, 5)) #second exp of inhom
		ax.plot(intlen * F,abs(Y), markersize=3, linestyle="-", marker="o")
		if not aFunc is None:
			for i in range(len(F)):
				print("f=%4.3f,aFunc=%4.3f, fftval(Y)=%4.3f" % (F[i], aFunc(F[i]), Y[i]))
			c = np.polyfit(aFunc(F), abs(Y) , 1)
			print("coef polyfit aFunc(F), Y degree 1")
			print(c)
			ax.plot(F,aFunc(F), markersize=3, linestyle="-", marker="o", color="r")
			


	def addFFTAxis(self, title, vals, aFunc=None):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		self.axes[title] = ax
		self.fftplot(ax, vals, aFunc)
		self.figures.append(fig)
		fig.subplots_adjust(right=0.8)
		plt.draw()
		plt.show(block=False)
		
	def updateFFTAxis(self, title, vals, aFunc=None):
		ax = self.axes[title]
		ax.cla()
		self.fftplot(ax,vals, aFunc)

	def addGraph(self, title, vals):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		self.addAxis(ax, title, vals)
		self.figures.append(fig)
		fig.subplots_adjust(right=0.8)
		plt.draw()
		plt.show(block=False)

	def plotAxisTwin(self, title, vals, newtitle):
		ax2 = self.axes[title].twinx()
		ax2.set_ylabel(newtitle)
		ax2.plot(self.z, vals, color=getRandomColor())
		plt.draw()
		plt.show(block=False)

	def plotAxis(self, title, vals, label=None):
		ax = self.axes[title]
		ax.plot(self.z, vals, color=getRandomColor(),label=label)
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small')
		plt.draw()
		plt.show(block=False)

	def updateValues(self, title, newValues):
		#print("updateValues %s" % title)
		shape = np.shape(newValues)
		#we can plot multiple graphs on the same axis : example numerical and analytical: see addAxis before!!
		if(len(shape)==1):
			self.lines[title].set_ydata(newValues)
			#print(" ".join(map(str, newValues)))
		elif(len(shape)==2):
			#print(" ".join(map(str, newValues[0])))
			nlines = shape[0]
			if(hasattr(self.lines[title], "__len__") and len(self.lines[title])==nlines):
				for i in range(0, nlines):
					self.lines[title][i].set_ydata(newValues[i])
		if(not ylim):
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
		import time
		time.sleep(10)
		#pass







