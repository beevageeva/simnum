import matplotlib
from mpl_toolkits.mplot3d import Axes3D
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import sys, os

from notifier_params import fullscreenMainfigure, projections

saveImages = False


	
#because methods are different in python3 and python2
def testKeyInDict(key, dictionary):
	import sys	
	if (sys.version_info[0]==2):
		return dictionary.has_key(key)
	else:
		return key in dictionary	


class VisualPlot:


	def addProjAxes(self, arrayToAppendAxes, title, vals, n , i, subplotNumber, colspan=False):
		if hasattr(self, "dim0ProjIndex"):
			plt.figure(2)
			if(colspan):
				ax = plt.subplot2grid((n,2), (i,subplotNumber), colspan=3)
			else:
				ax = plt.subplot2grid((n,2), (i,subplotNumber))
			if(vals.ndim == 3):
				values = vals[self.dim0ProjIndex, :,:]
			else:
				values = vals[self.dim0ProjIndex, :, :,subplotNumber]
			if not testKeyInDict(title, self.markPoints):
				self.markPoints[title]= {"dim0"}

			self.addAxisProj(ax, title, values)
			arrayToAppendAxes.append(ax)
		if hasattr(self, "dim1ProjIndex"):
			if(vals.ndim == 3):
				values = vals[:,self.dim1ProjIndex,:]
			else:
				values = vals[:,self.dim1ProjIndex,:, subplotNumber]
			plt.figure(3)
			if(colspan):
				ax = plt.subplot2grid((n,2), (i,subplotNumber), colspan=3)
			else:
				ax = plt.subplot2grid((n,2), (i,subplotNumber))
			self.addAxisProj(ax, title, values)
			arrayToAppendAxes.append(ax)
		if hasattr(self, "dim2ProjIndex"):
			if(vals.ndim == 3):
				values = vals[:,:,self.dim2ProjIndex]
			else:
				values = vals[:,:,self.dim2ProjIndex, subplotNumber]
			plt.figure(4)
			if(colspan):
				ax = plt.subplot2grid((n,2), (i,subplotNumber), colspan=3)
			else:
				ax = plt.subplot2grid((n,2), (i,subplotNumber))
			self.addAxisProj(ax, title, values)
			arrayToAppendAxes.append(ax)


	
	def addAxisProj(self, ax, title, vals):
		ax.set_xlabel("z")
		ax.set_ylabel(title)
		ax.grid(True)
		ax.imshow(vals)
		ax.relim()
		ax.autoscale_view(True,True,True)
	#		if fullscreenMainfigure:
	#			ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	#		else:
	#			ax.legend()


	def __init__(self, z, titles, iniValues):
		if saveImages:
			from common import createFolder	
			self.dirname = createFolder("outImages")
		#if for example velocity has 2 components I have to have 2 subplots
		self.z = z
		fig = plt.figure(1)
		self.figures = [fig]
		if projections:
			self.markPoints = {}
			from common import getZIndex0, getZIndex1, getZIndex2
			if testKeyInDict("dim0", projections):
				fig = plt.figure(2)
				self.dim0ProjIndex = getZIndex0(projections["dim0"])
				fig.suptitle("Dim0 z1=%4.3f" % projections["dim0"])
				self.figures.append(fig)
			if testKeyInDict("dim1", projections):
				fig = plt.figure(3)
				fig.suptitle("Dim1 z0=%4.3f" % projections["dim1"])
				self.dim1ProjIndex = getZIndex1(projections["dim1"])
				self.figures.append(fig)
			if testKeyInDict("dim2", projections):
				fig = plt.figure(4)
				fig.suptitle("Dim2 z0=%4.3f" % projections["dim2"])
				self.dim2ProjIndex = getZIndex2(projections["dim2"])
				self.figures.append(fig)
		self.axes = {}
		n = len(titles)
		for i in range(0, len(titles)):
			vals = iniValues[i]
			title = titles[i]
			plt.figure(1)
			if(vals.ndim == 3):	
				ax = plt.subplot2grid((n,2), (i,0), colspan=3, projection='3d')
				self.addAxis(ax, title, vals)
				self.axes[title] = [ax]
				self.addProjAxes(self.axes[title], title, vals, n, i, 0, True)

			elif(vals.ndim == 4):	
				ax = plt.subplot2grid((n,2), (i,0), projection='3d')
				self.addAxis(ax, ("%s dim 0" % title), vals[:,:,:,0])
				self.axes[title] = [[ax]]
				self.addProjAxes(self.axes[title][0], title, vals, n, i, 0)
				plt.figure(1)
				ax2 = plt.subplot2grid((n,2), (i,1), projection='3d')
				self.addAxis(ax2, ("%s dim 1" % title), vals[:,:,:,1])
				self.axes[title].append([ax2])
				self.addProjAxes(self.axes[title][1], title, vals, n, i, 1)
				plt.figure(1)
				ax2 = plt.subplot2grid((n,2), (i,2), projection='3d')
				self.addAxis(ax2, ("%s dim 2" % title), vals[:,:,:,2])
				self.axes[title].append([ax2])
				self.addProjAxes(self.axes[title][2], title, vals, n, i, 2)
			else:
				#print("Dim invalid %d" % vals.ndim)
				sys.exit(0)
		self.plotTitle = self.figures[0].suptitle("Time 0")
		plt.figure(1)
		wm = plt.get_current_fig_manager()
		if fullscreenMainfigure:
			#I think this only works with TkAgg backend
			wm.full_screen_toggle()
		else:
			wm.window.wm_geometry("800x900+50+50")
		plt.draw()
		plt.show(block=False)

	def afterInit(self):
		import time
		time.sleep(10)
		#save initial figures to files
		if saveImages:
			numFig = 0
			for fig in self.figures:
				os.mkdir(os.path.join(self.dirname, "Fig%d" % numFig))
				fig.savefig(os.path.join(self.dirname, "Fig%d" % numFig , "img000000.png"))
				numFig +=1

		
	def addAxis(self, ax, title, vals):
		ax.set_xlabel("x")
		ax.set_ylabel("y")
		ax.set_zlabel("z")
		ax.set_title(title)
		ax.grid(True)
		#ax.plot_surface(vals, cmap=cm.coolwarm,linewidth=0, antialiased=False)
		ax.plot_surface(vals[0], vals[1], vals[2])
		#ax.plot_surface(vals[0][0][:,0], vals[1][:,0][:,0], vals[2][:,0][0])
		#ax.view_init(0, 90)
		#ax.view_init(45, 45)
		ax.relim()
		ax.autoscale_view(True,True,True)
	#		if fullscreenMainfigure:
	#			ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	#		else:
	#			ax.legend()
	

	def updateAxis(self, ax, title, vals):	
		ax.cla()
		ax.set_xlabel("x")
		ax.set_ylabel("y")
		ax.set_zlabel("z")
		ax.set_title(title)
		#ax.plot_surface(vals, cmap=cm.coolwarm,linewidth=0, antialiased=False)
		ax.plot_surface(vals[0], vals[1], vals[2])
		ax.grid(True)
		ax.relim()
		ax.autoscale_view(True,True,True)

	def updateAxisProj(self, ax, title, vals):	
		ax.cla()
		ax.set_title(title)
		ax.imshow(vals)
		ax.grid(True)
		ax.relim()
		ax.autoscale_view(True,True,True)



	def updateProjAxis(self, axesArray, title, vals, index=None):
		if hasattr(self, "dim0ProjIndex"):
			ni = 2
			if(vals.ndim == 3):
				values = vals[self.dim0ProjIndex, :, :]
			else:
				values = vals[self.dim0ProjIndex, :, :, index]	
			#print("dim0")
			#print(" ".join(map(str, values)))
			self.updateAxisProj(axesArray[1], title, values)
		else:
			ni = 1
		if hasattr(self, "dim1ProjIndex"):
			if(vals.ndim == 3):
				values = vals[:,self.dim1ProjIndex,:]
			else:
				values = vals[:,self.dim1ProjIndex, :, index]
			self.updateAxisProj(axesArray[ni], title, values)
			ni+=1
		if hasattr(self, "dim2ProjIndex"):
			if(vals.ndim == 3):
				values = vals[:,:, self.dim1ProjIndex]
			else:
				values = vals[:,:,self.dim2ProjIndex, index]
			self.updateAxisProj(axesArray[ni], title, values)



	def updateValues(self, title, vals):
		#print("updateValues %s" % title)
		ax = self.axes[title]
		if(vals.ndim == 3):
			self.updateAxis(ax[0], title, vals)
			self.updateProjAxis(ax, title, vals)
		else:
			self.updateAxis(ax[0][0], ("%s dim0" % title), vals[:,:,:,0])
			self.updateProjAxis(ax[0], ("%s dim0" %title), vals, 0)
			self.updateAxis(ax[1][0], ("%s dim1" %title), vals[:,:,:,1])
			self.updateProjAxis(ax[1], ("%s dim1" %title), vals, 1)
			self.updateAxis(ax[2][0], ("%s dim2" %title), vals[:,:,:,2])
			self.updateProjAxis(ax[2], ("%s dim2" %title), vals, 2)

			
		
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






