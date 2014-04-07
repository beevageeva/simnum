import matplotlib
from mpl_toolkits.mplot3d import Axes3D
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import sys, os

from notifier_params import fullscreenMainfigure

saveImages = False


	


class VisualPlot:



	def __init__(self, z, titles, iniValues):
		if saveImages:
			from common import createFolder	
			self.dirname = createFolder("outImages")
		#if for example velocity has 2 components I have to have 2 subplots
		self.z = z
		fig = plt.figure()
		self.figures = [fig]
		self.axes = {}
		n = len(titles)
		for i in range(0, len(titles)):
			vals = iniValues[i]
			title = titles[i]
			if(vals.ndim == 2):	
				#ax = fig.add_subplot(len(titles), 1, nplotted,  projection='3d')
				ax = plt.subplot2grid((n,2), (i,0), colspan=2, projection='3d')
				self.addAxis(ax, title, vals)
				self.axes[title] = ax
			elif(vals.ndim == 3):	
				#ax = fig.add_subplot(len(titles), 2, nplotted,  projection='3d')
				ax = plt.subplot2grid((n,2), (i,0), projection='3d')
				self.addAxis(ax, ("%s dim 0" % title), vals[:,:,0])
				
				#ax2 = fig.add_subplot(len(titles), 2, nplotted + 1,  projection='3d')
				ax2 = plt.subplot2grid((n,2), (i,1), projection='3d')
				self.addAxis(ax2, ("%s dim 1" % title), vals[:,:,1])
				self.axes[title] = [ax, ax2]
			else:
				print("Dim invalid %d" % vals.ndim)
				sys.exit(0)
		fax  =  self.axes[self.axes.keys()[0]]
		if(hasattr(fax, "__len__")):
			fax = fax[0]	
		self.plotTitle = fax.set_title("Time 0")
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
		ax.set_xlabel("z")
		ax.set_ylabel(title)
		ax.grid(True)
		#ax.plot_surface(self.z[0], self.z[1], vals, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
		#ax.plot_surface(self.z[0], self.z[1], vals)
		ax.plot_wireframe(self.z[0], self.z[1], vals)
		ax.view_init(0, 90)
		ax.relim()
		ax.autoscale_view(True,True,True)
	#		if fullscreenMainfigure:
	#			ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	#		else:
	#			ax.legend()
	

	def updateAxis(self, ax, vals):	
		ax.cla()
		#ax.plot_surface(self.z[0], self.z[1], vals, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
		#ax.plot_surface(self.z[0], self.z[1], vals)
		ax.plot_wireframe(self.z[0], self.z[1], vals)
		ax.relim()
		ax.autoscale_view(True,True,True)




	def updateValues(self, title, vals):
		ax = self.axes[title]
		if(vals.ndim == 2):
			self.updateAxis(ax, vals)
		else:
			self.updateAxis(ax[0], vals[:,:,0])
			self.updateAxis(ax[1], vals[:,:,1])

			
		
	def afterUpdateValues(self, newTime):
		#self.plotTitle.set_text("Time %4.4f" % newTime)
		print("Time %4.4f" % newTime)
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
		import time
		time.sleep(5)

	def finish(self):
		pass







