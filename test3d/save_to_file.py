import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import sys, os



def getStringFromFloatArray(arr):
	res = ""
	first = True
	for aval in arr:
		if first:
			first = False
		else:
			res += " " 
		res += "%E" % aval
	return res
	


class SaveToFile:


	def __init__(self, z, titles, iniValues):
		from common import createFolder
		self.dirname = createFolder("outFiles")
		self.z = z
		self.files = {}
		for i in range(0, len(titles)):
			self.addGraph(titles[i], iniValues[i])

	def afterInit(self):
		pass

	def addGraph(self, title, vals):
		filename = os.path.join(self.dirname, title)	
		self.files[title] = open(filename, "w")
		self.files[title].write(getStringFromFloatArray(self.z))
		self.files[title].write("\n")
		self.updateValues(title, vals)	

	def markPoint(self, title, pointName, value):
		self.files[title].write("Mark point %s with value %4.3f" % (pointName, value))
		self.files[title].write("\n")
		


	def updateValues(self, title, vals):
		shape = np.shape(vals)
		#we can plot multiple graphs on the same axis : example numerical and analytical
		if(len(shape)==1):
			self.files[title].write(getStringFromFloatArray(vals))
			self.files[title].write("\n")
		elif(len(shape)==2):
			for i in range(0, shape[0]):
				self.files[title].write(getStringFromFloatArray(vals[i]))
				self.files[title].write("\n")
		
	def afterUpdateValues(self, newTime):
		pass


	def finish(self):
		for title in self.files.keys():
			self.files[title].close()





