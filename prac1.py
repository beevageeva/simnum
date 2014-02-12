import os, math
import numpy as np
import matplotlib.pyplot as plt
import scipy.misc 
from matplotlib.widgets import Slider



class DiscreteSlider(Slider):
    """A matplotlib slider widget with discrete steps."""

    def set_val(self, val):
        discrete_val = int(val)
        # We can't just call Slider.set_val(self, discrete_val), because this 
        # will prevent the slider from updating properly (it will get stuck at
        # the first step and not "slide"). Instead, we'll keep track of the
        # the continuous value as self.val and pass in the discrete value to
        # everything else.
        xy = self.poly.xy
        xy[2] = discrete_val, 1
        xy[3] = discrete_val, 0
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt % discrete_val)
        if self.drawon: 
            self.ax.figure.canvas.draw()
        self.val = val
        if not self.eventson: 
            return
        for cid, func in self.observers.items():
            func(discrete_val)







a = -4
b = 8


def getDx(numberPoints):
	return float(b-a) / numPoints

def sech(x):
	return 1.0/np.cosh(x)


def ad(x):
	return 0.25 * sech(5 + 0.5 * x) * (2 * sech(x) ** 2 * np.sin(x) + (2 * np.cos(x) - np.sin(x) * np.tanh(5 + 0.5*x))*(1+np.tanh(x)) ) 


def func(x):
	return 0.5 * np.sin(x) * (np.tanh(x) + 1) * 1.0 / (np.cosh(0.5 * (x+10)))



def maxDif(numPoints, toPlot= False):
	x = np.linspace(a, b , numPoints)
	dx = getDx(numPoints)
	if toPlot:
		xd = []
		yd = []
		yan = []
	maxErr = 0
	for i in range(0,len(x)-1):
		xm = 0.5 * (x[i+1] + x[i])	
		nr = (func(x[i+1]) - func(x[i]))/dx
		ar = ad(xm)
		if(toPlot):
			xd.append(xm)
			yd.append(nr)
			yan.append(ar)
		err = abs(nr - ar)
		if(err>maxErr):
			maxErr = err
	if(toPlot):
		#print("maxDif")
		#print(xd)
		return [xd,yd,yan,maxErr]
	else:
		return maxErr

def maxDif2(numPoints, toPlot=False):
	dx = getDx(numPoints)
	maxErr = 0
	x = np.linspace(a, b , numPoints)
	if toPlot:
		xd = []
		yd = []
		yan = []
	for i in range(2,len(x)-2):
		nr = (8*(func(x[i+1])- func(x[i-1])) + (func(x[i-2]) - func(x[i+2])))/(12 * dx)
		ar = ad(x[i])
		if(toPlot):
			yd.append(nr)
			yan.append(ar)
		err = abs(nr - ar)
		if(err>maxErr):
			maxErr = err
	if(toPlot):
		return [x[2:numPoints -2],yd,yan,maxErr]
	else:
		return maxErr

#CONFIG
baseLog = 10
derivFunc = maxDif
#ENDCONFIG


#ax1 = plt.subplot(221)
#ax2 = plt.subplot(223)
#ax3 = plt.subplot(122)
ax1 = plt.subplot(221)
ax2 = plt.subplot(223)
ax3 = plt.subplot(222)
ax4 = plt.subplot(224)
ax1.set_xlabel("x")
ax2.set_xlabel("x")
ax1.set_ylabel("y(numeric)")
ax2.set_ylabel("y(analitic)")
ax3.set_xlabel("numPoints(log base %d)" % baseLog)
ax3.set_ylabel("absMaxErr(log base %d)" % baseLog)
ax4.set_xlabel("dx")
ax4.set_ylabel("absMaxErr")

ax1.grid(True)
ax2.grid(True)
ax3.grid(True)

npArray = [16,32,64,128,256,512,1028]
#I don't know the base, I have to use log from math package(for each element)
#xerror=np.log(npArray, baseLog)
xerror = []
yerror=[]
dxe = []
ye = []
for numPoints in [16,32,64,128,256,512,1028]:
	xerror.append(math.log(numPoints, baseLog))
	dxe.append(getDx(numPoints))
	if(numPoints==64):
		res = derivFunc(numPoints, True)
		#print("x=")
		#print(res[0])
		ax1.plot(res[0], res[1],  marker='o', linestyle='-', color="r")
		ax2.plot(res[0], res[2],  marker='o', linestyle='-', color="r")
		yval = res[3]
	else:
		res = derivFunc(numPoints)
		yval = res
	yerror.append(math.log(yval, baseLog))
	ye.append(yval)
print("polyfit coef a(n), a(n-1), ..a(0) in polynom a(n) * x**n + .. a(1)x + a(0)")
print("for log NumPoints and log AbsMaxErr:  degree 1")
print(np.polyfit(xerror,yerror,1))


ax3.plot(xerror, yerror,  marker='o', linestyle='-', color="r")
ax4.plot(dxe, ye, "ro")
lpolyax4, = ax4.plot(dxe,ye, "g-")

def drawPlot4(m):
	c = np.polyfit(dxe,ye,m)
	print("polyfit coef a(n), a(n-1), ..a(0) in polynom a(n) * x**n + .. a(1)x + a(0)")
	print("for dx and AbsMaxErr:  degree %d" % m)
	print(c)
	polyCurve = np.poly1d(c)
	ty = polyCurve(dxe) #teoretic y
	lpolyax4.set_ydata(ty)
	ax4.relim()
	ax4.autoscale_view(True,True,True)
	plt.draw()


varHash = {'m':1} # I have to keep track of the current slider value in order not to update if value is not changing
	#I want a discrete slider so several values of the slider will be converted to the same integer value
	#I have to use the hash because m is changed in the listener function and m variable would be local to this function
	#I don't use a global variable m. In python 3 there is the nonlocal statement which causes the listed identifiers to refer to previously bound variables in the nearest enclosing scope.
drawPlot4(1)	
axSlider = plt.axes([0.25, 0.01, 0.65, 0.03], axisbg='white')
	#mSlider =  Slider(axSlider, 'Degree', 1, 23, valinit=1, valfmt='%d')#max degree 23
	#the difference in using Slider or DiscreteSlider is only in the visual update of the slider bar
	#if discrete it will update only with portions corresponding to 1, I still have to keep track of the value
	#because the actual value of slider.val will always be the float
mSlider =  DiscreteSlider(axSlider, 'Degree', 1, 23, valinit=1, valfmt='%d')#max degree 23

def sliderChanged(val):
		print("SLIDER CHANGED EVENT " + str(val))
		intVal = int(val)
		#I don't want an update if values are equal (as integers)
		if(intVal!=varHash["m"]):
			print("Integer Value of slider CHANGED set slider to %d" % intVal)
			varHash["m"] = intVal
			drawPlot4(intVal)
		else:
			print("BUT integer VALUE did NOT change")

mSlider.on_changed(sliderChanged) 	



plt.draw()
plt.show()



