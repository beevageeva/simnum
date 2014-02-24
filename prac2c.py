import os, math, sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.misc 
from matplotlib.widgets import Slider, Button



#NMAX=100
NMAX=300
TIMEMAX=20

#nslider = True
nslider = False

#method = "forward"
#method = "centered"  #this is unstable with function a(x) not constant
method="backward"  #this works


xL=1.2

#nint=64
nint=256
#nint=1024    
#nint=64 and method=centered is not working see Courant-Friedrichs-Lewy (CFL) condition dt <= dx / a 
#http://www.physics.udel.edu/~jim/PHYS460_660_13S/Advection/advection.htm





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

def a(x):
	a0 = 0.2
	b0 = 6
	return a0 * (1 + b0 * np.power(np.cos((4 * math.pi * x)/(9 * xL )),2))

def getDx(numInt):
	return float(2.0 * xL) / numInt

def sech(x):
	return 1.0/np.cosh(x)

def getDt(numInt):
	x = np.linspace(-xL, xL , numInt + 1)
	return np.min((0.98 / numInt) * abs(2.0 * xL / a(x)))

def func_t0(x):
	return np.power(np.cos(6 * math.pi * x / 5) ,2) / np.cosh(5 * x**2)

def plotFunc_t0(nint):
#	x = np.linspace(-xL, xL , nint + 1)
#	y = func_t0(x)
	ax.plot(x, func_t0(x),  marker='o', linestyle='-', color="r")
	plt.draw()	
	plt.show()	


def plotFunc_a():
	ax.plot(x, a(x),  marker='o', linestyle='-', color="r")
	plt.draw()	
	plt.show()	


#return x in [-a,a] assuming period 2a
def getPeriodicX(xval, a):
	while (xval < -a):
		xval+=2*a
	while (xval > a):
		xval-=2*a
	return xval

def getPeriodicXInt(xval, a, b):
	p = b - a
	while (xval < a):
		xval+=p
	while (xval > b):
		xval-=p
	return xval

def getPeriodicX2(xval, a):
	k = int(xval / (2 * a))
	res =  xval - 2 * a * k
	if(res>a):
		res -= 2*a
	if(res<-a):
		res += 2*a
	return res


def getPeriodicX2Int(xval, a, b):
	p = b - a
	k = int((xval-a)/p)
	res = xval - k * p
	if(res < a):
		res+=p
	if(res > a):
		res-=p
	return res


##t_0=0
##t = n * dt 
#def calcFunc_tRecursive(n, nint):
#	#print("calcFunc_t: n= %d, nint=%d" % (n, nint))
#	x = np.linspace(-xL, xL , nint + 1)
#	if(n==0):
#		return func_t0(x)
#	dx = getDx(nint)
#	dt = getDt(nint)
#	uAnt = calcFunc_t(n-1, nint)
#	res = []
#	for i in range(0, nint):
#		val = uAnt[i] - a * (uAnt[i+1] - uAnt[i] ) * dt / dx
#		res.append(val)
#	res.append(res[0])
#	return res



if method == "forward":
	#t_0=0
	#t = n * dt 
	def calcFunc_t(n, nint):
		#print("calcFunc_t: n= %d, nint=%d" % (n, nint))
		x = np.linspace(-xL, xL , nint + 1)
		uAnt = func_t0(x)
		dx = getDx(nint)
		dt = getDt(nint)
		res = uAnt
		for j in range(0, n):
			res = []
			for i in range(0, nint):
				val = uAnt[i] - a(x[i]) * (uAnt[i+1] - uAnt[i] ) * dt / dx
				res.append(val)
			res.append(res[0])
			uAnt = res
		return res

elif method == "backward":
	#t_0=0
	#t = n * dt 
	def calcFunc_t(n, nint):
		#print("calcFunc_t: n= %d, nint=%d" % (n, nint))
		x = np.linspace(-xL, xL , nint + 1)
		uAnt = func_t0(x)
		dx = getDx(nint)
		dt = getDt(nint)
		res = uAnt
		for j in range(0, n):
			res = []
			for i in range(1, nint + 1):
				val = uAnt[i] - a(x[i]) * (uAnt[i] - uAnt[i - 1] ) * dt / dx
				res.append(val)
			res.insert(0, res[nint-1])
			uAnt = res
		return res

elif method == "centered":
	#t_0=0
	#t = n * dt 
	def calcFunc_t(n, nint):
		#print("calcFunc_t: n= %d, nint=%d" % (n, nint))
		x = np.linspace(-xL, xL , nint + 1)
		uAnt = func_t0(x)
		dx = getDx(nint)
		dt = getDt(nint)
		res = uAnt
		for j in range(0, n):
			res = []
			for i in range(1, nint):
				val = uAnt[i] - 0.5 * a(x[i]) * (uAnt[i + 1] - uAnt[i - 1] ) * dt / dx
				res.append(val)
			res.insert(0, res[nint-2])
			res.append(res[1])
			uAnt = res
		return res

else:
	print("method unknown: %s" % method)


	



ax = plt.subplot(111)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("numInt=%d"%nint)
ax.grid(True)
x = np.linspace(-xL, xL , nint + 1)


def main():
	lnf, = ax.plot(x, func_t0(x),  marker='o', linestyle='-', color="r")
	varHash = {'n':0} # I have to keep track of the current slider value in order not to update if value is not changing
		#I want a discrete slider so several values of the slider will be converted to the same integer value
		#I have to use the hash because m is changed in the listener function and m variable would be local to this function
		#I don't use a global variable m. In python 3 there is the nonlocal statement which causes the listed identifiers to refer to previously bound variables in the nearest enclosing scope.
	axSlider = plt.axes([0.25, 0.01, 0.65, 0.03], axisbg='white')
	print("dt=%4.10f" % getDt(nint))
	
	def plotFunc_t(n, nint):
		lnf.set_ydata(calcFunc_t(n, nint))
		ax.relim()
		ax.autoscale_view(True,True,True)
		plt.draw()	
		plt.show()	

	def sliderChangedIntVal(intVal):
			#I don't want an update if values are equal (as integers)
			if(intVal!=varHash["n"]):
				print("Integer Value of slider CHANGED set slider to %d" % intVal)
				varHash["n"] = intVal
				plotFunc_t(intVal, nint)
			else:
				print("BUT integer VALUE did NOT change")
	
	
	def sliderChangedTime(val):
			print("SLIDER TIME CHANGED EVENT " + str(val))
			intVal = int(val / getDt(nint))
			sliderChangedIntVal(intVal)
	
	def sliderChangedN(val):
			print("SLIDER N CHANGED EVENT " + str(val))
			intVal = int(val)
			sliderChangedIntVal(intVal)
	
	
	if(nslider):
		mSlider =  DiscreteSlider(axSlider, 'N', 0, NMAX, valinit=0, valfmt='%d')
		mSlider.on_changed(sliderChangedN) 	
	else:	
		mSlider =  Slider(axSlider, 'Time', 0, TIMEMAX, valinit=0, valfmt='%4.3f' )#max degree 23
		mSlider.on_changed(sliderChangedTime) 	
	
	
	plt.draw()
	plt.show()

#plotFunc_t0(128)
main()
#plotFunc_a()

