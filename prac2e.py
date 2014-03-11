import os, math, sys
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scipy.misc 
from matplotlib.widgets import Slider, Button



#NMAX=100
NMAX=1000
TIMEMAX=100

plot_ao = True
plot_int = True


#nslider = True
nslider = False

#method="backward"  #this works


xL=1.2
#xL=30

#nint=64
#nint=256
nint=1024 
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



def getDx(numInt):
	return float(2.0 * xL) / numInt

def sech(x):
	return 1.0/np.cosh(x)

def getDt(numInt):
	x = np.linspace(-xL, xL , numInt + 1)
	return np.min((0.98 / numInt) * abs(2.0 * xL / func_t0(x)))

def func_t0(x):
	a = 0.05
	b = 0.65
	c = 0.25 
	#http://www.uwyo.edu/llee/papers/jcp_short_note_revision_v4.pdf
	#a = 0.5
	#b = 15
	#c = 0.01 
	return a * (np.tanh((x+b)/c) - np.tanh((x-b)/c)) 


def getTimeDisc(numInt):
	x = np.linspace(-xL, xL , numInt + 1)
	dx=	float(2.0 * xL) / numInt
	t = {}
	for i in range(0,len(x)-1):
		xm = 0.5 * (x[i+1] +x[i])	
		nr = (func_t0(x[i+1]) - func_t0(x[i]))/dx
		if(nr<0):
			t[-1.0 / nr] = xm
	return t
	



#return x in [-a,a] assuming period 2xL
def getPeriodicX(xval, xL):
	while (xval < -xL):
		xval+=2*xL
	while (xval > xL):
		xval-=2*xL
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



ax = plt.subplot(111)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("numInt=%d"%nint)
ax.grid(True)
x = np.linspace(-xL, xL , nint + 1)


def main():
	lnf, = ax.plot(x, func_t0(x),  marker='o', linestyle='-', color="r")
	#laf, = ax.plot(x, func_t0(x),  marker='o', linestyle='--', color="g")
	varHash = {'n':0} # I have to keep track of the current slider value in order not to update if value is not changing
		#I want a discrete slider so several values of the slider will be converted to the same integer value
		#I have to use the hash because m is changed in the listener function and m variable would be local to this function
		#I don't use a global variable m. In python 3 there is the nonlocal statement which causes the listed identifiers to refer to previously bound variables in the nearest enclosing scope.
	axSlider = plt.axes([0.25, 0.01, 0.65, 0.03], axisbg='white')
	print("dt=%4.10f" % getDt(nint))
	timeDiscOne = getTimeDisc(nint)

	#BACKWARD METHOD
	#t_0=0
	#t = n * dt 
	def calcFunc_t(n, nint, timeDisc,plot_int = False):
		print("calcFunc_t: n= %d, nint=%d" % (n, nint))
		x = np.linspace(-xL, xL , nint + 1)
		if(plot_int):
			markPoint = 0    #the maximum
			stepInterval = 10
			lmark = ax.vlines(markPoint, 0, 1, 'b')
			lmark2 = None
		uAnt = func_t0(x)
		dx = getDx(nint)
		dt = getDt(nint)
		res = uAnt
		for j in range(0, n):
			res = []
			#recalculate discontinuities:
			for k in timeDisc.iterkeys():
				v = timeDisc[k]
				timeDisc[k] = getPeriodicX2(v + uAnt[int((v+xL)/dx)] * dt, xL)
			v = None
			if (sys.version_info[0]==2):
				k = min(timeDisc.iterkeys())
			else:
				k = min(timeDisc.keys())
			if(k<=j*dt):
				v = timeDisc[k]
				print("DISC : time = %4.3f, x = %4.3f" % (k, v))
				del timeDisc[k]

			for i in range(1, nint + 1):
				val = uAnt[i] - uAnt[i] * (uAnt[i] - uAnt[i - 1] ) * dt / dx
				res.append(val)
			res.insert(0, res[nint-1])
			uAnt = res
			if(plot_int):
				markPoint = getPeriodicX2(markPoint + uAnt[int((markPoint+xL)/dx)] * dt, xL)
				if(j%stepInterval==0):
					#print("calcFunc_t INT: step= %d" % (j))
					ax.set_title("Time %4.3f" % ((j+1) * dt))
					lnf.set_ydata(res)
					#redraw mark point vlines

					lmark.remove()
					del lmark
					lmark = ax.vlines(markPoint, 0, 1, 'b')

					if lmark2:
						lmark2.remove()
						del lmark2
						lmark2 = None
					if(v):
						print("DISC-- : v = %4.3f" % v)
						lmark2 = ax.vlines(v, 0, 1, 'k', lw=5)


					ax.relim()
					ax.autoscale_view(True,True,True)
					plt.draw()
		if(plot_int):
			lmark.remove()
			del lmark
			if lmark2:
				lmark2.remove()
				del lmark2
					
		return res
	
	
	
	def plotFunc_t(n, nint):
		if(plot_int):
			calcFunc_t(n, nint,timeDiscOne.copy(), True)
		else:
			lnf.set_ydata(calcFunc_t(n, nint, timeDiscOne.copy()))
			#laf.set_ydata(anFunc_t(n, nint))
			ax.relim()
			ax.autoscale_view(True,True,True)
			plt.draw()	

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

main()

