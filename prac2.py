import os, math, sys
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scipy.misc 
from matplotlib.widgets import Slider, Button


plot_ao = True
#plot_ao = False

#nslider = True
nslider = False

method = "forward"
#method = "centered"
#method="backward"  


xL=1.2
#ex 2.a
a = -1
#nint = 512
#ex 2.b
#a = 1
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

def getDx(numInt):
	return float(2.0 * xL) / numInt

def sech(x):
	return 1.0/np.cosh(x)

def getDt(numInt):
	return (0.98 / numInt) * abs(2.0 * xL / a)

def func_t0(x):
	return np.power(np.cos(6 * math.pi * x / 5) ,2) / np.cosh(5 * x**2)

#def plotFunc_t0(nint):
#	x = np.linspace(-xL, xL , nint + 1)
#	y = func_t0(x)
#	ax.plot(x, y,  marker='o', linestyle='-', color="r")
#	plt.draw()	
#	plt.show()	


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
		print("calc num sol(forward) solution for n = %d, nint = %d" % (n , nint))	
		x = np.linspace(-xL, xL , nint + 1)
		uAnt = func_t0(x)
		dx = getDx(nint)
		dt = getDt(nint)
		res = uAnt
		for j in range(0, n):
			res = []
			for i in range(0, nint):
				val = uAnt[i] - a * (uAnt[i+1] - uAnt[i] ) * dt / dx
				res.append(val)
			res.append(res[0])
			uAnt = res
		return res

elif method == "backward":
	#t_0=0
	#t = n * dt 
	def calcFunc_t(n, nint):
		print("calc num sol(backward) solution for n = %d, nint = %d" % (n , nint))	
		x = np.linspace(-xL, xL , nint + 1)
		uAnt = func_t0(x)
		dx = getDx(nint)
		dt = getDt(nint)
		res = uAnt
		for j in range(0, n):
			res = []
			for i in range(1, nint + 1):
				val = uAnt[i] - a * (uAnt[i] - uAnt[i - 1] ) * dt / dx
				res.append(val)
			res.insert(0, res[nint-1])
			uAnt = res
		return res

elif method == "centered":
	#t_0=0
	#t = n * dt 
	def calcFunc_t(n, nint):
		print("calc num sol(centered) solution for n = %d, nint = %d" % (n , nint))	
		x = np.linspace(-xL, xL , nint + 1)
		uAnt = func_t0(x)
		dx = getDx(nint)
		dt = getDt(nint)
		res = uAnt
		for j in range(0, n):
			res = []
			for i in range(1, nint):
				val = uAnt[i] - 0.5 * a * (uAnt[i + 1] - uAnt[i - 1] ) * dt / dx
				res.append(val)
			res.insert(0, res[nint-2])
			res.append(res[1])
			uAnt = res
		return res

else:
	print("method unknown: %s" % method)



def anFunc_t(n, nint):
	#print("anFunc_t: n= %d, nint=%d" % (n, nint))
	#t = n * dt
	#u(x,t) = u_0(x - a*t)
	print("calc an solution for n = %d, nint = %d" % (n , nint))	
	dt = getDt(nint)
	x = np.linspace(-xL, xL , nint + 1)
	#periodic boundary condition:
	#u(xL, t) = u(0, t) for all t EQUIV u(n*xL + x, t) = u(x,t)
	#we must haxe all x in [-xL , xL] when calculating func
	#ft0xarg = x - a * n * dt
	ft0xarg = []
	for xval in  x - a * n * dt:
		#TODO getPeriodicX	
		#ft0xarg.append(getPeriodicX(xval, xL))
		ft0xarg.append(getPeriodicX2(xval, xL))
	ft0xarg = np.array(ft0xarg)
	res = func_t0(ft0xarg)
	#print("newxarg")
	#print(ft0xarg)
	#print("funct0")
	#print(res)
	return res
	
	

def plotFunc_t(n, nint):
	lnf.set_ydata(calcFunc_t(n, nint))
	laf.set_ydata(anFunc_t(n, nint))
	ax.relim()
	ax.autoscale_view(True,True,True)
	plt.draw()	
	plt.show()	


print("display nint = %d" % nint)
ax = plt.subplot(111)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("numInt=%d"%nint) 
ax.grid(True)

x = np.linspace(-xL, xL , nint + 1)
lnf, = ax.plot(x, func_t0(x),  marker='o', linestyle='-', color="r")
laf, = ax.plot(x, func_t0(x),  marker='o', linestyle='--', color="g")
varHash = {'n':0} # I have to keep track of the current slider value in order not to update if value is not changing
	#I want a discrete slider so several values of the slider will be converted to the same integer value
	#I have to use the hash because m is changed in the listener function and m variable would be local to this function
	#I don't use a global variable m. In python 3 there is the nonlocal statement which causes the listed identifiers to refer to previously bound variables in the nearest enclosing scope.
axSlider = plt.axes([0.25, 0.01, 0.65, 0.03], axisbg='white')
print("dt=%4.10f" % getDt(nint))

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
	mSlider =  DiscreteSlider(axSlider, 'N', 0, 100, valinit=0, valfmt='%d')
	mSlider.on_changed(sliderChangedN) 	
else:	
	mSlider =  Slider(axSlider, 'Time', 0, 100, valinit=0, valfmt='%4.3f' )#max degree 23
	mSlider.on_changed(sliderChangedTime) 	

axObutt = plt.axes([0.93, 0.05, 0.05, 0.05])
obutt = Button(axObutt, 'AO')

def showApproxOrder(event):
	if(varHash["n"]==0):
		print("fot t = 0 functions should be identical")
		return
	baseLog = 10
	xe = []
	ye = []
	#for ni in [16, 32, 64, 128,256,512,1024]:
	#for ni in [256,512,1024,2048,4096]:
	for ni in [64,128,256,512,1024]:
		#step size is DEPENDENT of nint!!!!
		n = int(float(varHash["n"] * ni) / nint)
		print("Calculate approx order for n = %d,  ni=%d" % (n, ni))
		xe.append(math.log(ni, baseLog))
		y1 = calcFunc_t(n, ni)
		y2 = anFunc_t(n, ni)
		err = np.max(np.absolute(np.subtract(y1, y2)))	
		ye.append(math.log(err, baseLog))
		#plot function
		if plot_ao:
			ax.set_title("numInt=%d"%ni) 
			x = np.linspace(-xL, xL , ni + 1)
			lnf.set_xdata(x) 
			laf.set_xdata(x) 
			lnf.set_ydata(y1) 
			laf.set_ydata(y2) 
			ax.relim() 
			ax.autoscale_view(True,True,True) 
			plt.draw()  
			if (sys.version_info[0]==2):
				c = raw_input("press n to continue: ")
			else:
				c = input("press n to continue: ")
			while(c!="n"):
				print("You pressed >%s<" % c)
				if (sys.version_info[0]==2):
					c = raw_input("press n to continue: ")
				else:
					c = input("press n to continue: ")
		#end plot function
	#replot for display nint
	if plot_ao:
		n = varHash["n"]
		x = np.linspace(-xL, xL , nint + 1)
		lnf.set_xdata(x)
		laf.set_xdata(x)
		ax.set_title("numInt=%d"%nint)
		lnf.set_ydata(calcFunc_t(n, nint))
		laf.set_ydata(anFunc_t(n, nint))
		ax.relim()
		ax.autoscale_view(True,True,True)
	#end replot for display nint


	cpf = np.polyfit(xe,ye,1)
	print("t = dt * n (n = %d), fitting %4.10fx+%4.10f" % (n, cpf[0], cpf[1]))
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax1.set_xlabel("log(ni)")
	ax1.set_ylabel("log(maxerr)")
	ax1.set_title("baseLog=%d" % baseLog)
	ax1.grid(True)
	ax1.plot(xe, ye, marker='o', linestyle='-', color="r")
	plt.draw()
	plt.show()


obutt.on_clicked(showApproxOrder)


plt.draw()
plt.show()

#plotFunc_t(2,64)

