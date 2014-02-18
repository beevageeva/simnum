import os, math, sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.misc 
from matplotlib.widgets import Slider



a = -4
b = 8


def getDx(numInt):
	return float(b-a) / numInt

def sech(x):
	return 1.0/np.cosh(x)


def ad(x):
	return 0.25 * sech(5 + 0.5 * x) * (2 * sech(x) ** 2 * np.sin(x) + (2 * np.cos(x) - np.sin(x) * np.tanh(5 + 0.5*x))*(1+np.tanh(x)) ) 


def func(x):
	return 0.5 * np.sin(x) * (np.tanh(x) + 1) * 1.0 / (np.cosh(0.5 * (x+10)))



def maxDif(numInt, toPlot= False):
	x = np.linspace(a, b , numInt + 1)
	dx = getDx(numInt)
	if toPlot:
		xd = []
		yd = []
		yan = []
	maxErr = 0
	for i in range(0,len(x)-1):
		xm = 0.5 * (x[i+1] +x[i])	
		nr = (func(x[i+1]) - func(x[i]))/dx
		ar = ad(xm)
		if(toPlot):
			xd.append(xm)
			yd.append(nr)
			yan.append(ar)
		err = abs( np.float128(nr) -np.float128(ar))
		if(err>maxErr):
			maxErr = err
	print("MAX ERROR %4.20f" % maxErr)
	if(toPlot):
		return [xd,yd,yan,maxErr]
	else:
		return maxErr

def maxDif2(numInt, toPlot=False):
	dx = getDx(numInt)
	x = np.linspace(a, b , numInt + 1)
	maxErr = 0
	if(toPlot):
		yd = []
		yan = []
	for i in range(2,len(x)-2):
		nr = (8*(func(x[i+1])- func(x[i-1])) + (func(x[i-2]) - func(x[i+2])))/(12 * dx)
		ar = ad(x[i])
		if(toPlot):
			yd.append(np.float128(nr))
			yan.append(np.float128(ar))
		err = abs( np.float128(nr) -np.float128(ar))
		if(err>maxErr):
			maxErr = err
	print("ERROR %4.10f" % (maxErr) )
	if(toPlot):
		return [x[2:numInt -1],np.array(yd),np.array(yan), maxErr]
	else:
		return maxErr

def getValuesDer(numInt, ax=None):
	dx = getDx(numInt)
	x = np.linspace(a, b , numInt + 1)
	y = []
	for i in range(2,len(x)-2):
		nr = np.float128( (8*(func(x[i+1])- func(x[i-1])) + (func(x[i-2]) - func(x[i+2])))/(12 * dx))
		y.append(nr)
	#print("numInt=%d,x=" % numInt)
	#print(x[2:numInt -1])
	#print("y=")
	#print(y)
	if(ax):
		ax.plot(x[2:numInt -1], y,  marker='o', linestyle='-', color="r")
	return y




#CONFIG
baseLog = 10
derivFunc = maxDif
if(len(sys.argv)>1):
	if(sys.argv[1] == "2"):
		derivFunc = maxDif2
	if(len(sys.argv)>2):
		try:
			baseLog = int(sys.argv[2])
		except ValueError:
			print("second argument for baseLog is not number: %s using default 10" % sys.argv[2])
#ENDCONFIG



def makePlot():
	ax1 = plt.subplot(221)
	ax2 = plt.subplot(223)
	ax3 = plt.subplot(122)
	ax1.set_xlabel("x")
	ax2.set_xlabel("x")
	ax1.set_ylabel("y(numeric)")
	ax2.set_ylabel("y(analitic)")
	ax3.set_xlabel("numInt(log base %d)" % baseLog)
	ax3.set_ylabel("absMaxErr(log base %d)" % baseLog)
	ax1.set_title(derivFunc.func_name)
	
	ax1.grid(True)
	ax2.grid(True)
	ax3.grid(True)
	
	#npArray = [16,32,64,128,256,512,1024]
	niArray = [128,256,512,1028,512,1024, 2048, 4096, 8192]
	#I don't know the base, I have to use log from math package(for each element)
	#xerror=np.log(npArray, baseLog)
	xerror = []
	yerror=[]
	for numInt in niArray:
		xerror.append(math.log(numInt, baseLog))
		if(numInt==npArray[0]):
			res = derivFunc(numInt, True)
			#print("x=")
			#print(res[0])
			ax1.plot(res[0], res[1],  marker='o', linestyle='-', color="r")
			ax2.plot(res[0], res[2],  marker='o', linestyle='-', color="r")
			yval = res[3]
			testyval = np.max(np.absolute(np.subtract(res[1], res[2])))
			print("yval from result %4.10f , yval from maxNumpy %4.10f" % (yval, testyval))
		else:
			res = derivFunc(numInt, False)
			yval = res
				
		#yval = np.max(np.absolute(np.subtract(res[1], res[2])))
		print("yval = %4.10f , LOG = %4.10f" % (yval,math.log(yval, baseLog) ))
		yerror.append(math.log(yval, baseLog))
	clog = np.polyfit(xerror,yerror,1)
	print("yerror")
	print(yerror)
	print("polyfit coef a(n), a(n-1), ..a(0) in polynom a(n) * x**n + .. a(1)x + a(0)")
	print("for log NumPoints and log AbsMaxErr:  degree 1")
	polyCurve = np.poly1d(clog)
	yerror_curve = polyCurve(xerror)
	print(clog)
	
	ax3.plot(xerror, yerror, "ro")
	ax3.plot(xerror,yerror_curve, "g-")
	
	plt.draw()
	plt.show()


#m = log((f(2N) - f(N))/ (f(4N) - f(2N)))  / log(2)
#m = log((f(2N) - f(N))/ (f(4N) - f(2N)))  / log(2)
def calculateM(N, ax=None):
	ax1=ax2=ax3=None
	if(not ax is None):
		ax1 = ax[0]
		ax2 = ax[1]
		ax3 = ax[2]
	y1 = getValuesDer(N, ax1)
	y2 = getValuesDer(2 * N, ax2)
	y3 = getValuesDer(4*N, ax3)
	#print("f(N)")
	#print(y1)
	#print("f(2N)")
	#print(y2)
	#print("f(4N)")
	#print(y3)
	
	marray = []
	for i in range(0,len(y1)):
		t1 = np.float128(y2[2+2*i] - y1[i])
		t2 = np.float128(y3[6+4*i] - y2[2+2*i])
		if(t2==0):
			mval = 0
		else:
			mval =abs(t1/t2)
			#TODO append only if defined?
			marray.append(mval)
		#print("i=%d,f2N=%4.10f,fN=%4.10f,f4N=%4.10f, mval=%4.10f" % (i,y2[2+2*i], y1[i], y3[6+4*i], mval))
		#marray.append(mval)
	marrayMax = np.max(marray)
	indMax =  np.argmax(marray)
	print("log2 of Max value is %4.10f in position %d" % (math.log( np.max(marray), 2), indMax))
		
	print("log2 of Mean value is %4.10f" % (math.log( np.mean(marray), 2)))
	x = np.linspace(a, b , N + 1)[2:N -1]
	print("MAX: x=%4.10f, func(x)=%4.10f" % (x[indMax], func(x[indMax])))	
	print("ARGSORT: ")
	argIndS  =  np.argsort(marray)[::-1]
	print(argIndS)
	print("XARGSORT")
	print(x[argIndS])
	print("first 6")
	for i in range(0,6):
		print(x[argIndS][i])
#	if(not ax is None):
#		ax1.vlines(x[argIndS][0:5])
	
	

def makePlot2():
	ax1 = plt.subplot(311)
	ax2 = plt.subplot(312)
	ax3 = plt.subplot(313)
	ax1.set_xlabel("x")
	ax2.set_xlabel("x")
	ax3.set_xlabel("x")
	ax1.set_ylabel("y(N)")
	ax2.set_ylabel("y(2N)")
	ax3.set_xlabel("y(3N)")
	
	ax1.grid(True)
	ax2.grid(True)
	ax3.grid(True)
	#for i in  [8,16,32]:
	for i in  [16,32,128,256,512,1028,512,1024, 2048]:
		print("N=%d" % i)
		if(i == 16):
			calculateM(i, [ax1, ax2, ax3])
		else:
			calculateM(i)
	plt.draw()	
	plt.show()	
	


makePlot2()

