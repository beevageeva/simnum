import os, math
import numpy as np
import matplotlib.pyplot as plt
import scipy.misc 


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
baseLog = 2
derivFunc = maxDif
#ENDCONFIG


ax1 = plt.subplot(221)
ax2 = plt.subplot(223)
ax3 = plt.subplot(122)
ax1.set_xlabel("x")
ax2.set_xlabel("x")
ax1.set_ylabel("y(numeric)")
ax2.set_ylabel("y(analitic)")
ax3.set_xlabel("numPoints(log base %d)" % baseLog)
ax3.set_ylabel("absMaxErr(log base %d)" % baseLog)

ax1.grid(True)
ax2.grid(True)
ax3.grid(True)

npArray = [16,32,64,128,256,512,1028]
#I don't know the base, I have to use log from math package(for each element)
#xerror=np.log(npArray, baseLog)
xerror = []
yerror=[]
for numPoints in [16,32,64,128,256,512,1028]:
	xerror.append(math.log(numPoints, baseLog))
	if(numPoints==64):
		res = derivFunc(numPoints, True)
		#print("x=")
		#print(res[0])
		ax1.plot(res[0], res[1],  marker='o', linestyle='-', color="r")
		ax2.plot(res[0], res[2],  marker='o', linestyle='-', color="r")
		yerror.append(res[3])
	else:
		res = derivFunc(numPoints)
		yerror.append(res)
print("polyfit coef a(n), a(n-1), ..a(0) in polynom a(n) * x**n + .. a(1)x + a(0)")
print("degree 1")
print(np.polyfit(xerror,yerror,1))

ax3.plot(xerror, yerror,  marker='o', linestyle='-', color="r")

plt.draw()
plt.show()



