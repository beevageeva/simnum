import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft,fftfreq#forFourierTransform
from common import getZArray
from sound_wave_packet_params import k0, W, zc
from math import pi, exp, sqrt
from constants import z0, zf

from sound_wave import getSoundWavePacketFunction

def ftAnMat(k):
	return (np.exp(-(zc**2/W**2) + ((2.0*1j)*k0*pi*z0)/(z0 - zf))*(np.exp((W**2*((2.0*zc)/W**2 - (2.0*1j)*pi*(k + k0/(z0 - zf)))**2)/4.0) + np.exp(((4.0*1j)*k0*pi*z0)/(-z0 + zf) + (W**2*((2.0*zc)/W**2 - (2.0*1j)*pi*(k + k0/(-z0 + zf)))**2)/4.0))*sqrt(pi))/(2.0*sqrt(W**(-2)))

def ftAnNotW(k):
	z0C = complex(z0)
	zfC = complex(zf)
	WC = complex(W)
	k0C = complex(k0)
	t1 = np.exp(-(zc**2/WC**2) + ((complex(2.0)*1j)*k0C*pi*z0C)/(z0C - zfC))
	t2 = np.exp((WC**2*((complex(2.0)*zc)/WC**2 - (complex(2.0)*1j)*pi*(k + k0C/(z0C - zfC)))**2)/complex(4.0))
	t3 = np.exp(((complex(4.0)*1j)*k0C*pi*z0C)/(-z0C + zfC) + (WC**2*((complex(2.0)*zc)/WC**2 - (complex(2.0)*1j)*pi*(k + k0C/(-z0C + zfC)))**2)/complex(4.0))
	print("t1=%e, t2=%e, t3=%e" % (t1,t2,t3))
	return (t1*(t2 + t3)*sqrt(pi))/(2.0*sqrt(W**(-2)))

def ftAn(k):
	real = (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(z0 - zf))**2))/4)*sqrt(pi)*np.cos(2.0*pi*zc*(k + k0/(z0 - zf)))*np.cos((2.0*k0*pi*z0)/(z0 - zf)))/(2.0*sqrt(W**(-2))) + (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(-z0 + zf))**2))/4)*sqrt(pi)*np.cos((2.0*k0*pi*z0)/(z0 - zf))*np.cos((4.0*k0*pi*z0)/(-z0 + zf) - 2.0*pi*zc*(k + k0/(-z0 + zf))))/(2.0*sqrt(W**(-2))) + (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(z0 - zf))**2))/4)*sqrt(pi)*np.sin(2.0*pi*zc*(k + k0/(z0 - zf)))*np.sin((2.0*k0*pi*z0)/(z0 - zf)))/(2.0*sqrt(W**(-2))) - (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(-z0 + zf))**2))/4)*sqrt(pi)*np.sin((2.0*k0*pi*z0)/(z0 - zf))*np.sin((4.0*k0*pi*z0)/(-z0 + zf) - 2.0*pi*zc*(k + k0/(-z0 + zf))))/(2.0*sqrt(W**(-2)))  
	imag = (-(np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(z0 - zf))**2))/4)*sqrt(pi)*np.cos((2.0*k0*pi*z0)/(z0 - zf))*np.sin(2.0*pi*zc*(k + k0/(z0 - zf))))/(2.0*sqrt(W**(-2))) + (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(z0 - zf))**2))/4)*sqrt(pi)*np.cos(2.0*pi*zc*(k + k0/(z0 - zf)))*np.sin((2.0*k0*pi*z0)/(z0 - zf)))/(2.0*sqrt(W**(-2))) + (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(-z0 + zf))**2))/4)*sqrt(pi)*np.cos((4.0*k0*pi*z0)/(-z0 + zf) - 2.0*pi*zc*(k + k0/(-z0 + zf)))*np.sin((2.0*k0*pi*z0)/(z0 - zf)))/(2.0*sqrt(W**(-2))) + (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(-z0 + zf))**2))/4)*sqrt(pi)*np.cos((2.0*k0*pi*z0)/(z0 - zf))*np.sin((4.0*k0*pi*z0)/(-z0 + zf) - 2.0*pi*zc*(k + k0/(-z0 + zf))))/(2.0*sqrt(W**(-2))))
	#print("real=%e, imag=%e" % (real, imag))
	return np.sqrt(real**2+imag**2)





z = getZArray()
	
vals = getSoundWavePacketFunction(k0, zc, W)(z)
newvals = np.zeros(len(z)-1)
for i in range(len(z)-1):
	newvals[i] = 0.5 * (vals[i] + vals[i+1])

numPoints = len(z) - 1
dz = z[1] - z[0]
print("NUMPOINTS ")
print(numPoints)
print("DZ ")
print(dz)
print("len newvals = %d equal NUMPOINTS = %d" % (len(newvals), numPoints))
plt.grid(True)
#Y=fft(vals)/numPoints * (zf -z0)
Y=fft(newvals) / numPoints
F=fftfreq(numPoints, dz)
print("max freq")
print(F[np.argmax(abs(Y))])

#print("F=")
#print(F)
#plt.plot(2.0 * pi * F  ,abs(Y), markersize=3, linestyle="-", marker="o")
plt.xlim(-60, 60)
f1 = abs(Y)
f2 = ftAn(F)
print(np.max(f2) / np.max(f1))
coef = zf - z0
plt.plot(F  , f1 * coef, markersize=3, linestyle="None", marker="o")
#for f in newF:
#	print("f = %e, valf = %e" % (f, ftAn(f)))

plt.plot(F , f2 , markersize=3, linestyle="None", marker="o", color='r')
#plt.plot(F ,abs(Y), markersize=3, linestyle="-", marker="o")
#plt.plot(intlen * F,abs(wFFTAn(F)), markersize=3, linestyle="o", marker="o", color="r")
plt.draw()
plt.show()
