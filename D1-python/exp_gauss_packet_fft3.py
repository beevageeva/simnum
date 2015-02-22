import numpy as np
from constants import z0, zf
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from math import pi, sqrt, exp


from sound_wave_packet_params import k0,zc,W

print("sound wva epacket paranms z0=%1.3f, zf = %1.3f, k0 = %d, zc = %1.3f, W = %1.3f"  % (z0, zf, k0,zc,W))
print("sound wva epacket paranms %1.3f, %1.3f, %1.3f, %1.3f, %1.3f"  % (zc, z0, zf,k0,W))


def gaussPacketFunction(z):
  return np.exp(-(z - zc)**2 / W**2) *  np.cos( 2.0 * pi * k0 * (z-z0) / (zf - z0))


def fftMathNoAssumptions(k):
	from math import sqrt	
	real = (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(z0 - zf))**2))/4)*sqrt(pi)*np.cos(2.0*pi*zc*(k + k0/(z0 - zf)))*np.cos((2.0*k0*pi*z0)/(z0 - zf)))/(2.0*sqrt(W**(-2))) + (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(-z0 + zf))**2))/4)*sqrt(pi)*np.cos((2.0*k0*pi*z0)/(z0 - zf))*np.cos((4.0*k0*pi*z0)/(-z0 + zf) - 2.0*pi*zc*(k + k0/(-z0 + zf))))/(2.0*sqrt(W**(-2))) + (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(z0 - zf))**2))/4)*sqrt(pi)*np.sin(2.0*pi*zc*(k + k0/(z0 - zf)))*np.sin((2.0*k0*pi*z0)/(z0 - zf)))/(2.0*sqrt(W**(-2))) - (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(-z0 + zf))**2))/4)*sqrt(pi)*np.sin((2.0*k0*pi*z0)/(z0 - zf))*np.sin((4.0*k0*pi*z0)/(-z0 + zf) - 2.0*pi*zc*(k + k0/(-z0 + zf))))/(2.0*sqrt(W**(-2)))  
		imag = (-(np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(z0 - zf))**2))/4)*sqrt(pi)*np.cos((2.0*k0*pi*z0)/(z0 - zf))*np.sin(2.0*pi*zc*(k + k0/(z0 - zf))))/(2.0*sqrt(W**(-2))) + (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(z0 - zf))**2))/4)*sqrt(pi)*np.cos(2.0*pi*zc*(k + k0/(z0 - zf)))*np.sin((2.0*k0*pi*z0)/(z0 - zf)))/(2.0*sqrt(W**(-2))) + (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(-z0 + zf))**2))/4)*sqrt(pi)*np.cos((4.0*k0*pi*z0)/(-z0 + zf) - 2.0*pi*zc*(k + k0/(-z0 + zf)))*np.sin((2.0*k0*pi*z0)/(z0 - zf)))/(2.0*sqrt(W**(-2))) + (np.exp(-(zc**2/W**2) + (W**2.0*((4.0*zc**2)/W**4 - 4.0*pi**2.0*(k + k0/(-z0 + zf))**2))/4)*sqrt(pi)*np.cos((2.0*k0*pi*z0)/(z0 - zf))*np.sin((4.0*k0*pi*z0)/(-z0 + zf) - 2.0*pi*zc*(k + k0/(-z0 + zf))))/(2.0*sqrt(W**(-2))))
	return np.sqrt(real**2+imag**2)


def fftMathExplExp(m):
	#(Sqrt[Pi]*W)/E^((Pi*(k0^2*Pi*W^2 - 2*k0*(m*Pi*W^2 + I*(z0 - zc)*(z0 - zf)) + m*(m*Pi*W^2 - (2*I)*zc*(z0 - zf))))/(z0 - zf)^2)
	return (sqrt(pi)*W)/np.exp((pi*(k0**2*pi*W**2 - 2*k0*(m*pi*W**2 + 1j*(z0 - zc)*(z0 - zf)) + m*(m*pi*W**2 - (2*1j)*zc*(z0 - zf))))/(z0 - zf)**2)



def fftMathExplExpCE(m):
#	(Sqrt[Pi]*W*Cos[(Pi*(-2*k0*(z0 - zc)*(z0 - zf) - 2*m*zc*(z0 - zf)))/(z0 - zf)^2])/
#   E^((Pi*(k0^2*Pi*W^2 - 2*k0*m*Pi*W^2 + m^2*Pi*W^2))/(z0 - zf)^2) - 
#  (I*Sqrt[Pi]*W*Sin[(Pi*(-2*k0*(z0 - zc)*(z0 - zf) - 2*m*zc*(z0 - zf)))/(z0 - zf)^2])/
#   E^((Pi*(k0^2*Pi*W^2 - 2*k0*m*Pi*W^2 + m^2*Pi*W^2))/(z0 - zf)^2)

	real = (sqrt(pi)*W*np.cos((pi*(-2*k0*(z0 - zc)*(z0 - zf) - 2*m*zc*(z0 - zf)))/(z0 - zf)**2))/np.exp((pi*(k0**2*pi*W**2 - 2*k0*m*pi*W**2 + m**2*pi*W**2))/(z0 - zf)**2)  
	imag = (-sqrt(pi)*W*np.sin((pi*(-2*k0*(z0 - zc)*(z0 - zf) - 2*m*zc*(z0 - zf)))/(z0 - zf)**2))/np.exp((pi*(k0**2*pi*W**2 - 2*k0*m*pi*W**2 + m**2*pi*W**2))/(z0 - zf)**2)
	return np.sqrt(real**2 + imag**2)


def fftMathExpl(m):
#((E^(((4*I)*k0*Pi*z0*(zc + zf))/(z0 - zf)^2) + E^((4*k0*Pi*(m*Pi*W^2 + I*(z0^2 + zc*zf)))/(z0 - zf)^2))*Sqrt[Pi]*W)/
#  (2*E^((Pi*(k0^2*Pi*W^2 + m*(m*Pi*W^2 - (2*I)*zc*(z0 - zf)) + 2*k0*(m*Pi*W^2 + I*(z0 + zc)*(z0 + zf))))/(z0 - zf)^2))
	return ((np.exp(((4*1j)*k0*pi*z0*(zc + zf))/(z0 - zf)**2) + np.exp((4*k0*pi*(m*pi*W**2 + 1j*(z0**2 + zc*zf)))/(z0 - zf)**2))*sqrt(pi)*W)/(2*np.exp((pi*(k0**2*pi*W**2 + m*(m*pi*W**2 - (2*1j)*zc*(z0 - zf)) + 2*k0*(m*pi*W**2 + 1j*(z0 + zc)*(z0 + zf))))/(z0 - zf)**2))


def fftMath(k):
#((E^(-((Pi*(k0^2*Pi*W^2 - 2*k0*(k*Pi*W^2 - I*z0 + I*zc)*(z0 - zf) + k*(k*Pi*W^2 + (2*I)*zc)*(z0 - zf)^2))/(z0 - zf)^2)) + 
#      E^(-((Pi*(k0^2*Pi*W^2 + 2*k0*(k*Pi*W^2 - I*z0 + I*zc)*(z0 - zf) + k*(k*Pi*W^2 + (2*I)*zc)*(z0 - zf)^2))/(z0 - zf)^2)))*Sqrt[Pi]*W)/2
	return ((np.exp(-((pi*(k0**2*pi*W**2 - 2*k0*(k*pi*W**2 - 1j*z0 + 1j*zc)*(z0 - zf) + k*(k*pi*W**2 + (2*1j)*zc)*(z0 - zf)**2))/(z0 - zf)**2)) +  np.exp(-((pi*(k0**2*pi*W**2 + 2*k0*(k*pi*W**2 - 1j*z0 + 1j*zc)*(z0 - zf) + k*(k*pi*W**2 + (2*1j)*zc)*(z0 - zf)**2))/(z0 - zf)**2)))*sqrt(pi)*W)/2


def fftMathPart1(k):
#((E^(-((Pi*(k0^2*Pi*W^2 - 2*k0*(k*Pi*W^2 - I*z0 + I*zc)*(z0 - zf) + k*(k*Pi*W^2 + (2*I)*zc)*(z0 - zf)^2))/(z0 - zf)^2)) + 
#      E^(-((Pi*(k0^2*Pi*W^2 + 2*k0*(k*Pi*W^2 - I*z0 + I*zc)*(z0 - zf) + k*(k*Pi*W^2 + (2*I)*zc)*(z0 - zf)^2))/(z0 - zf)^2)))*Sqrt[Pi]*W)/2
	return ((np.exp(-((pi*(k0**2*pi*W**2 - 2*k0*(k*pi*W**2 - 1j*z0 + 1j*zc)*(z0 - zf) + k*(k*pi*W**2 + (2*1j)*zc)*(z0 - zf)**2))/(z0 - zf)**2)))*sqrt(pi)*W)/2

def fftMathPart2(k):
	return (( np.exp(-((pi*(k0**2*pi*W**2 + 2*k0*(k*pi*W**2 - 1j*z0 + 1j*zc)*(z0 - zf) + k*(k*pi*W**2 + (2*1j)*zc)*(z0 - zf)**2))/(z0 - zf)**2)))*sqrt(pi)*W)/2



from common import getZArray

z = getZArray()
dz = z[1] - z[0]
numPoints = len(z)
from scipy.fftpack import fft,fftfreq#forFourierTransform
Y=fft(gaussPacketFunction(z))/numPoints
F=fftfreq(numPoints, dz)
#plt.plot(F, abs(Y), markersize=3, linestyle="None", marker="o", color="r")




FM = np.arange(-80,80)
#plt.plot(FM, abs(fftMathExpl(FM)), markersize=3, linestyle="None", marker="o", color="r")
#plt.plot(FM, abs(fftMathExpl(FM)), markersize=3, linestyle="None", marker="o", color="g")

plt.plot(FM, abs(fftMath(FM/(zf-z0)) ), markersize=3, linestyle="None", marker="o", color="r")
plt.plot(FM, abs(fftMathPart1(FM/(zf-z0)) ), markersize=3, linestyle="None", marker="o", color="g")
plt.plot(FM, abs(fftMathPart2(FM/(zf-z0)) ), markersize=3, linestyle="None", marker="o", color="b")


#plot complex numbers!!!
#plt.plot(FM, (fftMath(FM/(zf-z0)) ), markersize=3, linestyle="None", marker="o", color="r")
#plt.plot(FM, (fftMathPart1(FM/(zf-z0)) ), markersize=3, linestyle="None", marker="o", color="g")
#plt.plot(FM, (fftMathPart2(FM/(zf-z0)) ), markersize=3, linestyle="None", marker="o", color="b")
#print("max frequency %d" %  F[np.argmax(abs(Y))])

plt.draw()
plt.show()


