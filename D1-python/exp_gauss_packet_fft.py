import numpy as np
from constants import z0, zf
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from math import pi


from sound_wave_packet_params import k0,zc,W

print("sound wva epacket paranms z0=%1.3f, zf = %1.3f, k0 = %d, zc = %1.3f, W = %1.3f"  % (z0, zf, k0,zc,W))
print("sound wva epacket paranms %1.3f, %1.3f, %1.3f, %1.3f, %1.3f"  % (zc, z0, zf,k0,W))


def gaussPacketFunction(z):
  return np.exp(-(z - zc)**2 / W**2) *  np.cos( 2.0 * pi * k0 * (z-z0) / (zf - z0))


def fftMathematica(k):
#FourierTransform[
# Exp[-(z - zc)^2/W^2] Cos[2.0 Pi k0 (z - z0)/(zf - z0)], z, k, 
# FourierParameters -> {0, -2 Pi}]	
#
#
#
#1/(2 Sqrt[1/W^2]) E^(-(zc^2/W^2) - ((0. + 6.28319 I) k0 z0)/(
#  z0 - 1. zf)) (E^(
#   1/4 (2 I k \[Pi] W - (2 zc)/W - ((0. + 6.28319 I) k0 W)/(
#      z0 - 1. zf))^2) + E^(
#   1/4 (2 I k \[Pi] W - (2 zc)/W + ((0. + 6.28319 I) k0 W)/(
#       z0 - 1. zf))^2 + ((0. + 12.5664 I) k0 z0)/(
#    z0 - 1. zf))) Sqrt[\[Pi]]


#6.28319 = 2 pi
#12.5664 = 4 pi

	return (1.0/(2.0 * np.sqrt(1.0/W**2)) *  np.exp(-(zc**2/W**2) - ( 2.0 * pi * 1j*  k0 *z0)/(z0 - zf)) *
 ( np.exp(1.0/4.0 * (2.0 * 1j * k * pi * W - (2.0 * zc)/W - (2.0 * pi * 1j *  k0 * W)/(z0 -zf))**2) + 
	np.exp(1.0/4.0 * (2.0 * 1j * k * pi * W - (2.0 * zc)/W + (2.0 * pi * 1j * k0 * W)/(z0 - zf))**2 + (4.0 * pi * 1j * k0 * z0)/(z0 - zf))) * 
	np.sqrt(pi) )

#no parameters
#FourierTransform[Exp[-(z - zc)^2/W^2] Cos[2.0 Pi k0 (z - z0)/(zf - z0)], z, k]
#1/(2 Sqrt[2] Sqrt[1/W^2]) E^(-(zc^2/W^2) - ((0. + 6.28319 I) k0 z0)/(
#  z0 - 1. zf)) (E^(
#   1/4 (I k W + (2 zc)/W + ((0. + 6.28319 I) k0 W)/(z0 - 1. zf))^2) + 
#   E^(1/4 (I k W + (2 zc)/W - ((0. + 6.28319 I) k0 W)/(
#       z0 - 1. zf))^2 + ((0. + 12.5664 I) k0 z0)/(z0 - 1. zf)))

#	return 1.0/(2.0 * np.sqrt(2) * np.sqrt(1.0/W**2)) *  np.exp(-(zc**2/W**2) - (2.0 * pi * 1j * k0 * z0)/(z0 - zf)) *
# ( np.exp(1.0/4.0 (1j * k* W + (2.0 * zc)/W + (2.0 * pi * 1j * k0 * W)/(z0 - zf))**2) + 
#   np.exp(1.0/4.0 (1j * k * W + (2.0 * zc)/W - (2.0 * pi * 1j *  k0 * W)/(z0 - zf))**2 +(4.0 * pi * 1j * k0 * z0)/(z0 - zf)))


from common import getZArray

z = getZArray()
dz = z[1] - z[0]
numPoints = len(z)
from scipy.fftpack import fft,fftfreq#forFourierTransform




Y=fft(gaussPacketFunction(z))/numPoints
#Y=fft(self.pres)
#F=fftfreq(numPoints, dz)
F=fftfreq(numPoints)
#plt.plot(F,abs(Y), markersize=3, linestyle="None", marker="o", color="r")
#mantissa*^n

np.set_printoptions(threshold='nan')
y = gaussPacketFunction(z)
for yval in y:
	print("%1.5f "%yval ),

#print(gaussPacketFunction(z))
#print("F")
#print(F)
#print("coef")
#print(abs(fftMathematica(F)))

#plt.plot(F,abs(fftMathematica(F)), markersize=3, linestyle="None", marker="o", color="r")

#x = np.arange(0,10000)
#plt.plot(x, abs(fftMathematica(x)), markersize=3, linestyle="None", marker="o", color="r")

#print("max frequency %d" %  F[np.argmax(abs(Y))])

plt.draw()
plt.show()


