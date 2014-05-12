import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft,fftfreq#forFourierTransform
from common import getZArray
from sound_wave_defined_params import wFFTAn, w, kf, z0, zf
from math import pi




z = getZArray()

vals = w(z)

numPoints = len(z)
dz = z[1] - z[0]
print("NUMPOINTS ")
print(numPoints)
print("DZ ")
print(dz)
plt.grid(True)
#Y=fft(vals)/numPoints * (zf -z0)
Y=fft(vals) / numPoints
F=fftfreq(numPoints, dz)
print("F=")
print(F)
intlen = z[z[len(z)-1]] - z[0]
#plt.plot(2.0 * pi * F  ,abs(Y), markersize=3, linestyle="-", marker="o")
plt.plot(F  ,abs(Y), markersize=3, linestyle="-", marker="o")
#plt.plot(F ,abs(Y), markersize=3, linestyle="-", marker="o")
#plt.plot(intlen * F,abs(wFFTAn(F)), markersize=3, linestyle="o", marker="o", color="r")
plt.draw()
plt.show()
