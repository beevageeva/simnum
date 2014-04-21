import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft,fftfreq#forFourierTransform
from common import getZArray
from sound_wave_defined_params import wFFTAn, w, kf, z0, zf



z = getZArray()

vals = w(z)

numPoints = len(z)
plt.grid(True)
#Y=fft(vals)/numPoints * (zf -z0)
Y=fft(vals) / numPoints
F=fftfreq(numPoints, z[0] - z[1])
plt.xlim(-80,80)
plt.plot(kf * F * 3 ,abs(Y), markersize=3, linestyle="-", marker="o")
#plt.plot(F ,abs(Y), markersize=3, linestyle="-", marker="o")
plt.plot(F,abs(wFFTAn(F)), markersize=3, linestyle="o", marker="o", color="r")
plt.draw()
plt.show()
