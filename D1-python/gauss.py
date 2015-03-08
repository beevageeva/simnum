import numpy as np
from constants import z0, zf, gamma
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from math import pi, sqrt, exp


from sound_wave_packet_params import  W, zc,k0, getSoundWaveFunction

from common import getZArray

z = getZArray()





#plt.plot(z, csder(z), markersize=3, linestyle="None", marker="o", color="r")
#plt.plot(z, cs(z), markersize=3, linestyle="None", marker="o", color="r")
#plt.plot(z, np.gradient(cs(z)), markersize=3, linestyle="None", marker="o", color="r")
plt.plot(z, getSoundWaveFunction(k0, zc, W)(z), markersize=3, linestyle="None", marker="o", color="r")



plt.draw()
plt.show()


