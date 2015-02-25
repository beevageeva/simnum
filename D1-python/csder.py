import numpy as np
from constants import z0, zf, gamma
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from math import pi, sqrt, exp


from sound_wave_params import rho00, rho01, ze, we, p00

print("sound wva epacket paranms z0=%1.3f, zf = %1.3f, rho00 = %1.3f, rho01=%1.3f,  ze = %1.3f, we = %1.3f"  % (z0, zf, rho00, rho01, ze, we))

densFunc = lambda z: rho00 + 0.5 * (rho01-rho00) * (1 + np.tanh((z-ze)/we))

def cs(z):
	return 1.0 / densFunc(z)


def  csder(z):
	return np.sqrt(gamma * p00) *  (-(-0.5*rho00 + 0.5*rho01)*(-np.tanh((z - ze)/we)**2 + 1)/(2*we*(rho00 + (-0.5*rho00 + 0.5*rho01)*(np.tanh((z - ze)/we) + 1))**(3/2)))



from common import getZArray

z = getZArray()





#plt.plot(z, csder(z), markersize=3, linestyle="None", marker="o", color="r")
#plt.plot(z, cs(z), markersize=3, linestyle="None", marker="o", color="r")
plt.plot(z, np.gradient(cs(z)), markersize=3, linestyle="None", marker="o", color="r")



plt.draw()
plt.show()


