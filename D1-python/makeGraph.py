
import numpy as np
from constants import z0, zf, gamma
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from math import pi, sqrt, exp


x = np.loadtxt("kcGraph")
y = np.loadtxt("cpGraph")
plt.plot(x, y)


plt.draw()
plt.show()
